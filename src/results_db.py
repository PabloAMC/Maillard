"""
src/results_db.py

SQLite-backed storage for Maillard reaction barriers and calculation provenance.
Provides a queryable alternative to ad-hoc JSON files.

Provenance policy (2026-08-27, audit remediation TIER 3 items 3.7/3.8/3.9)
------------------------------------------------------------------------
1. The database path is anchored to the repository root, so a lookup never
   depends on the process working directory.
2. Only *genuinely computed* rows (DFT / composite-DFT / MLP / semiempirical
   tiers) may override the curated FAST_BARRIERS registry in
   `src/barrier_constants.py`.  Rows whose method string is a non-computed
   provenance marker ("literature_heuristic" and friends) are book-keeping,
   not evidence: they must never silently displace an owner-approved curated
   constant.  (The shipped DB used to serve Thiol_Addition 15.0 kcal/mol
   against the curated 28.60 — a ~1e7x rate difference.)
3. Barrier rows carry a `thermo_version` tag.  Lookups are filtered to the
   current thermodynamic-treatment version, so rows produced before a thermo
   fix (e.g. the QRRHO free-rotor correction) cannot leak into published
   numbers just because nobody deleted them.
4. Method matching is by method *family* prefix, not string equality, so
   "wB97M-V/def2-tzvp" is recognised as wB97M-V (previously it fell through
   the preference list, letting xTB outrank DFT and handing DFT rows the
   larger default sigma).
"""

import json
import math
import sqlite3
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from rdkit import Chem

from src.barrier_constants import get_barrier  # noqa: E402

# --- Repository-anchored default location -----------------------------------
# Anchoring on __file__ (not the CWD) so `ResultsDB()` resolves to the same file
# whether it is constructed from the repo root, from scripts/, or from a test
# tmpdir.  The previous relative default silently created an empty DB next to
# whatever directory the process happened to start in.
REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB_PATH = REPO_ROOT / "results" / "maillard_results.db"

# --- Thermo-treatment versioning --------------------------------------------
# Bump CURRENT_THERMO_VERSION whenever a change to the thermochemistry pipeline
# (QRRHO treatment, standard state, frequency scaling, ...) makes previously
# stored ΔG‡ values non-comparable.  Rows written before the tag existed are
# backfilled to LEGACY_THERMO_VERSION by the schema migration and are therefore
# excluded from lookups without being destroyed.
CURRENT_THERMO_VERSION = "qrrho_fixed_2026_08"
LEGACY_THERMO_VERSION = "pre_qrrho_fix_2026_08"

# Temperature at which a stored ΔG‡ was evaluated.  298.15 K is the historic
# implicit assumption of every row written before the column existed.
DEFAULT_BARRIER_TEMPERATURE_K = 298.15

# Method families, best first.  Matching is by prefix/containment on this list,
# so "wB97M-V/def2-tzvp", "wB97M-V//r2SCAN" and "GFN2-xTB" all resolve.
_METHOD_FAMILY_PRIORITY: Tuple[str, ...] = ("wb97m-v", "r2scan-3c", "mace-off", "xtb", "hf")

# Expert uncertainty mapping (R.11), keyed by method *family*.
_METHOD_FAMILY_UNCERTAINTY_KCAL: Dict[str, float] = {
    "wb97m-v": 1.5,
    "r2scan-3c": 2.0,
    "mace-off": 2.5,
    "xtb": 3.0,
    "hf": 5.0,
}
_UNKNOWN_METHOD_UNCERTAINTY_KCAL = 3.5

# Provenance markers that identify a row as *not* the product of a calculation.
# Such rows are transcriptions of literature/heuristic values and must not
# outrank the curated registry.
_NON_COMPUTED_METHOD_TOKENS: Tuple[str, ...] = (
    "heuristic",
    "literature",
    "estimate",
    "estimated",
    "guess",
    "manual",
    "surrogate",
    "placeholder",
)


def canonicalize_smiles(smiles: str) -> str:
    """Return the RDKit-canonical SMILES, or the input unchanged if unparseable.

    Reaction identity in this DB is the JSON of the sorted SMILES lists, so a
    non-canonical spelling creates a row that can never be looked up again.
    Non-SMILES labels (used by some tests and by mechanism placeholders) are
    passed through untouched rather than rejected.
    """
    text = str(smiles)
    try:
        mol = Chem.MolFromSmiles(text)
    except Exception:
        return text
    if mol is None:
        return text
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        return text


def _canonical_side(values: Iterable[str]) -> List[str]:
    return sorted(canonicalize_smiles(value) for value in values)


def method_family(method: Optional[str]) -> Optional[str]:
    """Map a free-form method string onto a known method family, or None."""
    if method is None:
        return None
    norm = str(method).strip().lower().replace("_", "-")
    for family in _METHOD_FAMILY_PRIORITY:
        if norm.startswith(family):
            return family
    for family in _METHOD_FAMILY_PRIORITY:
        if family in norm:
            return family
    return None


def is_computed_method(method: Optional[str]) -> bool:
    """True only for rows produced by an actual calculation tier."""
    norm = "" if method is None else str(method).strip().lower()
    if any(token in norm for token in _NON_COMPUTED_METHOD_TOKENS):
        return False
    return method_family(norm) is not None


def method_uncertainty_kcal(method: Optional[str]) -> float:
    family = method_family(method)
    if family is None:
        return _UNKNOWN_METHOD_UNCERTAINTY_KCAL
    return _METHOD_FAMILY_UNCERTAINTY_KCAL.get(family, _UNKNOWN_METHOD_UNCERTAINTY_KCAL)


def _coerce_barrier_value(value: Any) -> Optional[float]:
    """Reject non-numeric / non-finite ΔG‡ values.

    SQLite is dynamically typed: a REAL column happily stores the string
    'n/a' or NULL.  A silently-coerced barrier is worse than a missing one.
    """
    if value is None or isinstance(value, bool):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


class ResultsDB:
    def get_best_barrier(
        self,
        reactants: List[str],
        products: List[str],
        family: str = "unknown",
    ) -> Tuple[float, str, float]:
        """
        Single source of truth for barrier lookups.

        1. Query the DB for an exact reaction match (current thermo version,
           converged rows only, best method family first).
        2. Accept that row ONLY if it is a genuinely computed result.  A
           non-computed provenance row ("literature_heuristic", ...) is not
           allowed to override the curated FAST_BARRIERS registry, which is the
           owner-approved source of truth for FAST-mode barriers.
        3. Otherwise fall back to the curated constant for the family.

        Returns (barrier_kcal, source_string, uncertainty_kcal)
        """
        res = self.find_barrier(reactants, products)
        if res is not None and res.get("is_computed"):
            return res["delta_g_kcal"], f"DB:{res['method']}", method_uncertainty_kcal(res["method"])

        # Curated registry (or its heuristic default) wins over anything that
        # was not actually computed.
        barrier_kcal, unc = get_barrier(family)
        return barrier_kcal, "Heuristic", unc

    def __init__(self, db_path: Optional[str] = None):
        self.db_path = Path(db_path) if db_path is not None else DEFAULT_DB_PATH
        if not self.db_path.is_absolute():
            # Explicit relative paths (tests, ad-hoc scripts) are resolved
            # against the repo root rather than the process CWD.
            self.db_path = (REPO_ROOT / self.db_path).resolve()
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_db()

    @contextmanager
    def _get_connection(self) -> Iterator[sqlite3.Connection]:
        """Commit-on-success, always-close connection scope."""
        conn = sqlite3.connect(self.db_path)
        try:
            yield conn
            conn.commit()
        finally:
            conn.close()

    def _init_db(self):
        """Build the relational schema (and migrate older files in place)."""
        with self._get_connection() as conn:
            cursor = conn.cursor()

            # Species table: unique SMILES
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS species (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    smiles TEXT UNIQUE NOT NULL,
                    inchi TEXT,
                    name TEXT
                )
            """)

            # Reactions table: relates a set of reactants to products
            # We store the sorted CANONICAL SMILES lists as JSON strings for
            # identity checking (see canonicalize_smiles).
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS reactions (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    family TEXT,
                    reactants_json TEXT NOT NULL,
                    products_json TEXT NOT NULL,
                    UNIQUE(reactants_json, products_json)
                )
            """)

            # Barriers table: specific calculation results.
            #
            # temperature_k: the temperature at which delta_g_kcal (a ΔG‡) was
            #   evaluated.  DOCUMENTATION ONLY for now — every current consumer
            #   feeds delta_g_kcal into an Arrhenius expression as if it were an
            #   Ea (i.e. treats it as temperature-independent).  That is a known
            #   approximation; rewiring the consumers is deliberately out of
            #   scope here, but the evaluation temperature is now recorded so
            #   the approximation is auditable instead of invisible.
            # thermo_version: thermodynamic-treatment lineage; lookups filter on
            #   it so pre-fix rows cannot resurface.
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS barriers (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    reaction_id INTEGER,
                    delta_g_kcal REAL,
                    method TEXT NOT NULL,
                    basis TEXT,
                    solvation TEXT,
                    cpu_time_sec REAL,
                    converged BOOLEAN,
                    timestamp DATETIME,
                    temperature_k REAL,
                    thermo_version TEXT,
                    FOREIGN KEY(reaction_id) REFERENCES reactions(id)
                )
            """)

            cursor.execute("""
                CREATE TABLE IF NOT EXISTS ml_adoption_decisions (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    candidate_id TEXT NOT NULL,
                    model_family TEXT NOT NULL,
                    model_name TEXT NOT NULL,
                    proposed_role TEXT NOT NULL,
                    decision TEXT NOT NULL,
                    benchmark_set_id TEXT NOT NULL,
                    coverage_ratio REAL,
                    rank_correlation REAL,
                    mean_abs_error_kcal REAL,
                    max_abs_error_kcal REAL,
                    stop_reasons_json TEXT,
                    rationale TEXT,
                    fallback_comparator TEXT,
                    benchmark_visible_gap TEXT,
                    approved_for_default BOOLEAN,
                    timestamp DATETIME
                )
            """)

            self._migrate_barrier_columns(cursor)

    @staticmethod
    def _migrate_barrier_columns(cursor: sqlite3.Cursor) -> None:
        """Additive schema migration for pre-existing database files.

        Existing rows predate both the thermo-version tag and the recorded
        evaluation temperature, so they are tagged as legacy (which excludes
        them from lookups) rather than assumed current.
        """
        columns = {row[1] for row in cursor.execute("PRAGMA table_info(barriers)")}

        if "temperature_k" not in columns:
            cursor.execute("ALTER TABLE barriers ADD COLUMN temperature_k REAL")
        if "thermo_version" not in columns:
            cursor.execute("ALTER TABLE barriers ADD COLUMN thermo_version TEXT")

        # Read first, write only if there is anything to backfill: this method
        # runs on every ResultsDB construction, and an unconditional UPDATE
        # would take a write lock each time.
        pending = cursor.execute(
            "SELECT count(*) FROM barriers WHERE temperature_k IS NULL OR thermo_version IS NULL"
        ).fetchone()[0]
        if pending:
            cursor.execute(
                "UPDATE barriers SET temperature_k = ? WHERE temperature_k IS NULL",
                (DEFAULT_BARRIER_TEMPERATURE_K,),
            )
            cursor.execute(
                "UPDATE barriers SET thermo_version = ? WHERE thermo_version IS NULL",
                (LEGACY_THERMO_VERSION,),
            )

    def _get_or_create_reaction(self, reactants: List[str], products: List[str], family: str = "unknown") -> int:
        """Get reaction ID, creating it if it doesn't exist."""
        # Canonicalize + sort so identity is independent of input order and of
        # SMILES spelling.
        r_json = json.dumps(_canonical_side(reactants))
        p_json = json.dumps(_canonical_side(products))

        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                "SELECT id FROM reactions WHERE reactants_json = ? AND products_json = ?",
                (r_json, p_json)
            )
            row = cursor.fetchone()
            if row:
                return row[0]

            cursor.execute(
                "INSERT INTO reactions (family, reactants_json, products_json) VALUES (?, ?, ?)",
                (family, r_json, p_json)
            )
            res = cursor.lastrowid
            assert res is not None
            return int(res)

    def add_barrier(self, reactants: List[str], products: List[str], delta_g_kcal: float,
                    method: str, family: str = "unknown", basis: Optional[str] = None,
                    solvation: Optional[str] = None, cpu_time: Optional[float] = None,
                    converged: bool = True, temperature_k: float = DEFAULT_BARRIER_TEMPERATURE_K,
                    thermo_version: str = CURRENT_THERMO_VERSION):
        """Add a calculation result to the database."""
        value = _coerce_barrier_value(delta_g_kcal)
        if value is None:
            raise ValueError(f"delta_g_kcal must be a finite number, got {delta_g_kcal!r}")

        reaction_id = self._get_or_create_reaction(reactants, products, family)

        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO barriers (
                    reaction_id, delta_g_kcal, method, basis, solvation,
                    cpu_time_sec, converged, timestamp, temperature_k, thermo_version
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                reaction_id, value, method, basis, solvation,
                cpu_time, converged, time.strftime('%Y-%m-%d %H:%M:%S'),
                float(temperature_k), str(thermo_version),
            ))

    def find_barrier(
        self,
        reactants: List[str],
        products: List[str],
        method_preference: Optional[Sequence[str]] = None,
        thermo_version: Optional[str] = CURRENT_THERMO_VERSION,
    ) -> Optional[Dict[str, Any]]:
        """
        Find the best available barrier for a reaction based on method priority.

        `method_preference` entries are resolved to method families, so both
        "wB97M-V" and "wB97M-V/def2-tzvp" select the same rows.
        `thermo_version=None` disables the version filter (inspection only).
        """
        preference = _METHOD_FAMILY_PRIORITY if method_preference is None else tuple(method_preference)
        preferred_families: List[str] = []
        for entry in preference:
            family = method_family(entry)
            if family is not None and family not in preferred_families:
                preferred_families.append(family)

        r_json = json.dumps(_canonical_side(reactants))
        p_json = json.dumps(_canonical_side(products))

        query = """
            SELECT b.delta_g_kcal, b.method, b.basis, b.solvation, b.timestamp,
                   b.temperature_k, b.thermo_version, b.id
            FROM barriers b
            JOIN reactions r ON b.reaction_id = r.id
            WHERE r.reactants_json = ? AND r.products_json = ? AND b.converged = 1
        """
        params: List[Any] = [r_json, p_json]
        if thermo_version is not None:
            query += " AND b.thermo_version = ?"
            params.append(str(thermo_version))
        # Deterministic ordering: newest converged row first, so insertion order
        # never decides which of several equally-preferred rows is served.
        query += " ORDER BY b.timestamp DESC, b.id DESC"

        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query, params)
            rows = cursor.fetchall()

        results: List[Dict[str, Any]] = []
        for row in rows:
            value = _coerce_barrier_value(row[0])
            if value is None:
                # Corrupt / non-numeric row: skip rather than crash a caller
                # that expects a float.
                continue
            results.append({
                "delta_g_kcal": value,
                "method": row[1],
                "basis": row[2],
                "solvation": row[3],
                "timestamp": row[4],
                "temperature_k": row[5],
                "thermo_version": row[6],
                "method_family": method_family(row[1]),
                "is_computed": is_computed_method(row[1]),
            })

        if not results:
            return None

        for family in preferred_families:
            for res in results:
                if res["method_family"] == family:
                    return res

        # No preferred family matched: return the newest row, still labelled
        # with its provenance so the caller can refuse it.
        return results[0]

    def list_all_barriers(self, thermo_version: Optional[str] = None) -> List[Dict[str, Any]]:
        """Utility for summary reports (unfiltered by default)."""
        query = """
            SELECT r.family, r.reactants_json, r.products_json, b.delta_g_kcal,
                   b.method, b.timestamp, b.temperature_k, b.thermo_version
            FROM barriers b
            JOIN reactions r ON b.reaction_id = r.id
        """
        params: List[Any] = []
        if thermo_version is not None:
            query += " WHERE b.thermo_version = ?"
            params.append(str(thermo_version))
        query += " ORDER BY b.id ASC"

        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query, params)
            rows = cursor.fetchall()

        keys = ["family", "reactants", "products", "barrier", "method", "time",
                "temperature_k", "thermo_version"]
        return [dict(zip(keys, row)) for row in rows]

    def add_ml_adoption_decision(self, decision: Dict[str, Any]):
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                """
                INSERT INTO ml_adoption_decisions (
                    candidate_id, model_family, model_name, proposed_role, decision,
                    benchmark_set_id, coverage_ratio, rank_correlation,
                    mean_abs_error_kcal, max_abs_error_kcal, stop_reasons_json,
                    rationale, fallback_comparator, benchmark_visible_gap,
                    approved_for_default, timestamp
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    str(decision.get("candidate_id", "unknown")),
                    str(decision.get("model_family", "unknown")),
                    str(decision.get("model_name", "unknown")),
                    str(decision.get("proposed_role", "unknown")),
                    str(decision.get("decision", "defer")),
                    str(decision.get("benchmark_set_id", "unknown")),
                    decision.get("coverage_ratio"),
                    decision.get("rank_correlation"),
                    decision.get("mean_abs_error_kcal"),
                    decision.get("max_abs_error_kcal"),
                    json.dumps(decision.get("stop_reasons", [])),
                    str(decision.get("rationale", "")),
                    str(decision.get("fallback_comparator", "unknown")),
                    str(decision.get("benchmark_visible_gap", "unknown")),
                    bool(decision.get("approved_for_default", False)),
                    time.strftime('%Y-%m-%d %H:%M:%S'),
                ),
            )

    def list_ml_adoption_decisions(self, candidate_id: Optional[str] = None) -> List[Dict[str, Any]]:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            if candidate_id:
                cursor.execute(
                    """
                    SELECT candidate_id, model_family, model_name, proposed_role, decision,
                           benchmark_set_id, coverage_ratio, rank_correlation,
                           mean_abs_error_kcal, max_abs_error_kcal, stop_reasons_json,
                           rationale, fallback_comparator, benchmark_visible_gap,
                           approved_for_default, timestamp
                    FROM ml_adoption_decisions
                    WHERE candidate_id = ?
                    ORDER BY id ASC
                    """,
                    (candidate_id,),
                )
            else:
                cursor.execute(
                    """
                    SELECT candidate_id, model_family, model_name, proposed_role, decision,
                           benchmark_set_id, coverage_ratio, rank_correlation,
                           mean_abs_error_kcal, max_abs_error_kcal, stop_reasons_json,
                           rationale, fallback_comparator, benchmark_visible_gap,
                           approved_for_default, timestamp
                    FROM ml_adoption_decisions
                    ORDER BY id ASC
                    """
                )
            rows = cursor.fetchall()

        payload: List[Dict[str, Any]] = []
        for row in rows:
            payload.append(
                {
                    "candidate_id": row[0],
                    "model_family": row[1],
                    "model_name": row[2],
                    "proposed_role": row[3],
                    "decision": row[4],
                    "benchmark_set_id": row[5],
                    "coverage_ratio": row[6],
                    "rank_correlation": row[7],
                    "mean_abs_error_kcal": row[8],
                    "max_abs_error_kcal": row[9],
                    "stop_reasons": json.loads(row[10] or "[]"),
                    "rationale": row[11],
                    "fallback_comparator": row[12],
                    "benchmark_visible_gap": row[13],
                    "approved_for_default": bool(row[14]),
                    "timestamp": row[15],
                }
            )
        return payload
