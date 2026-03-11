"""
src/results_db.py

SQLite-backed storage for Maillard reaction barriers and calculation provenance.
Provides a queryable alternative to ad-hoc JSON files.
"""

import sqlite3
import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

from src.barrier_constants import get_barrier  # noqa: E402

class ResultsDB:
    def get_best_barrier(self, reactants: List[str], products: List[str], family: str = "unknown") -> Tuple[float, str, float]:
        """
        Single source of truth for barrier lookups.
        1. Queries DB for exact match (using method priority).
        2. Falls back to heuristic constants if no DB entry exists.
        
        Returns (barrier_kcal, source_string, uncertainty_kcal)
        """
        res = self.find_barrier(reactants, products)
        if res:
            method = res['method'].lower()
            # Expert Uncertainty Mapping (R.11)
            method_unc_map = {
                "wb97m-v": 1.5,
                "r2scan-3c": 2.0,
                "mace-off": 2.5,
                "xtb": 3.0,
                "hf": 5.0
            }
            unc = method_unc_map.get(method, 3.5)
            return res["delta_g_kcal"], f"DB:{res['method']}", unc
        
        # Fallback to heuristic
        barrier_kcal, unc = get_barrier(family)
        return barrier_kcal, "Heuristic", unc

    def __init__(self, db_path: str = "results/maillard_results.db"):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_db()

    def _get_connection(self):
        return sqlite3.connect(self.db_path)

    def _init_db(self):
        """Build the relational schema."""
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
            # We store the sorted SMILES lists as JSON strings for identity checking
            cursor.execute("""
                CREATE TABLE IF NOT EXISTS reactions (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    family TEXT,
                    reactants_json TEXT NOT NULL,
                    products_json TEXT NOT NULL,
                    UNIQUE(reactants_json, products_json)
                )
            """)
            
            # Barriers table: specific calculation results
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
                    FOREIGN KEY(reaction_id) REFERENCES reactions(id)
                )
            """)
            conn.commit()

    def _get_or_create_reaction(self, reactants: List[str], products: List[str], family: str = "unknown") -> int:
        """Get reaction ID, creating it if it doesn't exist."""
        # Sort to ensure identity regardless of input order
        r_json = json.dumps(sorted(reactants))
        p_json = json.dumps(sorted(products))
        
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
            return cursor.lastrowid

    def add_barrier(self, reactants: List[str], products: List[str], delta_g_kcal: float,
                    method: str, family: str = "unknown", basis: Optional[str] = None,
                    solvation: Optional[str] = None, cpu_time: Optional[float] = None,
                    converged: bool = True):
        """Add a calculation result to the database."""
        reaction_id = self._get_or_create_reaction(reactants, products, family)
        
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO barriers (reaction_id, delta_g_kcal, method, basis, solvation, cpu_time_sec, converged, timestamp)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                reaction_id, delta_g_kcal, method, basis, solvation, 
                cpu_time, converged, time.strftime('%Y-%m-%d %H:%M:%S')
            ))
            conn.commit()

    def find_barrier(self, reactants: List[str], products: List[str], 
                     method_preference: List[str] = ["wB97M-V", "r2SCAN-3c", "mace-off", "xtb", "hf", "literature_heuristic"]) -> Optional[Dict[str, Any]]:
        """
        Find the best available barrier for a reaction based on method priority.
        """
        r_json = json.dumps(sorted(reactants))
        p_json = json.dumps(sorted(products))
        
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT b.delta_g_kcal, b.method, b.basis, b.solvation, b.timestamp
                FROM barriers b
                JOIN reactions r ON b.reaction_id = r.id
                WHERE r.reactants_json = ? AND r.products_json = ? AND b.converged = 1
            """, (r_json, p_json))
            
            rows = cursor.fetchall()
            if not rows:
                return None
            
            # results: list of dicts
            results = []
            for r in rows:
                results.append({
                    "delta_g_kcal": r[0],
                    "method": r[1],
                    "basis": r[2],
                    "solvation": r[3],
                    "timestamp": r[4]
                })
            
            # Sort by preference
            for pref in method_preference:
                for res in results:
                    if res["method"].lower() == pref.lower():
                        return res
            
            # Fallback to the first found if no preference matched
            return results[0]

    def list_all_barriers(self) -> List[Dict[str, Any]]:
        """Utility for summary reports."""
        with self._get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT r.family, r.reactants_json, r.products_json, b.delta_g_kcal, b.method, b.timestamp
                FROM barriers b
                JOIN reactions r ON b.reaction_id = r.id
            """)
            return [dict(zip(["family", "reactants", "products", "barrier", "method", "time"], row)) for row in cursor.fetchall()]
