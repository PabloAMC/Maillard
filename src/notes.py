"""Note headings: which odour note a compound belongs to (STATUS_AND_PLAN sec. 0 step 3).

A heading groups the per-compound outputs the way a formulator reads them. It is a declared grouping of the
descriptors already on file (``data/species/desirable_targets.yml`` and ``off_flavour_targets.yml``,
``sensory_desc``), not a measurement: no summed note score is computed (PREREG-level rule of sec. 0 step 3), and
a compound with no heading is listed under "other". Thresholds there are marked uncited, so no odour-activity weight
is applied here.
"""
from __future__ import annotations

from typing import Optional

from src import compound_keys

#: heading -> registry ids, in reading order
HEADINGS = {
    "meaty (thiols, sulfides, thiazoles)": (
        "2_methyl_3_furanthiol", "2_furfurylthiol", "bis_2_methyl_3_furyl_disulfide", "bis_2_furfuryl_disulfide",
        "2_methyl_3_furfuryldithio_furan", "hydrogen_sulfide", "dimethyl_disulfide", "dimethyl_trisulfide",
        "2_methylthiophene", "4_5_dihydro_2_methylthiazole", "2_pentyl_4_methylthiazole", "2_hexyl_4_methylthiazole"),
    "roasted, nutty (pyrazines, acetylthiazole)": ("2_acetylthiazole", "2_3_dimethylpyrazine", "2_5_dimethylpyrazine", "2_6_dimethylpyrazine",
                                   "2_ethyl_3_5_dimethylpyrazine", "2_acetyl_1_pyrroline"),
    "malty, potato (Strecker aldehydes, methional)": ("methional", "2_methylbutanal", "3_methylbutanal",
                                                      "2_methylpropanal", "phenylacetaldehyde"),
    "caramel, sweet (furanones, furfural)": ("hemf", "hdmf", "furfural"),
    "off-notes (beany, green, cardboard)": ("hexanal", "nonanal", "1_octen_3_ol", "2_pentylfuran", "1_hexanol",
                                            "3_isobutyl_2_methoxypyrazine"),
    "not an odorant (process and safety markers)": ("acrylamide", "hmf", "cml", "cel", "lysinoalanine",
                                                    "lanthionine", "phip", "meiqx"),
}
ORDER = list(HEADINGS) + ["other"]
_BY_ID = {cid: h for h, ids in HEADINGS.items() for cid in ids}


def note_of(compound: str) -> str:
    key: Optional[compound_keys.CompoundKey] = compound_keys.resolve(str(compound))
    return _BY_ID.get(key.id, "other") if key else "other"


def order_key(compound: str) -> int:
    return ORDER.index(note_of(compound))
