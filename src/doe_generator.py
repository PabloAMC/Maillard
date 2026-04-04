"""
doe_generator.py

S12.3: Model-Guided Active Learning (DoE Loop)
This module formalizes "Structural Gaps" into explicit Design of Experiments (DoE) workflows.
When the pipeline identifies unresolvable literature gaps in the process registry, 
it generates a formal lab protocol optimized for calibration gain.
"""

import json
import os
import logging

logger = logging.getLogger(__name__)

# Basic templated DoE protocols based on the gap type
DOE_TEMPLATES = {
    "missing_absolute_anchor": {
        "method": "SIDA GC-MS/O",
        "factors": ["Temperature (95C, 120C)", "Time (10m, 30m)"],
        "instructions": "Prepare standard aqueous target matrix. Add isotopically labeled internal standards (e.g. 13C-MFT, d3-methional) before extraction. Heat uniformly."
    },
    "blocking_benchmark_gap": {
        "method": "Multi-factorial Quantitative Headspace SIDA",
        "factors": ["Precursor Load (1x, 5x)", "Temperature (90C, 130C)", "Matrix (SPI, PPI)"],
        "instructions": "Standard PBMA formulation baseline. Use Safe+SPME extraction to capture both highly volatile (H2S) and semi-volatile (pyrazines) simultaneously with SIDA quantitation."
    },
    "missing_positive_flavor_anchor": {
        "method": "Targeted GC-MS (Furanone band)",
        "factors": ["Water Activity", "pH (5.5, 6.5)"],
        "instructions": "Focus resolution on HEMF and DMHF specifically. Ensure SPME fiber captures polar furanones effectively without early saturation by aldehydes."
    },
    "missing_kinetic_dataset": {
        "method": "Time-course LC-MS/MS",
        "factors": ["Time (0, 5, 10, 20, 60 min)", "Temp (80C, 100C, 120C)"],
        "instructions": "Sample dynamically during the heating phase to capture pseudo-first-order formation rates. Quench immediately in ice bath."
    },
    "missing_process_state_bundle": {
        "method": "Simultaneous Ellman/OPA and DSC",
        "factors": ["Heating time", "Moisture Content (Extrusion regime)"],
        "instructions": "Run Differential Scanning Calorimetry alongside Ellman's reagent for free thiols and OPA for primary amines on the exact same homogenized samples to eliminate lot-to-lot variance."
    }
}

def generate_active_learning_requests(gap_registry_path: str) -> dict:
    """Read the structural gaps and generate explicit DoE requests."""
    
    if not os.path.exists(gap_registry_path):
        logger.warning(f"Gap registry not found at {gap_registry_path}. Cannot generate DoE.")
        return {"active_learning_requests": []}
        
    with open(gap_registry_path, "r", encoding="utf-8") as f:
        gaps = json.load(f)
        
    requests = []
    
    for entry in gaps.get("entries", []):
        if entry.get("wet_lab_requirement") != "required":
            continue
            
        gap_type = entry.get("gap_type", "unknown")
        template = DOE_TEMPLATES.get(gap_type, {
            "method": "Standard GC-MS Profiling",
            "factors": ["Vary primary bottleneck"],
            "instructions": "Perform standard empirical screen to establish directional bounds."
        })
        
        request = {
            "request_id": f"doe_request_{entry['gap_id']}",
            "priority": "HIGH" if "blocking" in gap_type else "MEDIUM",
            "target_chemistry_family": entry.get("chemistry_family"),
            "observables": entry.get("observable_panel_tags", []),
            "justification": entry.get("why_literature_cannot_close_it"),
            "suggested_next_step": entry.get("cheapest_next_step"),
            "experimental_design": template
        }
        
        requests.append(request)
        
    return {"active_learning_requests": requests}


def export_doe_requests(gap_registry_path: str, output_path: str):
    """Generates and writes the DoE requests to a JSON artifact."""
    payload = generate_active_learning_requests(gap_registry_path)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        
    logger.info(f"Exported {len(payload['active_learning_requests'])} Active Learning DoE requests to {output_path}")
    return payload
