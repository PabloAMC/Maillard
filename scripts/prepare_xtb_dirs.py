import os
from pathlib import Path

# The remaining reactions to prepare for xTB path searches
REMAINING_REACTIONS = [
    "enolisation", 
    "strecker", 
    "cys_ribose", 
    "retro_aldol", 
    "dha", 
    "trapping", 
    "pyrazine"
]

def main():
    base_dir = Path("data/geometries/xtb_inputs")
    base_dir.mkdir(parents=True, exist_ok=True)
    
    script_content = """#!/bin/bash
# Auto-generated script to run xTB path search
# Generates the transition state guess. 

REACTION_NAME="{name}"
echo "Running xTB path search for $REACTION_NAME"

if [[ -f "reactant.xyz" && -f "product.xyz" ]]; then
    echo "Found atom-mapped geometries. Executing xTB path search..."
    # GFN2-xTB with implicit water (ALPB)
    xtb reactant.xyz --path product.xyz --gfn 2 --alpb water
else
    echo "ERROR: reactant.xyz and/or product.xyz missing in $(pwd)"
    echo "Did you run 'python scripts/generate_mapped_geometries.py'?"
    exit 1
fi
"""

    for reaction in REMAINING_REACTIONS:
        r_dir = base_dir / reaction
        r_dir.mkdir(parents=True, exist_ok=True)
        
        run_script = r_dir / "run_xtb.sh"
        with open(run_script, "w") as f:
            f.write(script_content.format(name=reaction))
            
        os.chmod(run_script, 0o755)
        
    print(f"Generated xTB setup directories in {base_dir}")

if __name__ == "__main__":
    main()
