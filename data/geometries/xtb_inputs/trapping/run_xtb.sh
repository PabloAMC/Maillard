#!/bin/bash
# Auto-generated script to run xTB path search
# Generates the transition state guess. 

REACTION_NAME="trapping"
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
