#!/bin/bash
echo "Running xTB path search for explicit water cluster"
xtb reactant.xyz --path product.xyz --gfn 2 --alpb water
