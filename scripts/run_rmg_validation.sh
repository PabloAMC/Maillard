#!/usr/bin/env bash
# run_rmg_validation.sh
# 
# Master script to execute the complete Phase 1B RMG-Py validation suite.
# This script:
#   1. Generates the input.py files for the three test cases
#   2. Runs RMG on each test case (if installed)
#   3. Validates the generated species_dictionary.txt against the target SMILES
#
# Usage: ./scripts/run_rmg_validation.sh

set -e

GREEN='\033[92m'
RED='\033[91m'
YELLOW='\033[93m'
NC='\033[0m'

echo -e "${GREEN}==> Generating RMG input files...${NC}"
python3 scripts/generate_rmg_inputs.py

if ! command -v rmg &> /dev/null; then
    echo -e "\n${YELLOW}WARNING: 'rmg' command not found in PATH.${NC}"
    echo "This indicates RMG-Py is not installed in the active conda environment."
    echo "The input scripts have been successfully generated in:"
    echo "  data/reactions/rmg_validation_cases/"
    echo ""
    echo "To run manually once RMG is installed:"
    echo "  rmg data/reactions/rmg_validation_cases/case_1_ribose_cys/input.py"
    echo "  python3 scripts/validate_rmg_output.py data/reactions/rmg_validation_cases/case_1_ribose_cys 'SCc1ccco1' 'Cc1occc1S'"
    exit 0
fi

echo -e "\n${GREEN}==> RMG installed. Running validations...${NC}"

# Case 1: Ribose + Cys -> FFT, MFT
echo -e "\n--- Case 1: Ribose + Cysteine ---"
rmg data/reactions/rmg_validation_cases/case_1_ribose_cys/input.py
python3 scripts/validate_rmg_output.py data/reactions/rmg_validation_cases/case_1_ribose_cys "SCc1ccco1" "Cc1occc1S"

# Case 2: Glucose + Glycine -> Furfural
echo -e "\n--- Case 2: Glucose + Glycine ---"
rmg data/reactions/rmg_validation_cases/case_2_glucose_gly/input.py
python3 scripts/validate_rmg_output.py data/reactions/rmg_validation_cases/case_2_glucose_gly "O=Cc1ccco1"

# Case 3: Ribose + Cys + Leu -> FFT, MFT, 3-methylbutanal, Pyrazine
echo -e "\n--- Case 3: Ribose + Cys + Leu ---"
rmg data/reactions/rmg_validation_cases/case_3_ribose_cys_leu/input.py
python3 scripts/validate_rmg_output.py data/reactions/rmg_validation_cases/case_3_ribose_cys_leu "SCc1ccco1" "Cc1occc1S" "CC(C)CC=O" "Cc1nccnc1C"

echo -e "\n${GREEN}==> All RMG validations completed successfully.${NC}"
