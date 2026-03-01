import os
import sys

# Ensure RMG-Py is in path
rmg_path = os.path.expanduser("~/Developer/Maillard/RMG-Py")
sys.path.append(rmg_path)

from rmgpy.data.kinetics import KineticsDatabase
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.data.rmg import RMGDatabase

def test_maillard_families():
    print("Loading RMG Database...")
    db_path = os.path.expanduser("~/Developer/Maillard/RMG-database")
    
    database = RMGDatabase()
    
    print("Loading kinetics families...")
    try:
        # Load only the kinetics families, specifically targeting our custom ones
        database.kinetics = KineticsDatabase()
        database.kinetics.load_families(
            path=os.path.join(db_path, 'input', 'kinetics', 'families'),
            families=[
                "Schiff_Base_Formation",
                "Amadori_Rearrangement",
                "Heyns_Rearrangement",
                "Enolisation",
                "Retro_Aldol_Fragmentation",
                "Strecker_Degradation",
                "Thiol_Addition",
                "Aminoketone_Condensation",
                "Cysteine_Degradation",
                "Lipid_Schiff_Base",
                "Lipid_Thiazole_Condensation",
                "Beta_Elimination",
                "DHA_Crosslinking"
            ]
        )
        print("Successfully loaded custom Maillard kinetics families!")
        
        for fam_name, fam in database.kinetics.families.items():
            print(f" - Validated Family: {fam_name}")
            
    except Exception as e:
        print(f"FAILED to load kinetics families. Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_maillard_families()
