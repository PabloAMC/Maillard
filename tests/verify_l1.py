import sys
from pathlib import Path

# Add src to path
sys.path.append(str(Path.cwd()))

def test_imports():
    try:
        from src import MaillardPipeline, FormulationResult, ReactionConditions, SmirksEngine, FormulationOptimizer
        print("✅ PUBLIC API IMPORT SUCCESSFUL")
    except ImportError as e:
        print(f"❌ PUBLIC API IMPORT FAILED: {e}")
        sys.exit(1)

def test_sugar_classifier():
    from src.pathway_extractor import Species
    from src.sugar_classifier import is_sugar, is_hexose, is_pentose
    
    glucose = Species(label="Glucose", smiles="C(C1C(C(C(C(O1)O)O)O)O)O")
    ribose = Species(label="Ribose", smiles="C1C(C(C(C(O1)O)O)O)O")
    
    is_glucose_sugar = is_sugar(glucose)
    is_glucose_hexose = is_hexose(glucose)
    is_ribose_pentose = is_pentose(ribose)
    
    print(f"Glucose is sugar: {is_glucose_sugar}")
    print(f"Glucose is hexose: {is_glucose_hexose}")
    print(f"Ribose is pentose: {is_ribose_pentose}")
    
    if is_glucose_hexose and is_ribose_pentose:
        print("✅ SUGAR CLASSIFIER WORKING")
    else:
        print("❌ SUGAR CLASSIFIER FAILED")
        sys.exit(1)

if __name__ == "__main__":
    test_imports()
    test_sugar_classifier()
