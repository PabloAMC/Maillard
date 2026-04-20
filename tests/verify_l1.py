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
    from src.sugar_classifier import _is_sugar, _is_hexose, _is_pentose
    
    glucose = Species(label="Glucose", smiles="C(C1C(C(C(C(O1)O)O)O)O)O")
    ribose = Species(label="Ribose", smiles="C1C(C(C(C(O1)O)O)O)O")
    
    is_glucose_sugar = _is_sugar(glucose)
    is_glucose_hexose = _is_hexose(glucose)
    is_ribose_pentose = _is_pentose(ribose)
    
    print(f"Glucose is sugar: {is_glucose_sugar}")
    print(f"Glucose is hexose: {is_glucose_hexose}")
    print(f"Ribose is pentose: {is_ribose_pentose}")
    
    if is_glucose_hexose and is_ribose_pentose:
        print("✅ SUGAR CLASSIFIER WORKING")
    else:
        print("❌ SUGAR CLASSIFIER FAILED")
        sys.exit(1)

def test_archive():
    archive_dir = Path("scripts/dev_archive")
    files = ["debug.py", "traceback.txt", "get_internal_bench_results.py"]
    missing = [f for f in files if not (archive_dir / f).exists()]
    
    if not missing:
        print("✅ ARCHIVE SUCCESSFUL")
    else:
        print(f"❌ ARCHIVE ERRORED - MISSING: {missing}")
        sys.exit(1)

if __name__ == "__main__":
    test_imports()
    test_sugar_classifier()
    test_archive()
