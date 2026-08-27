import re
import json
from pathlib import Path

def main():
    root = Path(__file__).resolve().parent
    path = root / "data" / "Gemini_Deep_Research" / "Gemini_deep_research_I.md"
    content = path.read_text(encoding="utf-8")
    
    # We will search for code blocks or JSON blocks
    # Specifically matches anything enclosed in { and } containing "analyte"
    matches = re.findall(r"\{[^{}]*?\"analyte\"[^{}]*?\}", content, re.DOTALL)
    print(f"Total blocks found: {len(matches)}")
    for idx, m in enumerate(matches, 1):
        print(f"\n--- Block {idx} ---")
        # clean markdown escapes
        cleaned = m.replace("\\_", "_").replace("\\", "")
        print(cleaned)
        try:
            parsed = json.loads(cleaned)
            print("Successfully parsed JSON!")
        except Exception as e:
            print(f"Parse error: {e}")

if __name__ == "__main__":
    main()
