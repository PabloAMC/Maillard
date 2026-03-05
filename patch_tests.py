import re

with open('tests/test_smirks_engine.py', 'r') as f:
    code = f.read()

# Add test_mass_conservation to all test classes
test_classes = re.findall(r'class Test\w+:', code)

for cls in test_classes:
    if 'TestOutputCompatibility' in cls: continue
    
    insertion = """
    def test_mass_conservation(self):
        for step in self.steps:
            assert_balanced(step)"""
            
    # Find the class and append the new test method
    # It's easier to just find the end of the class or add it right after setup_method
    
    pattern = rf'({cls}.*?def setup_method\(self\):.*?self\.steps = .*?\n)'
    replacement = rf'\1{insertion}\n'
    code = re.sub(pattern, replacement, code, flags=re.DOTALL)

with open('tests/test_smirks_engine.py', 'w') as f:
    f.write(code)
