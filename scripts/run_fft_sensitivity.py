import os
import pandas as pd
import subprocess

# Precursors: Ribose(0.1), Cys(0.05), Leu(0.05)
precursors = "ribose:0.1,cysteine:0.05,leucine:0.05"
ph_values = [4.5, 6.0, 7.5, 9.0]
output_dir = "results/sensitivity"
os.makedirs(output_dir, exist_ok=True)

results = []

for ph in ph_values:
    out_prefix = f"ribose_cys_leu_ph{ph}"
    print(f"--- Running pH {ph} ---")
    cmd = [
        ".venv/bin/python", "scripts/run_cantera_kinetics.py",
        "--precursors", precursors,
        "--temp", "150",
        "--time", "3600",
        "--from-smirks",
        "--no-gating",
        "--ph", str(ph),
        "--track", "H2S,furfural,2-furfurylthiol,pyruvaldehyde",
        "--output", f"{output_dir}/{out_prefix}"
    ]
    # Set PYTHONPATH
    env = os.environ.copy()
    env["PYTHONPATH"] = env.get("PYTHONPATH", "") + ":" + os.getcwd()
    
    subprocess.run(cmd, env=env, check=True)
    
    # Read final yields from CSV
    csv_path = f"{output_dir}/{out_prefix}_results.csv"
    df = pd.read_csv(csv_path)
    final = df.iloc[-1]
    
    results.append({
        "pH": ph,
        "furfural": final["furfural"],
        "2-furfurylthiol": final["2-furfurylthiol"],
        "pyruvaldehyde": final["pyruvaldehyde"],
        "H2S": final["H2S"]
    })

summary_df = pd.DataFrame(results)
print("\n--- SENSITIVITY SUMMARY ---")
print(summary_df)
summary_df.to_csv(f"{output_dir}/summary.csv", index=False)
