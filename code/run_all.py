"""Run the full replication pipeline in order."""
import subprocess, sys, pathlib, time

HERE = pathlib.Path(__file__).parent
STEPS = [
    ("build.py",              "Assemble the quarterly dataset from data/raw"),
    ("gh_critical_values.py", "Simulate Gregory-Hansen critical values (slow)"),
    ("econ.py",               "Unit roots, bounds tests, placebo, power analysis"),
    ("scenario.py",           "Calibrated stress scenario"),
]

for script, label in STEPS:
    print(f"\n{'='*70}\n{label}\n{'='*70}")
    t = time.time()
    r = subprocess.run([sys.executable, str(HERE / script)])
    if r.returncode:
        sys.exit(f"FAILED: {script}")
    print(f"-- {script} completed in {time.time()-t:.0f}s")

print("\nAll steps completed. Results are in output/.")
