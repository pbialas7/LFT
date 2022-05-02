import subprocess as proc
import time

L = 8

for m2 in [0.8, 0.9, 1.0, 1.1, 1.15, 1.2, 1.25,  1.3,  1.35, 1.4]:
    m2i = int(100 * m2)
    start = time.time()
    results = proc.run(
        ['./cmake-build-release/phi4', '--lambda', "1.0", '--kappa', "1.0", '--n-term', "1000", '--n-sweeps', "1000000",
         '--save-freq',
         "10", '--seed', "34", '--Lx', f"{L}", '--Ly', f"{L}", '--M2', f"{-m2}", "--output",
         f"data/k100_m{m2i:03d}_{L:02d}x{L:02d}.bin"])
    end = time.time()
    print(f"ellapsed time {end - start:.2f}s")
