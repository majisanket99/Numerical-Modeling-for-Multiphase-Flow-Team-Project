from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def load_slice(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return np.loadtxt(path)


def compute_metrics(unperturbed, perturbed):
    difference = perturbed[:, 1] - unperturbed[:, 1]
    mean_abs_difference = np.mean(np.abs(difference))
    max_abs_difference = np.max(np.abs(difference))
    return mean_abs_difference, max_abs_difference


results_dir = Path("results")
figures_dir = Path("figures")
figures_dir.mkdir(exist_ok=True)

time_steps = [1000, 2000, 3000]

mean_values = []
max_values = []

for nt in time_steps:
    unperturbed = load_slice(results_dir / f"slice_unperturbed_nt{nt}.dat")
    perturbed = load_slice(results_dir / f"slice_perturbed_nt{nt}.dat")

    mean_diff, max_diff = compute_metrics(unperturbed, perturbed)

    mean_values.append(mean_diff)
    max_values.append(max_diff)

    print(f"nt = {nt}")
    print(f"  Mean absolute difference: {mean_diff:.6f}")
    print(f"  Maximum absolute difference: {max_diff:.6f}")
    print()

labels = [f"nt = {nt}" for nt in time_steps]

plt.figure(figsize=(8, 5))

plt.plot(labels, mean_values, marker="o", linewidth=2, label="Mean absolute difference")
plt.plot(labels, max_values, marker="o", linewidth=2, label="Maximum absolute difference")

plt.xlabel("Simulation time step")
plt.ylabel("Difference magnitude")
plt.title("Perturbation Effect at Different Time Steps")
plt.legend()
plt.grid(True)
plt.tight_layout()

output_path = figures_dir / "difference_growth_metrics.png"
plt.savefig(output_path, dpi=300)
plt.show()

print(f"Saved figure to: {output_path}")
