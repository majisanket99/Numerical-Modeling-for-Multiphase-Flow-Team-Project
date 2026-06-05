from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def load_slice_file(path):
    """
    Load a slice.dat file from the Fortran solver.
    The file is expected to contain numerical columns.
    """
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    data = np.loadtxt(path)

    if data.ndim == 1:
        data = data.reshape(-1, 1)

    return data


def main():
    results_dir = Path("results")
    figures_dir = Path("figures")
    figures_dir.mkdir(exist_ok=True)

    unperturbed_path = results_dir / "slice_unperturbed.dat"
    perturbed_path = results_dir / "slice_perturbed.dat"

    unperturbed = load_slice_file(unperturbed_path)
    perturbed = load_slice_file(perturbed_path)

    print("Unperturbed shape:", unperturbed.shape)
    print("Perturbed shape:", perturbed.shape)

    plt.figure(figsize=(8, 5))

    # If there are at least two columns, plot column 1 against column 0.
    # Otherwise, plot the values against index.
    if unperturbed.shape[1] >= 2:
        plt.plot(
            unperturbed[:, 0],
            unperturbed[:, 1],
            label="Unperturbed",
            linewidth=2,
        )
    else:
        plt.plot(
            unperturbed[:, 0],
            label="Unperturbed",
            linewidth=2,
        )

    if perturbed.shape[1] >= 2:
        plt.plot(
            perturbed[:, 0],
            perturbed[:, 1],
            label="Perturbed",
            linewidth=2,
        )
    else:
        plt.plot(
            perturbed[:, 0],
            label="Perturbed",
            linewidth=2,
        )

    plt.xlabel("Grid coordinate along slice")
    plt.ylabel("Normalized slice quantity")
    plt.title("Flow Slice Comparison: Baseline vs Interactive Perturbation")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    output_path = figures_dir / "slice_comparison.png"
    plt.savefig(output_path, dpi=300)
    difference = perturbed[:, 1] - unperturbed[:, 1]
    mean_abs_difference = np.mean(np.abs(difference))
    max_abs_difference = np.max(np.abs(difference))

    print("Mean absolute difference:", mean_abs_difference)
    print("Maximum absolute difference:", max_abs_difference)
    plt.show()

    print(f"Saved figure to: {output_path}")


if __name__ == "__main__":
    main()
