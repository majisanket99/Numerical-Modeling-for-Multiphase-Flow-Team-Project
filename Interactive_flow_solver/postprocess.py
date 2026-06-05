from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def load_slice_file(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    data = np.loadtxt(path)

    if data.ndim == 1:
        data = data.reshape(-1, 1)

    return data


def plot_single(data, title, label, output_path):
    plt.figure(figsize=(8, 5))

    if data.shape[1] >= 2:
        plt.plot(data[:, 0], data[:, 1], label=label, linewidth=2)
    else:
        plt.plot(data[:, 0], label=label, linewidth=2)

    plt.xlabel("Grid coordinate along slice")
    plt.ylabel("Normalized slice quantity")
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def plot_comparison(data1, label1, data2, label2, title, output_path):
    plt.figure(figsize=(8, 5))

    if data1.shape[1] >= 2:
        plt.plot(data1[:, 0], data1[:, 1], label=label1, linewidth=2)
    else:
        plt.plot(data1[:, 0], label=label1, linewidth=2)

    if data2.shape[1] >= 2:
        plt.plot(data2[:, 0], data2[:, 1], label=label2, linewidth=2)
    else:
        plt.plot(data2[:, 0], label=label2, linewidth=2)

    plt.xlabel("Grid coordinate along slice")
    plt.ylabel("Normalized slice quantity")
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def print_difference_metrics(reference, test, name):
    if reference.shape[1] >= 2 and test.shape[1] >= 2:
        diff = test[:, 1] - reference[:, 1]
    else:
        diff = test[:, 0] - reference[:, 0]

    mean_abs_difference = np.mean(np.abs(diff))
    max_abs_difference = np.max(np.abs(diff))

    print(f"{name}")
    print(f"  Mean absolute difference: {mean_abs_difference:.6f}")
    print(f"  Maximum absolute difference: {max_abs_difference:.6f}")
    print()


def main():
    results_dir = Path("results")
    figures_dir = Path("figures")
    figures_dir.mkdir(exist_ok=True)

    # Load all four slice files
    unp_1000 = load_slice_file(results_dir / "slice_unperturbed_nt1000.dat")
    per_1000 = load_slice_file(results_dir / "slice_perturbed_nt1000.dat")
    unp_2000 = load_slice_file(results_dir / "slice_unperturbed_nt2000.dat")
    per_2000 = load_slice_file(results_dir / "slice_perturbed_nt2000.dat")

    print("Loaded files successfully:")
    print("  unperturbed_nt1000:", unp_1000.shape)
    print("  perturbed_nt1000  :", per_1000.shape)
    print("  unperturbed_nt2000:", unp_2000.shape)
    print("  perturbed_nt2000  :", per_2000.shape)
    print()

    # Single plots
    plot_single(
        unp_1000,
        "Unperturbed Case at nt = 1000",
        "Unperturbed nt=1000",
        figures_dir / "unperturbed_nt1000.png",
    )

    plot_single(
        per_1000,
        "Perturbed Case at nt = 1000",
        "Perturbed nt=1000",
        figures_dir / "perturbed_nt1000.png",
    )

    plot_single(
        unp_2000,
        "Unperturbed Case at nt = 2000",
        "Unperturbed nt=2000",
        figures_dir / "unperturbed_nt2000.png",
    )

    plot_single(
        per_2000,
        "Perturbed Case at nt = 2000",
        "Perturbed nt=2000",
        figures_dir / "perturbed_nt2000.png",
    )

    # Comparison plots
    plot_comparison(
        unp_1000,
        "Unperturbed nt=1000",
        per_1000,
        "Perturbed nt=1000",
        "Flow Slice Comparison at nt = 1000",
        figures_dir / "comparison_nt1000.png",
    )

    plot_comparison(
        unp_2000,
        "Unperturbed nt=2000",
        per_2000,
        "Perturbed nt=2000",
        "Flow Slice Comparison at nt = 2000",
        figures_dir / "comparison_nt2000.png",
    )

    # Print metrics
    print_difference_metrics(unp_1000, per_1000, "Comparison at nt = 1000")
    print_difference_metrics(unp_2000, per_2000, "Comparison at nt = 2000")

    print("Saved figures:")
    print("  figures/unperturbed_nt1000.png")
    print("  figures/perturbed_nt1000.png")
    print("  figures/unperturbed_nt2000.png")
    print("  figures/perturbed_nt2000.png")
    print("  figures/comparison_nt1000.png")
    print("  figures/comparison_nt2000.png")


if __name__ == "__main__":
    main()
