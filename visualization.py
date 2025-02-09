import json
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats


def geometric_mean(data):
    """Calculates the geometric mean of a list of numbers, ignoring zeros."""
    data = [x for x in data if x > 0]  # Exclude zero values
    return round(scipy.stats.mstats.gmean(data)) if data else 0



def compute_sigma_range(data, gm):
    """Computes the 3-sigma range for a geometric mean, avoiding zero/negative values."""
    data = np.array(data)
    valid_data = data[data > 0]  # Remove zeros and negative values

    if len(valid_data) == 0:
        return gm, gm  # If no valid data, return the geometric mean itself

    log_data = np.log(valid_data)  # Log-transform only valid data
    mu = np.mean(log_data)   # Mean of log data
    sigma = np.std(log_data) # Standard deviation of log data

    lower = np.exp(mu - 3 * sigma)  # Convert back to original scale
    upper = np.exp(mu + 3 * sigma)

    return lower, upper



def visualize_histograms(json_file):
    """Reads JSON file and generates histograms for ICL1, ICL2, ICL3, and CTER lengths with bin size of 5."""
    # Load JSON data
    with open(json_file, "r") as file:
        data = json.load(file)

    # Extract length data for each region
    len_ICL1 = [x for x in [entry["len_ICL1"] for entry in data] if x <= 100]
    len_ICL2 = [x for x in [entry["len_ICL2"] for entry in data] if x <= 100]
    len_ICL3 = [x for x in [entry["len_ICL3"] for entry in data] if x <= 100]
    len_CTER = [x for x in [entry["len_CTER"] for entry in data] if x <= 250]

    # Calculate geometric means
    gm_ICL1 = geometric_mean(len_ICL1)
    gm_ICL2 = geometric_mean(len_ICL2)
    gm_ICL3 = geometric_mean(len_ICL3)
    gm_CTER = geometric_mean(len_CTER)

    # Compute sigma ranges
    sigma_range_ICL1 = compute_sigma_range(len_ICL1, gm_ICL1)
    sigma_range_ICL2 = compute_sigma_range(len_ICL2, gm_ICL2)
    sigma_range_ICL3 = compute_sigma_range(len_ICL3, gm_ICL3)
    sigma_range_CTER = compute_sigma_range(len_CTER, gm_CTER)

    # Create histograms for each region
    plt.figure(figsize=(12, 8))
    plt.suptitle("Inner Membrane Region Lengths of annotated human GPCR transcripts", fontsize=18, fontweight='bold')
#    plt.suptitle("Inner Membrane Region Lengths of all human GPCR transcripts", fontsize=18, fontweight='bold')


    for i, (lengths, gm, sigma_range, color, title, size) in enumerate([
        (len_ICL1, gm_ICL1, sigma_range_ICL1, 'b', "ICL1 Lengths", 14),
        (len_ICL2, gm_ICL2, sigma_range_ICL2, 'g', "ICL2 Lengths", 14),
        (len_ICL3, gm_ICL3, sigma_range_ICL3, 'r', "ICL3 Lengths", 14),
        (len_CTER, gm_CTER, sigma_range_CTER, 'purple', "C-ter Lengths", 14)
    ]):
        plt.subplot(2, 2, i + 1)
        if sigma_range[0] < sigma_range[1]:
            plt.axvspan(max(0, sigma_range[0]), sigma_range[1], color='grey', alpha=0.2)
        bin_sizes = [50, 50, 150, 250]  # Define bins for each subplot
        plt.hist(lengths, bins=range(0, bin_sizes[i] + 1, 1), alpha=0.7, color=color, edgecolor='black')
        plt.axvline(gm, color='black', linestyle='dashed', linewidth=1.5)
        plt.text(gm + 5, plt.ylim()[1] * 0.9, f'{gm}', color='black', fontsize=14)
        plt.xlabel("Length of loop", fontsize=10)
        plt.ylabel("Number of transcripts", fontsize=10)
        plt.title(title, fontsize=size)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


# Example usage
json_file = "loop_lengths_full.json"  # Update the path if needed
visualize_histograms(json_file)
