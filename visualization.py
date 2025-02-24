import json
import re

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
from sklearn.mixture import GaussianMixture

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



def visualize_histograms(json_file, which):
    """Reads JSON file and generates histograms for ICL1, ICL2, ICL3, and CTER lengths with bin size of 5."""
    # Load JSON data
    with open(json_file, "r") as file:
        data = json.load(file)

    # Extract length data for each region
    len_ICL1 = [x for x in [entry["len_ICL1"] for entry in data] if x <= 100]
    len_ICL2 = [x for x in [entry["len_ICL2"] for entry in data] if x <= 100]
    len_ICL3 = [x for x in [entry["len_ICL3"] for entry in data] if x <= 100]
    len_CTER = [x for x in [entry["len_CTER"] for entry in data] if x <= 250]

    # Extract receptor name for each region
    for entry in data:
        name = re.search(r"_(\w+)_HUMAN", entry["header"]).group(1) if re.search(r"_(\w+)_HUMAN", entry["header"]) else entry["header"]

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
    plt.suptitle("Inner membrane region lengths of "+ which, fontsize=18, fontweight='bold')
#    plt.suptitle("Inner Membrane Region Lengths of all human GPCR transcripts", fontsize=18, fontweight='bold')

    for i, (lengths, gm, sigma_range, color, title, bin_max) in enumerate([
        (len_ICL1, gm_ICL1, sigma_range_ICL1, 'b', "ICL1 length", 100),
        (len_ICL2, gm_ICL2, sigma_range_ICL2, 'g', "ICL2 length", 100),
        (len_ICL3, gm_ICL3, sigma_range_ICL3, 'r', "ICL3 length", 100),
        (len_CTER, gm_CTER, sigma_range_CTER, 'purple', "C-ter length", 250)
    ]):
        plt.subplot(2, 2, i + 1)

        # Apply sigma range shading
        if sigma_range[0] < sigma_range[1]:
            plt.axvspan(max(0, sigma_range[0]), sigma_range[1], color='grey', alpha=0.2)

        plt.xlim(0, 140)
        """
        # Compute x-axis limits dynamically
        if lengths:
            min_val, max_val = min(lengths), max(lengths)
            avg_val = np.mean(lengths)
            x_min = max(0, min_val - (0.1 * avg_val))  # Ensure x_min is non-negative
            x_max = max_val + (0.1 * avg_val)
            plt.xlim(x_min, x_max)
        """

        # Create histogram and capture returned bars (patches)
        counts, bins, bars = plt.hist(lengths, bins=range(0, bin_max + 1, 1),
                                      color=color,  alpha=0.7, linewidth=0.2)

        # Draw a vertical dashed line at the geometric mean
        plt.axvline(gm, color='black', linestyle='dashed', linewidth=1.5)
        plt.text(gm + 5, plt.ylim()[1] * 0.9, f'{gm}', color='black', fontsize=12)

        """
        # Example: Annotate each bar with its index (or use your own list of names if applicable)
        for idx, bar in enumerate(bars):
            if bar.get_height() > 1:
                plt.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5,
                    str(idx),  # Replace this with a relevant label if needed
                    ha="center", va="bottom", fontsize="10", rotation=45
                )
        """

        # Labels
        plt.xlabel("Length of loop", fontsize=10)
        plt.ylabel("Number of transcripts", fontsize=10)
        plt.title(title, fontsize=14)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plotname = which.replace(" ", "_")
    plt.savefig(f"plots/ICL_histograms_{plotname}.png", dpi=300, bbox_inches="tight")  # Save the figure as a PNG file
    plt.show()


def visualize_ratio_histogram(json_file):
    """Plots the histogram of 'ratio_inner_to_total' with reduced bin size and fits a Gaussian Mixture.
    The overall GMM density and individual Gaussian components are overlayed on the histogram.
    """
    # Load JSON data
    with open(json_file, "r") as file:
        data = json.load(file)

    # Extract ratio values, ensuring the key exists
    ratios = [entry["ratio_inner_to_total"] for entry in data if "ratio_inner_to_total" in entry]

    # Convert to numpy array and reshape for the GMM (expects 2D array)
    ratios = np.array(ratios).reshape(-1, 1)

    bins = 100
    n_components = 4

    # Plot the histogram with a reduced bin size and normalize to density
    plt.figure(figsize=(8, 6))
    counts, bin_edges, _ = plt.hist(ratios.flatten(), bins=bins, density=True,
                                    color="lightgray", edgecolor="black", alpha=0.7, linewidth=0.2)
    plt.xlabel("Loop length ratio (inner/total)", fontsize=12)
    plt.ylabel("Number of GPCRs", fontsize=12)
    plt.title("Distribution of intracellular loop length ratios", fontsize=14)

    # Fit a Gaussian Mixture Model to the ratio data
    gmm = GaussianMixture(n_components=n_components, random_state=0)
    gmm.fit(ratios)

    # Generate a fine x-axis grid spanning the range of the data for plotting densities
    x = np.linspace(ratios.min(), ratios.max(), 1000).reshape(-1, 1)
    logprob = gmm.score_samples(x)
    pdf = np.exp(logprob)

    # Plot the overall GMM density on top of the histogram
    plt.plot(x, pdf, '-k', label="Overall GMM")

    # Plot individual Gaussian components
    def gaussian(x, mu, sigma):
        return (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

    for i in range(n_components):
        weight = gmm.weights_[i]
        mean = gmm.means_[i, 0]
        sigma = np.sqrt(gmm.covariances_[i][0][0])
        component_pdf = weight * gaussian(x, mean, sigma)
        # Label using the peak maximum value (i.e. the mean)
        plt.plot(x, component_pdf, '--', label=f"Peak: {mean:.2f}")

    plt.legend()
    plt.tight_layout()
    plt.savefig(f"plots/loop_ratio_histogram.png", dpi=300, bbox_inches="tight")  # Save the figure as a PNG file
    plt.show()

# Source data and define receptor families
json_file = "loop_lengths_full.json"  # Update the path if needed
uniprot_codes = {"human adrenoceptors": [
    "P35348",  # α1A (ADRA1A)
    "P35368",  # α1B (ADRA1B)
    "P25100",  # α1D (ADRA1D)
    "P08913",  # α2A (ADRA2A)
    "P18089",  # α2B (ADRA2B)
    "P18825",  # α2C (ADRA2C)
    "P08588",  # β1 (ADRB1)
    "P07550",  # β2 (ADRB2)
    "P13945"   # β3 (ADRB3)
    ],
    "human opsins": [
        "P08100",  # Rhodopsin (RHO),
        "Q9H1Y3",  # Opsin-3 (OPN3), encephalopsin, panopsin
        "Q9UHM6",  # Melanopsin (OPN4), melanopsin
        "Q6U736",  # Opsin-5 (OPN5) , GPR136, neuropsin, TP13
        "P03999",  # Short-wave-sensitive opsin 1 (OPN1SW) - blue
        "P04001",  # Medium-wave-sensitive opsin 1 (OPN1MW) - green
        "P04000",  # Long-wave-sensitive opsin 1 (OPN1LW) - red
        "O14718"   # RRH visual pigment-like receptor protein
    ],
    "human cannabinoid receptors": [
        "P21554",  # Cannabinoid receptor 1 (CB1)
        "P34972",  # Cannabinoid receptor 2 (CB2)
        "Q9Y2T6"   # GPR55 (proposed CB3)
    ],
    "human prostaglandin receptors": [
        "P34995",  # Prostaglandin E2 receptor EP1 subtype (PTGER1)
        "P43116",  # Prostaglandin E2 receptor EP2 subtype (PTGER2)
        "P43115",  # Prostaglandin E2 receptor EP3 subtype (PTGER3)
        "P35408",  # Prostaglandin E2 receptor EP4 subtype (PTGER4)
        "Q13258",  # Prostaglandin D2 receptor (PTGDR)
        "Q9Y5Y4",  # Prostaglandin D2 receptor 2 (PTGDR2)
        "P43088",  # Prostaglandin F2-alpha receptor (PTGFR)
        "P43119",  # Prostacyclin receptor (PTGIR)
        "Q9P2B2"   # Prostaglandin F2 receptor negative regulator (PTGFRN)
    ],
    "human histamine receptors": [
        "P35367",  # Histamine H1 receptor (HRH1)
        "P25021",  # Histamine H2 receptor (HRH2)
        "Q9Y5N1",  # Histamine H3 receptor (HRH3)
        "Q9H3N8"   # Histamine H4 receptor (HRH4)
    ],
    "human serotonin receptors": [
        "P08908",  # 5-HT1A receptor (HTR1A)
        "P28221",  # 5-HT1B receptor (HTR1B)
        "P28566",  # 5-HT1E receptor (HTR1E)
        "P30939",  # 5-HT1F receptor (HTR1F)
        "P28223",  # 5-HT2A receptor (HTR2A)
        "P41595",  # 5-HT2B receptor (HTR2B)
        "P28222",  # 5-HT2C receptor (HTR2C)
        "P46098",  # 5-HT3A receptor (HTR3A)
        "Q9Y5N1",  # 5-HT3B receptor (HTR3B)
        "Q13639",  # 5-HT4 receptor (HTR4)
        "P47898",  # 5-HT5A receptor (HTR5A)
        "P50406",  # 5-HT6 receptor (HTR6)
        "P34969"   # 5-HT7 receptor (HTR7)
    ],
    "human muscarinic receptors": [
        "P11229",  # Muscarinic acetylcholine receptor M1 (CHRM1)
        "P08172",  # Muscarinic acetylcholine receptor M2 (CHRM2)
        "P20309",  # Muscarinic acetylcholine receptor M3 (CHRM3)
        "P08173",  # Muscarinic acetylcholine receptor M4 (CHRM4)
        "P08912"   # Muscarinic acetylcholine receptor M5 (CHRM5)
    ],
    "human dopamine receptors": [
        "P21728",  # D1 receptor (DRD1)
        "P14416",  # D2 receptor (DRD2)
        "P35462",  # D3 receptor (DRD3)
        "P21917",  # D4 receptor (DRD4)
        "P21918"   # D5 receptor (DRD5)
    ],
    "human endothelin receptors": [
        "P25101",  # Endothelin-1 receptor (EDNRA)
        "P24530"   # Endothelin receptor type B (EDNRB)
    ]}

# Choose which receptors to plot
which_receptors = "human opsins"
with open(json_file, "r", encoding="utf-8") as receptors:
    data = json.load(receptors)
up_codes = set(uniprot_codes.get(which_receptors, []))
filtered_data = [
    entry for entry in data if any(code in entry["header"] for code in up_codes)
]

# Open data file
with open("filtered_data.json", "w", encoding="utf-8") as f:
    json.dump(filtered_data, f, indent=4)

# Plot the loop histograms for selected receptors in JSON file
visualize_histograms("filtered_data.json", which_receptors)

# Plot the histogram for ratio_inner_to_total using the full JSON file
visualize_ratio_histogram(json_file)