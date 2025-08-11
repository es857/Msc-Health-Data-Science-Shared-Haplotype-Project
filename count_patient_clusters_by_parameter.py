import pandas as pd
import glob

# Patient → (Gene, Ancestry)
patient_regions = {
    "WG0091": ("Chr19_haplotype", "European"),
    "WG0094": ("Chr19_haplotype", "European"),
    "WG0347": ("Chr19_haplotype", "European"),
    "WG0600": ("Chr19_haplotype", "European"),
    "WG1537": ("Chr19_haplotype", "European"),
    "WG0225": ("ZNF808", "European"),
    "WG0872": ("HK1", "European"),
    "WG0631": ("FICD", "Pakistani"),
    "WG0153": ("FICD", "Pakistani"),
    "WG0154": ("FICD", "Pakistani"),
    "WG0161": ("FICD", "Pakistani"),
    "WG1726": ("FICD", "Pakistani"),
    "WG0867": ("FICD", "Pakistani"),
    "WG1322": ("RNU4ATAC", "Pakistani"),
    "WG0596": ("RNU6ATAC", "Pakistani"),
    "WG0512": ("NARS2", "Iranian"),
    "WG1068": ("PAX4", "Iranian"),
    "WG0364": ("NARS2", "Iranian"),
    "WG0361": ("TARS2", "Iranian"),
    "WG0513": ("ZNF808", "Saudi"),
    "WG0367": ("PDIA6", "Saudi"),
    "WG1094": ("NARS2", "Saudi"),
    "WG0878": ("PAX4", "Indian"),
    "WG0158": ("TARS2", "Indian"),
    "WG1078": ("EIF2B1", "Indian"),
    "WG0363": ("TARS2", "Turkish"),
    "WG0117": ("YIPF5", "Turkish"),
    "WG1758": ("FICD", "Sudanese"),
    "WG0718": ("NARS2", "Sudanese"),
    "WG0366": ("ZNF808", "Sudanese")
}

def process_file(file_path):
    # Read TSV file
    df = pd.read_csv(file_path, sep="\t", header=None)
    
    # Name columns
    df.columns = [
        "chrom", "start", "end", "num_vars", "num_inds", "inds_str",
        "mode", "err_left", "err_right", "n_ultra_rare", "ultra_rare_label"
    ]
    
    # Split patient list and assign match status
    def get_match_status(inds):
        genes = set()
        for ind in inds:
            gene = patient_regions.get(ind, ("UNKNOWN",))[0]
            genes.add(gene)
        return "TRUE_MATCH" if len(genes) == 1 else "NO_MATCH"
    
    df["inds_list"] = df["inds_str"].apply(lambda x: [i.strip() for i in x.split(",") if i.strip()])
    df["match_status"] = df["inds_list"].apply(get_match_status)
    
    # Calculate percentage true match
    total_entries = len(df)
    true_match_count = (df["match_status"] == "TRUE_MATCH").sum()
    true_match_percent = (true_match_count / total_entries) * 100 if total_entries > 0 else 0
    
    return {
        "file_name": file_path,
        "total_entries": total_entries,
        "true_match_count": true_match_count,
        "true_match_percent": true_match_percent
    }

# Find all matching files
file_pattern = "shared_regions_dist*_vars*_rare*_no_gene_filter.tsv"
matching_files = glob.glob(file_pattern)

# Process all files and collect results
results = []
for file_path in matching_files:
    result = process_file(file_path)
    results.append(result)

# Create summary DataFrame
summary_df = pd.DataFrame(results)

# Print and save results
print("\nSummary of true matches across all files:")
print(summary_df[["file_name", "total_entries", "true_match_count", "true_match_percent"]].to_string(index=False))

# Save detailed results to CSV
summary_df.to_csv("shared_regions_match_summary.csv", index=False)
print("\nDetailed results saved to 'shared_regions_match_summary.csv'")