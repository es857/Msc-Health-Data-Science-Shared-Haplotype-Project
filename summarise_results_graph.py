import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob


INPUT_FILES = "shared_regions_dist*_vars*_rare*_no_gene_filter.tsv"
GNOMAD_THRESHOLD = 0.0001 
def load_data():
    files = glob(INPUT_FILES)
    data = []
    
    for f in files:
        parts = f.replace("shared_regions_", "").replace(".tsv", "").split("_")
        params = {
            'cluster_dist': int(parts[0][4:]),
            'min_variants': int(parts[1][4:]),
            'min_ultra_rare': int(parts[2][4:])
        }
        
        # Count regions 
        with open(f) as fd:
            regions = sum(1 for line in fd if line.strip())
        
        data.append({**params, 'regions': regions})
    
    return pd.DataFrame(data)

def plot_results(df):
    plt.figure(figsize=(10, 6))
    sns.set_style("whitegrid")
    
    # Create lineplot with min_ultra_rare as hue
    ax = sns.lineplot(
        data=df,
        x='cluster_dist',
        y='regions',
        hue='min_ultra_rare',
        style='min_variants',
        markers=True,
        palette="viridis",
        linewidth=2.5
    )
    

    ax.axhline(y=0, color='red', linestyle='--', alpha=0.3)
    
    plt.title(f"Neonatal Diabetes: Shared Ultra-Rare Regions (AF ≤ {GNOMAD_THRESHOLD})")
    plt.xlabel("Cluster Distance (bp)")
    plt.ylabel("Regions Identified")
    plt.legend(title="Min Ultra-Rare Variants\n(per region)")
    
    # Annotate biological relevance
    plt.text(0.5, 0.95, 
             "Higher min_ultra_rare → Fewer but more confident regions",
             transform=ax.transAxes, ha='center', fontsize=10)
    
    plt.tight_layout()
    plt.savefig("neonatal_dm_regions.png", dpi=300)

if __name__ == "__main__":
    df = load_data()
    plot_results(df)
