import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import re

# Load all TSV files matching pattern
files = glob.glob("shared_regions_dist*_vars*_rare*_no_gene_filter.tsv")

# Pattern to extract parameters
pattern = re.compile(r"dist(\d+)_vars(\d+)_rare([\d.]+)")

results = []

for file in files:
    match = pattern.search(file)
    if match:
        dist, vars_, rare = match.groups()
        df = pd.read_csv(file, sep="\t", header=None)
        df.columns = [
            "chrom", "start", "end", "n_variants", "n_individuals", "individuals",
            "mode", "err_left", "err_right", "n_ultra_rare", "filter"
        ]
        
        total = len(df)
        true_matches = 0

        for inds in df["individuals"]:
            ids = inds.split(",")
            prefixes = set(id[:5] for id in ids)
            if len(prefixes) == 1:
                true_matches += 1

        results.append({
            "dist": int(dist),
            "n_vars": int(vars_),
            "rare": float(rare),
            "total_regions": total,
            "true_matches": true_matches,
            "percent_true": 100 * true_matches / total if total > 0 else 0
        })

res_df = pd.DataFrame(results)

# --- Figure 1: Bar Plot ---
plt.figure(figsize=(10, 6))
sns.barplot(data=res_df, x="n_vars", y="percent_true", hue="dist")
plt.title("Percentage of True Matches by Parameters")
plt.ylabel("True Matches (%)")
plt.xlabel("Min Variants")
plt.legend(title="Cluster Distance")
plt.tight_layout()
plt.savefig("figure1_barplot.png")

# --- Figure 2: Sensitivity vs Specificity ---
plt.figure(figsize=(8, 6))
sns.scatterplot(data=res_df, x="percent_true", y="true_matches", hue="n_vars", palette="tab10", style="dist")
plt.title("Sensitivity vs Specificity")
plt.xlabel("Specificity (% True Matches)")
plt.ylabel("Sensitivity (True Match Count)")
plt.legend(title="Min Variants / Distance")
plt.tight_layout()
plt.savefig("figure2_tradeoff.png")

# --- Figure 3: Heatmap ---
heatmap_data = res_df.pivot_table(values="percent_true", index="n_vars", columns="dist")
plt.figure(figsize=(8, 6))
sns.heatmap(heatmap_data, annot=True, fmt=".1f", cmap="YlGnBu")
plt.title("Heatmap of True Match %")
plt.xlabel("Cluster Distance")
plt.ylabel("Min Variants")
plt.tight_layout()
plt.savefig("figure3_heatmap.png")

g = sns.FacetGrid(res_df, col="dist", hue="n_vars", palette="tab20", height=5)
g.map(sns.scatterplot, "percent_true", "true_matches", s=80)
g.add_legend(title="Min Variants")
g.set_axis_labels("Specificity (% True Matches)", "Sensitivity (True Match Count)")
plt.tight_layout()
plt.savefig("figure4_tradeoff.png")

plt.figure(figsize=(10, 6))
sns.scatterplot(
    data=res_df,
    x="percent_true",
    y="true_matches",
    hue="n_vars",
    size="dist",
    sizes=(50, 200),
    palette="plasma",
    alpha=0.8
)
plt.title("Parameter Impact on Detection Performance")
plt.xlabel("Specificity (% True Matches)")
plt.ylabel("Sensitivity (True Match Count)")
plt.legend(bbox_to_anchor=(1.05, 1))
plt.tight_layout()
plt.savefig("figure5_tradeoff.png")


# --- Figure 1: 3D Scatter Plot (Distance vs Variants vs Rare) ---
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Create color mapping based on percent_true
norm = plt.Normalize(res_df['percent_true'].min(), res_df['percent_true'].max())
colors = plt.cm.viridis(norm(res_df['percent_true']))

sc = ax.scatter(
    xs=res_df['dist'],
    ys=res_df['n_vars'],
    zs=res_df['rare'],
    c=colors,
    s=res_df['true_matches']*5,  # Size by true matches
    alpha=0.8,
    depthshade=False
)

ax.set_xlabel('Cluster Distance (bp)')
ax.set_ylabel('Min Variants')
ax.set_zlabel('Ultra-Rare Threshold')
ax.set_title('Parameter Space Exploration (Color: % True, Size: Sensitivity)')

# Add colorbar
cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='viridis'), ax=ax, shrink=0.5)
cbar.set_label('% True Matches')

plt.tight_layout()
plt.savefig("figure1_3d_parameter_space.png", dpi=300)


