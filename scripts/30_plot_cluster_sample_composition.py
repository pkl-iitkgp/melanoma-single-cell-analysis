from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

input_file = Path("results/annotation/cluster_sample_composition.csv")
output_file = Path("figures/clustering/cluster_sample_composition.png")

df = pd.read_csv(input_file)
df["cluster"] = df["cluster"].astype(str)
df = df.sort_values("cluster", key=lambda s: s.astype(int))

fig, ax = plt.subplots(figsize=(11, 5.5))
bars = ax.bar(
    df["cluster"],
    df["dominant_sample_fraction"] * 100,
    edgecolor="black",
    linewidth=0.4,
)
ax.axhline(80, linestyle="--", linewidth=1, label="80% reference")
ax.set_xlabel("Leiden cluster")
ax.set_ylabel("Cells from dominant tumor sample (%)")
ax.set_title("Tumor-sample concentration within each Leiden cluster")
ax.set_ylim(0, 108)
ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()

output_file.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output_file, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {output_file}")
