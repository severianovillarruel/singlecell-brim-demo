import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

# --- setup ---
root = Path(__file__).resolve().parents[1]
data_dir = root / "data"
fig_dir = root / "figures"
data_dir.mkdir(parents=True, exist_ok=True)
fig_dir.mkdir(parents=True, exist_ok=True)

# --- synthetic data ---
np.random.seed(42)
n_cells, n_genes = 600, 400
conditions = np.random.choice(["CR_PR", "SD_PD"], size=n_cells)
gene_names = [f"Gene{i:03d}" for i in range(n_genes)]
cell_ids = [f"Cell{i:04d}" for i in range(n_cells)]

baseline = np.random.lognormal(mean=1.4, sigma=0.6, size=n_genes)
cluster_ids = np.random.randint(0, 5, size=n_cells)
cluster_effects = np.random.normal(0, 0.35, size=(5, n_genes))

de_idx = np.random.choice(n_genes, size=50, replace=False)
cond_effect = np.zeros(n_genes)
cond_effect[de_idx[:25]] = 0.6     # up in CR_PR
cond_effect[de_idx[25:]] = -0.6    # down in CR_PR

counts = np.zeros((n_cells, n_genes), float)
for i in range(n_cells):
    mu = baseline + cluster_effects[cluster_ids[i]]
    mu = mu + (cond_effect if conditions[i] == "CR_PR" else -cond_effect) * 0.5
    counts[i] = np.random.lognormal(mean=mu, sigma=0.5)

# scale and Poisson to look like UMIs
counts = np.random.poisson(np.clip(counts / counts.mean() * 5, 0, None))
counts_df = pd.DataFrame(counts, index=cell_ids, columns=gene_names)
meta_df = pd.DataFrame({"cluster": cluster_ids, "condition": conditions}, index=cell_ids)

counts_df.to_csv(data_dir / "synthetic_scRNA_counts.csv")
meta_df.to_csv(data_dir / "metadata.csv")

# --- embedding + clustering ---
emb = PCA(n_components=2, random_state=42).fit_transform(np.log1p(counts))
km = KMeans(n_clusters=5, n_init=10, random_state=42).fit(emb)
meta_df["km_cluster"] = km.labels_
meta_df.to_csv(data_dir / "metadata.csv")

plt.figure()
plt.scatter(emb[:, 0], emb[:, 1], s=8)
plt.title("PCA of synthetic scRNA-seq")
plt.savefig(fig_dir / "pca.png", bbox_inches="tight", dpi=160)
plt.close()

plt.figure()
plt.scatter(emb[:, 0], emb[:, 1], s=8)
plt.title("PCA with k-means clusters (overlay)")
plt.savefig(fig_dir / "pca_kmeans.png", bbox_inches="tight", dpi=160)
plt.close()

# --- simple DE ranking (Welch-style t-stat without SciPy) ---
lib = counts_df.sum(axis=1).values[:, None]
cpm = (counts_df.values / np.clip(lib, 1, None)) * 1e6
logx = np.log1p(cpm)

mask_cr = (meta_df["condition"].values == "CR_PR")
mask_sd = (meta_df["condition"].values == "SD_PD")
A, B = logx[mask_cr], logx[mask_sd]

na, nb = A.shape[0], B.shape[0]
ma, mb = A.mean(axis=0), B.mean(axis=0)
va, vb = A.var(axis=0, ddof=1), B.var(axis=0, ddof=1)
t_stat = (ma - mb) / np.sqrt(va / na + vb / nb + 1e-12)

de = pd.DataFrame({"gene": gene_names, "t_stat": t_stat}).sort_values("t_stat", ascending=False)
de.to_csv(data_dir / "de_results_global.csv", index=False)

print("✅ Generated data/, figures/, and DE table.")
