#!/usr/bin/env python3
"""
Cluster one cell type across all brains: scale -> PCA -> Harmony -> UMAP -> Leiden.

Replaces the eight Interactive_*.py variants with a single parameterized script.
Reads the reduced-set CSVs and samples.csv manifest written by make_cafe_csvs.py,
so group assignment comes from the manifest rather than from parsing filenames.

Usage
-----
    python3 cluster_cells.py DATA_DIR CELLTYPE [options]

    python3 cluster_cells.py "csv files for analysis_8-26-26" Microglia
    python3 cluster_cells.py "csv files for analysis_8-26-26" Tcells --resolution 0.3

Options
-------
    --resolution FLOAT   Leiden resolution                (default 0.4)
    --neighbors INT      UMAP/graph n_neighbors           (default 30)
    --min-dist FLOAT     UMAP min_dist                    (default 0.5)
    --metric STR         neighbor metric                  (default euclidean)
    --n-pcs INT          PCA components, capped at n_markers-1 (default 10)
    --no-harmony         skip batch correction
    --groups A,B         restrict to these groups (e.g. ICV-C,ICV-T)
    --outdir PATH        output directory (default DATA_DIR/CELLTYPE/output)

Run inside an environment with scanpy, anndata, harmonypy, igraph, seaborn:
    python cluster_cells.py ...
"""
import argparse
import sys
import warnings
from pathlib import Path
from time import time

import matplotlib
matplotlib.use("Agg")            # headless: write PDFs, never open a window

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
import anndata
import scanpy as sc
import harmonypy as hm

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=100, facecolor="white")

# Columns that are metadata / annotation, never clustering features.
NON_FEATURE = {"cell_id", "global_x", "global_y", "Subtype", "CD8_CD4_ratio"}

# Stable colors so every cell type's plots read the same way.
GROUP_COLORS = {
    "ICV-C": "#9ecae1",   # light blue
    "ICV-T": "#2171b5",   # dark blue
    "IP-C":  "#fcae91",   # light red
    "IP-T":  "#cb181d",   # dark red
}


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("data_dir")
    p.add_argument("celltype")
    p.add_argument("--resolution", type=float, default=0.4)
    p.add_argument("--neighbors", type=int, default=30)
    p.add_argument("--min-dist", type=float, default=0.5)
    p.add_argument("--metric", default="euclidean")
    p.add_argument("--n-pcs", type=int, default=10)
    p.add_argument("--no-harmony", action="store_true")
    p.add_argument("--groups", default=None)
    p.add_argument("--outdir", default=None)
    return p.parse_args(argv)


def load(data_dir, celltype, groups_filter):
    data_dir = Path(data_dir)
    cdir = data_dir / celltype
    if not cdir.is_dir():
        sys.exit(f"ERROR: {cdir} does not exist. Run make_cafe_csvs.py first.")

    manifest_fp = data_dir / "samples.csv"
    if not manifest_fp.exists():
        sys.exit(f"ERROR: manifest not found at {manifest_fp}")
    manifest = pd.read_csv(manifest_fp)
    manifest = manifest[manifest["CellType"] == celltype]
    if groups_filter:
        manifest = manifest[manifest["Group"].isin(groups_filter)]
    if manifest.empty:
        sys.exit(f"ERROR: no samples for CellType={celltype} in {manifest_fp}")

    print("=" * 60)
    print(f"LOADING {celltype}")
    print("=" * 60)

    dfs = []
    for row in manifest.itertuples():
        fp = cdir / row.File
        if not fp.exists():
            print(f"  WARNING: {fp.name} listed in manifest but missing -- skipped")
            continue
        df = pd.read_csv(fp)
        df["SampleID"] = row.BrainID
        df["Group"] = row.Group
        df["Route"] = row.Route
        df["Condition"] = row.Condition
        dfs.append(df)
        print(f"  {row.BrainID:12s} {len(df):7,d} cells   Group={row.Group}")

    if not dfs:
        sys.exit("ERROR: no CSVs could be loaded.")

    all_df = pd.concat(dfs, ignore_index=True)
    print(f"\nTotal cells: {len(all_df):,} across {all_df['SampleID'].nunique()} brains")
    print(f"Groups: {all_df['Group'].value_counts().sort_index().to_dict()}")
    return all_df


def build_adata(all_df):
    meta_cols = ["SampleID", "Group", "Route", "Condition"]
    passthrough = [c for c in ("Subtype", "CD8_CD4_ratio", "global_x", "global_y")
                   if c in all_df.columns]
    marker_cols = [c for c in all_df.columns
                   if c not in NON_FEATURE and c not in meta_cols]

    if not marker_cols:
        sys.exit("ERROR: no feature columns found.")

    # Drop rows with non-finite features -- PCA cannot handle them.
    feat = all_df[marker_cols].apply(pd.to_numeric, errors="coerce")
    good = np.isfinite(feat.to_numpy(dtype=np.float64)).all(axis=1)
    n_drop = int((~good).sum())
    if n_drop:
        print(f"  Dropping {n_drop:,} cells with missing/non-finite features")
        all_df = all_df[good].reset_index(drop=True)
        feat = feat[good].reset_index(drop=True)

    print(f"Features ({len(marker_cols)}): {marker_cols}")

    obs = all_df[meta_cols + passthrough].copy().reset_index(drop=True)
    for c in meta_cols + [p for p in passthrough if p == "Subtype"]:
        obs[c] = obs[c].astype("category")

    adata = anndata.AnnData(
        X=feat.to_numpy(dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=marker_cols),
    )
    adata.obs_names = [str(i) for i in range(adata.n_obs)]
    return adata


def embed(adata, args, outdir):
    adata.write_h5ad(outdir / "adata_raw.h5ad")

    print("\nSCALING AND PCA")
    adata.obs["batch"] = adata.obs["SampleID"].values
    adata.raw = adata.copy()                     # unscaled values for dotplots
    sc.pp.scale(adata, max_value=10)

    n_pcs = min(args.n_pcs, adata.n_vars - 1, adata.n_obs - 1)
    sc.pp.pca(adata, n_comps=n_pcs, svd_solver="auto")
    print(f"  PCA -> {adata.obsm['X_pca'].shape}")

    n_batches = adata.obs["batch"].nunique()
    if args.no_harmony:
        print("\nSkipping Harmony (--no-harmony)")
    elif n_batches < 2:
        print(f"\nSkipping Harmony (only {n_batches} batch)")
    else:
        print(f"\nHARMONY BATCH CORRECTION over {n_batches} brains")
        t0 = time()
        ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "batch",
                            max_iter_harmony=20)
        adata.obsm["X_pca_original"] = adata.obsm["X_pca"].copy()
        adata.obsm["X_pca_harmony"] = ho.Z_corr.T
        adata.obsm["X_pca"] = adata.obsm["X_pca_harmony"]
        print(f"  Harmony applied ({time()-t0:.1f}s)")
        adata.write_h5ad(outdir / "adata_harmony.h5ad")

    # n_neighbors must stay below the population size (matters for T cells)
    nn = min(args.neighbors, max(2, adata.n_obs - 1))
    if nn != args.neighbors:
        print(f"  Reducing n_neighbors {args.neighbors} -> {nn} (only {adata.n_obs} cells)")

    print(f"\nNEIGHBORS / UMAP / LEIDEN  (res={args.resolution}, nn={nn}, "
          f"min_dist={args.min_dist}, metric={args.metric})")
    sc.pp.neighbors(adata, n_neighbors=nn, metric=args.metric, use_rep="X_pca")
    sc.tl.umap(adata, min_dist=args.min_dist, spread=1.0)
    sc.tl.leiden(adata, resolution=args.resolution, flavor="igraph", n_iterations=2)
    n_clusters = adata.obs["leiden"].nunique()
    print(f"  -> {n_clusters} clusters")
    return adata


def _clean(ax):
    ax.set_xlabel(""); ax.set_ylabel("")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)


def plots(adata, celltype, outdir):
    print("\nPLOTS")
    groups = [g for g in GROUP_COLORS if g in set(adata.obs["Group"])]
    adata.obs["Group"] = adata.obs["Group"].cat.set_categories(groups)
    adata.uns["Group_colors"] = [GROUP_COLORS[g] for g in groups]

    size = max(4, min(30, 300000 // max(adata.n_obs, 1)))

    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color="leiden", ax=ax, show=False, title="",
               legend_loc="on data", legend_fontsize=8, size=size, alpha=0.8)
    _clean(ax); plt.tight_layout()
    fig.savefig(outdir / "umap_clusters.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color="Group", ax=ax, show=False, title="",
               size=size, alpha=0.8)
    _clean(ax); plt.tight_layout()
    fig.savefig(outdir / "umap_group.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(14, 10))
    sc.pl.umap(adata, color="SampleID", ax=ax, show=False, title="",
               size=max(3, size // 2), alpha=0.8)
    _clean(ax); plt.tight_layout()
    fig.savefig(outdir / "umap_sample.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # clusters split by group -- the main visual comparison
    fig, axes = plt.subplots(1, len(groups), figsize=(7 * len(groups), 7),
                             squeeze=False)
    for ax, g in zip(axes[0], groups):
        sub = adata[adata.obs["Group"] == g]
        sc.pl.umap(sub, color="leiden", ax=ax, show=False,
                   title=f"{celltype} — {g}  (n={sub.n_obs:,})",
                   legend_loc="on data", legend_fontsize=8, size=size, alpha=0.7)
        _clean(ax)
    plt.tight_layout()
    fig.savefig(outdir / "umap_clusters_by_group.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)

    if "Subtype" in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.umap(adata, color="Subtype", ax=ax, show=False, title="",
                   size=size, alpha=0.9)
        _clean(ax); plt.tight_layout()
        fig.savefig(outdir / "umap_subtype.pdf", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # marker dotplot -- needs >1 cluster for a dendrogram
    markers = list(adata.var.index)
    if adata.obs["leiden"].nunique() > 1:
        try:
            sc.tl.dendrogram(adata, groupby="leiden")
            sc.pl.dotplot(adata, var_names=markers, groupby="leiden",
                          use_raw=True, standard_scale="var",
                          dendrogram=True, show=False)
            fig = plt.gcf()
            fig.savefig(outdir / "dotplot_markers.pdf", dpi=150, bbox_inches="tight")
            plt.close(fig)
        except Exception as e:
            print(f"  dotplot skipped: {e}")
    print(f"  wrote plots to {outdir}")


def tables(adata, outdir):
    print("\nTABLES")
    obs = adata.obs
    freq = obs.groupby(["SampleID", "Group", "leiden"], observed=True)\
              .size().reset_index(name="count")
    total = obs.groupby(["SampleID", "Group"], observed=True)\
               .size().reset_index(name="total")
    freq = freq.merge(total, on=["SampleID", "Group"])
    freq["frequency"] = 100 * freq["count"] / freq["total"]
    freq.to_csv(outdir / "cluster_freq_long.csv", index=False)

    order = sorted(freq["leiden"].unique(), key=int)

    fig, ax = plt.subplots(figsize=(14, 8))
    sns.barplot(data=freq, x="leiden", y="frequency", hue="Group", order=order,
                palette=GROUP_COLORS, ax=ax, capsize=0.1,
                err_kws={"linewidth": 1.5})
    ax.set_xlabel("Cluster", fontsize=14)
    ax.set_ylabel("Frequency (% of brain)", fontsize=14)
    ax.set_title("Cluster frequency by group", fontsize=16)
    ax.tick_params(axis="both", labelsize=12)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    plt.tight_layout()
    fig.savefig(outdir / "barplot_cluster_freq_by_group.pdf", dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    for val, name in (("count", "cluster_count_all_samples.csv"),
                      ("frequency", "cluster_freq_all_samples.csv")):
        piv = freq.pivot_table(index="leiden", columns="SampleID",
                               values=val, fill_value=0)
        piv.index = piv.index.astype(int)
        piv = piv.sort_index()
        piv.index.name = "Cluster"
        piv.to_csv(outdir / name)

    # median expression per cluster, on the unscaled values
    markers = list(adata.var.index)
    X = adata.raw.X if adata.raw is not None else adata.X
    if sp.issparse(X):
        X = X.toarray()
    expr = pd.DataFrame(np.asarray(X), columns=markers, index=obs.index)
    expr["leiden"] = obs["leiden"].astype(int).values
    expr["SampleID"] = obs["SampleID"].values
    expr["Group"] = obs["Group"].values

    expr.groupby("leiden", observed=True)[markers].median()\
        .to_csv(outdir / "median_expr_by_cluster.csv")
    expr.groupby(["leiden", "Group"], observed=True)[markers].median()\
        .to_csv(outdir / "median_expr_by_cluster_group.csv")
    expr.groupby(["leiden", "SampleID"], observed=True)[markers].median()\
        .to_csv(outdir / "median_expr_by_cluster_sample.csv")

    if "Subtype" in obs.columns:
        obs.groupby(["Group", "SampleID", "Subtype"], observed=True)\
           .size().reset_index(name="count")\
           .to_csv(outdir / "tcell_subtype_counts.csv", index=False)

    print(f"  wrote tables to {outdir}")


def main(argv):
    args = parse_args(argv)
    groups_filter = args.groups.split(",") if args.groups else None

    outdir = Path(args.outdir) if args.outdir \
        else Path(args.data_dir) / args.celltype / "output"
    outdir.mkdir(parents=True, exist_ok=True)

    all_df = load(args.data_dir, args.celltype, groups_filter)
    adata = build_adata(all_df)
    adata = embed(adata, args, outdir)
    plots(adata, args.celltype, outdir)
    tables(adata, outdir)

    adata.write_h5ad(outdir / "adata_clustered.h5ad")
    print(f"\nDone. {adata.n_obs:,} cells, {adata.n_vars} features, "
          f"{adata.obs['leiden'].nunique()} clusters")
    print(f"All outputs under: {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
