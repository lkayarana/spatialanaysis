import argparse
import json
import os
import tarfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt


# ----------------------------
# Helpers
# ----------------------------
def extract_tar_gz(tar_path: Path, outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as tf:
        tf.extractall(outdir)
    return outdir


def find_one(root: Path, filename: str) -> Path:
    matches = list(root.rglob(filename))
    if not matches:
        raise FileNotFoundError(f"Could not find {filename} under {root}")
    if len(matches) > 1:
        # pick the shortest path (closest match)
        matches = sorted(matches, key=lambda p: len(str(p)))
    return matches[0]


def read_10x_h5_matrix(h5_path: Path):
    """
    Reads 10x Genomics HDF5 matrix (v3-style).
    Returns:
      X: sparse CSC (genes x barcodes) as stored in file -> convert to cells x genes
      barcodes: list[str]
      features: pd.DataFrame with columns: id, name, feature_type
    """
    with h5py.File(h5_path, "r") as f:
        if "matrix" not in f:
            raise ValueError("Unsupported 10x H5: missing 'matrix' group.")
        g = f["matrix"]

        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        shape = g["shape"][:]  # (n_features, n_barcodes)

        X = csc_matrix((data, indices, indptr), shape=shape)

        barcodes = [b.decode("utf-8") for b in g["barcodes"][:]]

        fg = g["features"]
        feat_id = [b.decode("utf-8") for b in fg["id"][:]]
        feat_name = [b.decode("utf-8") for b in fg["name"][:]]
        feat_type = [b.decode("utf-8") for b in fg["feature_type"][:]]

        features = pd.DataFrame(
            {"id": feat_id, "name": feat_name, "feature_type": feat_type}
        )

    return X, barcodes, features


def get_gene_index(features: pd.DataFrame, gene: str):
    # 10x stores gene symbol in "name"
    hits = np.where(features["name"].values == gene)[0]
    if len(hits) == 0:
        return None
    return int(hits[0])


def sparse_col_to_dense(X_csc, col_idx: int) -> np.ndarray:
    # X is (cells x genes) CSC? We'll store in CSR? We'll manage carefully.
    # We will keep X as CSC for efficient gene (column) slicing if X is cells x genes CSC.
    col = X_csc[:, col_idx]
    return np.asarray(col.toarray()).ravel()


def ensure_outdir(outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)


# ----------------------------
# Main analysis
# ----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Path to P1_8um_min.tar.gz")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--periphery_um", type=float, default=50.0, help="Periphery radius in microns")
    ap.add_argument("--tumor_q", type=float, default=0.90, help="Quantile threshold for tumor_score")
    ap.add_argument("--mac_thresh", type=float, default=1.0, help="Threshold for macrophage marker (log1p counts)")
    args = ap.parse_args()

    tar_path = Path(args.input).resolve()
    outdir = Path(args.outdir).resolve()
    ensure_outdir(outdir)

    workdir = outdir / "_extracted"
    extract_tar_gz(tar_path, workdir)

    h5_path = find_one(workdir, "filtered_feature_bc_matrix.h5")
    pos_path = find_one(workdir, "tissue_positions.parquet")
    scale_path = find_one(workdir, "scalefactors_json.json")

    # Load matrix
    X_gxB, barcodes, features = read_10x_h5_matrix(h5_path)

    # Convert to cells x genes (barcodes x genes) as CSC for fast gene slicing
    X = X_gxB.T.tocsc()

    # Load positions
    pos = pd.read_parquet(pos_path)
    # Pick coordinate columns (10x often uses pxl_* instead of x/y)
    XCOL = "pxl_col_in_fullres"
    YCOL = "pxl_row_in_fullres"
    
    # 10x column for tissue membership
    TISSUECOL = "in_tissue"

    # Expect columns: barcode, x, y, tissue (as you saw)
    pos = pos.set_index("barcode")
    pos = pos.loc[barcodes]  # align

    # Basic library size and log1p normalization (very lightweight)
    # We need log1p counts per 10k similar to common practice.
    libsize = np.asarray(X.sum(axis=1)).ravel()
    libsize_safe = np.where(libsize == 0, 1, libsize)
    scale = 1e4 / libsize_safe
    # log1p(normalized) for selected genes only (avoid densifying whole matrix)
    # We'll compute gene vectors on demand.

    def gene_log1p(gene: str) -> np.ndarray:
        idx = get_gene_index(features, gene)
        if idx is None:
            return np.zeros(X.shape[0], dtype=float)
        raw = sparse_col_to_dense(X, idx)
        norm = raw * scale
        return np.log1p(norm)

    # Tumor score using CEACAM5/6 if present
    tumor_score = np.zeros(X.shape[0], dtype=float)
    used_tumor_markers = []
    for g in ["CEACAM5", "CEACAM6"]:
        if get_gene_index(features, g) is not None:
            tumor_score += gene_log1p(g)
            used_tumor_markers.append(g)

    if not used_tumor_markers:
        raise RuntimeError("Neither CEACAM5 nor CEACAM6 found in matrix. Cannot define tumor bins.")

    thr = np.quantile(tumor_score, args.tumor_q)
    is_tumor = tumor_score > thr

    # Distances and periphery (<=50um from any tumor bin)
    xy = pos[[XCOL, YCOL]].values.astype(float)
    tumor_xy = xy[is_tumor]
    tree = cKDTree(tumor_xy)
    dist, _ = tree.query(xy, k=1)

    region = np.where(is_tumor, "tumor", np.where(dist <= args.periphery_um, "periphery", "tissue"))

    # Macrophage bins using C1QC marker if present
    c1qc = gene_log1p("C1QC")
    is_mac = c1qc > args.mac_thresh

    # Periphery macrophages subtype split by SPP1 vs SELENOP
    spp1 = gene_log1p("SPP1")
    selenop = gene_log1p("SELENOP")

    periph_mac_mask = (region == "periphery") & is_mac
    mac_subtype = np.full(X.shape[0], "", dtype=object)
    mac_subtype[periph_mac_mask] = np.where(
        spp1[periph_mac_mask] > selenop[periph_mac_mask],
        "SPP1+",
        "SELENOP+",
    )

    # ----------------------------
    # Tables
    # ----------------------------
    df_bins = pd.DataFrame(
        {
            "barcode": barcodes,
            "x": pos[XCOL].values,
            "y": pos[YCOL].values,
            "in_tissue": pos[TISSUECOL].astype(bool).values,
            "tumor_score": tumor_score,
            "is_tumor": is_tumor,
            "dist_to_tumor_um": dist,
            "region": region,
            "C1QC_log1p": c1qc,
            "is_macrophage": is_mac,
            "SPP1_log1p": spp1,
            "SELENOP_log1p": selenop,
            "macrophage_subtype": mac_subtype,
        }
    )

    # Region counts
    tbl_region = df_bins["region"].value_counts().rename_axis("region").reset_index(name="n_bins")

    # Periphery vs tissue macrophage prevalence
    def pct_true(mask):
        return 100.0 * np.mean(mask) if len(mask) else 0.0

    tbl_mac_prevalence = pd.DataFrame(
        [
            {
                "region_group": "periphery",
                "n_bins": int(np.sum(region == "periphery")),
                "n_macrophage_bins": int(np.sum((region == "periphery") & is_mac)),
                "pct_macrophage_bins": pct_true(is_mac[region == "periphery"]),
            },
            {
                "region_group": "tissue",
                "n_bins": int(np.sum(region == "tissue")),
                "n_macrophage_bins": int(np.sum((region == "tissue") & is_mac)),
                "pct_macrophage_bins": pct_true(is_mac[region == "tissue"]),
            },
        ]
    )

    # Macrophage subtype counts (periphery macrophages only)
    sub = df_bins[df_bins["macrophage_subtype"].isin(["SPP1+", "SELENOP+"])]
    tbl_subtypes = sub["macrophage_subtype"].value_counts().rename_axis("macrophage_subtype").reset_index(name="n_bins")

    # Median distances by region
    tbl_dist = df_bins.groupby("region")["dist_to_tumor_um"].median().reset_index(name="median_dist_um")

    # DE-style stats on key genes between SPP1+ and SELENOP+ periphery macrophages
    genes_for_stats = ["SPP1", "SELENOP", "APOC1", "LPL", "CHI3L1", "STAB1", "SLC40A1", "MMP7", "IL1RN", "FN1"]
    stats_rows = []
    spp_mask = sub["macrophage_subtype"].values == "SPP1+"
    sel_mask = sub["macrophage_subtype"].values == "SELENOP+"

    for g in genes_for_stats:
        vec = gene_log1p(g)
        x = vec[periph_mac_mask][spp_mask] if np.any(periph_mac_mask) else np.array([])
        y = vec[periph_mac_mask][sel_mask] if np.any(periph_mac_mask) else np.array([])
        if len(x) >= 10 and len(y) >= 10:
            stat, p = mannwhitneyu(x, y, alternative="two-sided")
            stats_rows.append(
                {
                    "gene": g,
                    "median_SPP1plus": float(np.median(x)),
                    "median_SELENOPplus": float(np.median(y)),
                    "MWU_p_value": float(p),
                }
            )

    tbl_gene_stats = pd.DataFrame(stats_rows).sort_values("MWU_p_value")

    # Save Excel
    excel_path = outdir / "P1_spatial_report_outputs.xlsx"
    with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
        tbl_region.to_excel(writer, sheet_name="Region_counts", index=False)
        tbl_mac_prevalence.to_excel(writer, sheet_name="Macrophage_prevalence", index=False)
        tbl_subtypes.to_excel(writer, sheet_name="Macrophage_subtypes", index=False)
        tbl_dist.to_excel(writer, sheet_name="Median_distance_um", index=False)
        tbl_gene_stats.to_excel(writer, sheet_name="Subtype_gene_stats", index=False)
        df_bins.to_excel(writer, sheet_name="Bin_level_table", index=False)

    # ----------------------------
    # Figures
    # ----------------------------
    # Fig 1: Tumor & periphery map
    plt.figure(figsize=(6, 6))
    for lab, color, s in [("tissue", "lightgray", 2), ("periphery", "dodgerblue", 2), ("tumor", "red", 2)]:
        m = df_bins["region"].values == lab
        plt.scatter(df_bins.loc[m, "x"], df_bins.loc[m, "y"], s=s, c=color, label=lab)
    plt.gca().invert_yaxis()
    plt.legend(markerscale=3)
    plt.title(f"Tumor & Periphery (≤{args.periphery_um} µm), tumor markers={'+'.join(used_tumor_markers)}")
    plt.tight_layout()
    plt.savefig(outdir / "Fig1_tumor_periphery.png", dpi=250)
    plt.close()

    # Fig 2: Macrophage subtypes spatial
    plt.figure(figsize=(6, 6))
    m_all = (df_bins["region"].values == "periphery") & df_bins["is_macrophage"].values
    plt.scatter(df_bins.loc[m_all, "x"], df_bins.loc[m_all, "y"], s=3, c="lightgray", label="Periphery macrophage bins")

    for lab, color in [("SELENOP+", "green"), ("SPP1+", "orange")]:
        m = df_bins["macrophage_subtype"].values == lab
        plt.scatter(df_bins.loc[m, "x"], df_bins.loc[m, "y"], s=6, c=color, label=lab)

    plt.gca().invert_yaxis()
    plt.legend(markerscale=3)
    plt.title("Macrophage subtypes in tumor periphery (marker-based)")
    plt.tight_layout()
    plt.savefig(outdir / "Fig2_macrophage_subtypes.png", dpi=250)
    plt.close()

    # Fig 3: Distance distribution (periphery vs tissue)
    plt.figure(figsize=(6, 4))
    per = df_bins.loc[df_bins["region"] == "periphery", "dist_to_tumor_um"].values
    tis = df_bins.loc[df_bins["region"] == "tissue", "dist_to_tumor_um"].values
    plt.hist(per, bins=40, alpha=0.7, label="periphery")
    plt.hist(tis, bins=40, alpha=0.7, label="tissue")
    plt.xlabel("Distance to nearest tumor bin (µm)")
    plt.ylabel("Number of bins")
    plt.title("Distance distribution: periphery vs tissue")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "Fig3_distance_hist.png", dpi=250)
    plt.close()

    # Fig 4: Key gene boxplots for subtypes (log1p normalized)
    key_genes = ["SPP1", "SELENOP", "STAB1", "APOC1"]
    data = []
    labels = []
    for g in key_genes:
        vec = gene_log1p(g)[periph_mac_mask]
        x = vec[spp_mask]
        y = vec[sel_mask]
        data.extend([x, y])
        labels.extend([f"{g}\nSPP1+", f"{g}\nSELENOP+"])

    plt.figure(figsize=(10, 5))
    plt.boxplot(data, labels=labels, showfliers=False)
    plt.ylabel("log1p(normalized counts)")
    plt.title("Key genes in periphery macrophage subtypes (marker-based)")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(outdir / "Fig4_key_gene_boxplots.png", dpi=250)
    plt.close()

    # Save a short JSON summary for the report
    summary = {
        "input_tar": str(tar_path),
        "tumor_markers_used": used_tumor_markers,
        "tumor_score_quantile_threshold": float(args.tumor_q),
        "periphery_um": float(args.periphery_um),
        "macrophage_marker": "C1QC",
        "macrophage_threshold_log1p": float(args.mac_thresh),
        "n_bins_total": int(df_bins.shape[0]),
        "n_bins_tumor": int(np.sum(df_bins["region"].values == "tumor")),
        "n_bins_periphery": int(np.sum(df_bins["region"].values == "periphery")),
        "n_bins_tissue": int(np.sum(df_bins["region"].values == "tissue")),
        "n_periphery_macrophage_bins": int(np.sum(periph_mac_mask)),
        "n_SPP1plus_bins": int(np.sum(df_bins["macrophage_subtype"].values == "SPP1+")),
        "n_SELENOPplus_bins": int(np.sum(df_bins["macrophage_subtype"].values == "SELENOP+")),
        "outputs": {
            "excel": str(excel_path),
            "figures": [
                "Fig1_tumor_periphery.png",
                "Fig2_macrophage_subtypes.png",
                "Fig3_distance_hist.png",
                "Fig4_key_gene_boxplots.png",
            ],
        },
    }
    with open(outdir / "run_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\nDONE ✅")
    print(f"Excel: {excel_path}")
    print(f"Figures saved in: {outdir}")
    print("Summary JSON: run_summary.json")


if __name__ == "__main__":
    main()
