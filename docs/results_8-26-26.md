> ## ⚠️ SUPERSEDED — do not use these results
>
> This analysis used an astrocyte segmentation that was later shown to select
> **neuropil rather than cells**: the objects it produced were GFAP-*depleted*
> relative to surrounding tissue (0.81× the tissue mean). Section 3.1 below —
> "Astrocytes — reactive shift, reproducible across both routes" — does **not**
> survive correction.
>
> Use [`results_8-28-26.md`](results_8-28-26.md) instead. This file is kept only
> to document what was originally reported and why it changed.

# COMET ICV + IP analysis — 8 brains, 4 cell types

**Run date:** 2026-08-27 · **Panel:** 28-channel COMET · **Model:** CD8-knockout mouse

Segmentation and phenotyping of 8 whole-slide mouse brain images across two
delivery routes (ICV, IP) and two conditions (control, treated), covering
astrocytes, microglia, neurons and T cells.

---

## 1. Samples

| Route | Control | Treated |
|---|---|---|
| ICV | ICV_C1_3, ICV_C3_3 | ICV_T1_3, ICV_T2_3 |
| IP | IP_C1_3, IP_C2_3 | IP_T2_3, IP_T3_3 |

533 tiles (4000 px, 100 px overlap) → CellProfiler on Slurm → per-brain merge
with overlap-zone deduplication (3–7% dropped per class).

**1,038,059 cells** total — 474,175 ICV, 563,884 IP.

| | Astrocytes | Microglia | Neurons | T cells |
|---|---|---|---|---|
| cells analysed | 120,352 | 218,135 | 460,410 | 652 |
| clusters | 9 | 10 | 10 | 9 |

## 2. What changed from the previous (ICV-only) run

- **Extended from 4 ICV brains to all 8**, adding the IP cohort as a second route.
- **Neurons and T cells can now be phenotyped.** Previously they carried no
  expression or morphology data at all — their output was 12 columns of position
  and parentage. The `MeasureObjectIntensity` module measuring the phenotypic
  panel ran only on microglia and astrocytes; `Cells` (the parent of neurons and
  T cells) measured only CD3/CD4/CD8/NeuN. Extending that module's image list to
  all 15 markers fixed it. **No threshold, filter, gate or segmentation setting
  was altered** — object counts are bit-identical to the previous run, verified
  brain by brain.
- **Reduced feature set.** Clustering uses 6 morphology + 11 expression features
  rather than the ~90 AreaShape columns CellProfiler emits (Zernike, Hu moments,
  inertia tensors dropped).

Features: `Area`, `MeanRadius`, `MajorAxisLength`, `Eccentricity`, `FormFactor`,
`ObjectSkeleton_NumberTrunks` + `APOE, Arg1, CD11b, CD206, CD44, CD68, CD74,
GFAP, Iba1, MAA, SNA` (neurons and T cells additionally carry CD3, CD4, CD8, NeuN).

Clustering: scale → PCA → Harmony (batch = brain) → UMAP → Leiden (res 0.4).

## 3. Results

### 3.1 Astrocytes — reactive shift, reproducible across both routes

The clearest finding. Two clusters move in opposite directions, **independently
in ICV and IP**:

| cluster | marker/morphology profile | ICV C→T | IP C→T |
|---|---|---|---|
| **5** | GFAP 2.2×, FormFactor 2.2×, trunks 0.56× — hypertrophic, compact, few processes | **+2.9** | **+3.4** |
| **4** | low across the board (CD68 0.48×, CD74 0.56×, CD11b 0.74×) — quiescent | **−3.1** | −1.4 |

(values are percentage-point change in cluster frequency, treated minus control;
marker values are fold-change vs the across-cluster median)

Quiescent astrocytes decline and GFAP-high reactive astrocytes expand with
treatment. That this reproduces in two independently processed cohorts is the
main reason to take it seriously.

### 3.2 T cells — the CD8 cluster expands with treatment

Unsupervised clustering placed most CD8 cells in a single cluster, and that
cluster tracks treatment:

| cluster 8 | ICV-C | ICV-T | IP-C | IP-T |
|---|---|---|---|---|
| cells | 5 | **16** | 4 | **26** |

Cluster 8 contains 17 of the 36 CD8 cells in the dataset (33% CD8 vs 6% overall).
Its profile is CD8 2.6×, GFAP 2.7×, NeuN 0.6×. Clustering had no access to
subtype labels, so this is an independent recovery of the treatment effect.

### 3.3 Microglia — modest, mostly ICV

| cluster | profile | ICV C→T | IP C→T |
|---|---|---|---|
| 7 | CD206 4.5×, CD68 2.0×, CD74 1.7× — alternatively activated / perivascular-like | +2.5 | +0.4 |
| 9 | GFAP 3.5×, CD44 3.0× | +1.1 | +0.6 |

### 3.4 Neurons — no consistent effect

Largest shifts are ±1–2 percentage points and do not reproduce between routes.

## 4. Two methodological findings that affect interpretation

### 4.1 The IP batch is dimmer with compressed dynamic range

At the 99th percentile of tissue pixels, IP runs 0.5–0.65× of ICV for most
markers (CD206 0.22×, CD44 0.18×), while median/background runs *higher*. The
channel map itself was verified to transfer correctly between cohorts (myeloid
triad Iba1/CD11b/CD68 co-localises in both; secondary-only controls dark in both;
per-plane coverage rank agreement ρ = 0.84).

Because segmentation and gating run per tile *before* Harmony, **Harmony corrects
the clustering but not the counting**. Within-route control-vs-treated contrasts
are sound; ICV-vs-IP comparisons of absolute counts are not.

Reassuringly, the major cell classes are stable across cohorts as a fraction of
all cells — astrocytes 8.3–14.7%, microglia 18.8–24.4%, neurons 41.9–47.8%, with
no route pattern. Segmentation is robust to the intensity difference.

### 4.2 T-cell gating was re-derived (changed from the pipeline)

**The absolute CD3 gate produced unusable counts.** The pipeline calls T cells on
`CD3 MeanIntensity >= 0.022`. Because CD3 intensity shifts between batches, that
fixed cutoff lands at a different point in each brain's distribution — CD3 p99
across the 8 brains spans 0.0174–0.0271, straddling it. Counts swung **78-fold**:

| brain | % over gate | T cells | % of all cells |
|---|---|---|---|
| ICV_C1_3 | 0.18% | 94 | 0.085% |
| IP_C2_3 (**control**) | 3.22% | 2,606 | 2.19% |
| IP_T2_3 | 8.52% | 8,925 | **6.26%** |
| IP_T3_3 | 0.15% | 190 | 0.118% |

A control brain exceeded every treated brain; two brains of the same cohort and
treatment differed 47-fold; and 6.26% of all brain cells being T cells is not
biologically plausible.

**Fix:** `CD3 >= 2.0 × (that brain's median CD3)`, leaving NeuN (≤0.050), Area
(100–950) and FormFactor (≥0.70) at their validated values. This is the same
per-brain normalisation the pipeline already applies to CD8/CD4 subtyping, for
the same reason. Between-brain range drops 78× → 5.6×, giving **44–197 T cells
per brain, 652 total**.

**CD8 ratio set per cohort: ICV 1.0, IP 1.5.** The pipeline's single cutoff of
1.0 was validated on ICV, and that validation holds — across both ICV control
brains *no cell exceeds ratio 1.0 at all*, while ICV treated brains carry cells
up to 4.70. ICV needs no correction.

IP is different: its ratio distribution is compressed, IP control brains carry
cells up to 2.67, and treated cells pile up just above 1.0. Every CD8 false
positive in the dataset is in the IP cohort, almost all in IP_C2_3.

Applying IP's stricter cutoff to ICV would discard real CD8 cells for no
specificity gain, so the cutoff is per cohort. Sensitivity/specificity across
the options (CD3 gate fixed at 2.0×):

| CD8 rule | control false positives | treated CD8 |
|---|---|---|
| global 1.0 | 7 | 45 |
| global 1.25 | 4 | 38 |
| global 1.5 | 2 | 33 |
| global 2.0 | 1 | 26 |
| **ICV 1.0 / IP 1.5 (used)** | **2** | **34** |

Two alternatives were tested and rejected. Gating on **CD8 brightness relative to
each brain's own CD8 median** performs far worse — at 2× it yields 443 control
false positives against 356 in treated, and even at 5× its specificity is below
the ratio's. Requiring **both** a high ratio and CD8 brightness is redundant,
changing treated counts by one cell at equal specificity. CD4 and CD8 share
bleedthrough characteristics, so the ratio cancels what absolute intensity
cannot — which is why the pipeline used a ratio to begin with.

### Validation against the CD8-knockout design

These are CD8-knockout mice: CD8 T cells should be absent except where adoptive
transfer took, while CD4 T cells (CD3+) persist.

| | control brains | treated brains |
|---|---|---|
| CD8 counts | 0, 0, 0, **2** | 7, 11, 8, 8 |

Three of four controls at exactly zero; the fourth (IP_C2_3) carries 2.
**652 T cells total: 616 CD4, 36 CD8** — the CD4-dominant profile the genotype
predicts, at the expected order of magnitude.

**Ceiling on CD8 counts.** CD8 cells are a subset of the T-cell pool (44–197 per
brain), so realistic numbers here are single digits to mid-teens per brain.
Loosening the CD3 gate to enlarge that pool does not help — it moves treated CD8
from 26 to only 30 while returning T-cell counts to implausible values (at 1.6×,
ICV_C1_3 balloons to 1,562 T cells). The constraint is signal in the images, not
the threshold.

> **Open item:** the 2 CD8 calls in control IP_C2_3 survive the cutoff and warrant
> a look in the image before control CD8 is reported as zero. IP_C2_3 has been the
> one persistently noisy control at every threshold tested.

## 5. Limitations

- **n = 2 brains per group.** All differences here are descriptive. No
  significance testing was performed — per-cell tests would be pseudoreplication,
  and brain-level tests at n=2 vs 2 have no power. The astrocyte result carries
  weight because it reproduces across independent routes, not from a p-value.
- **Cross-route count comparisons are confounded** by acquisition batch (§4.1).
- **T-cell counts depend on the gating choice** in §4.2 and are not directly
  comparable to the previous ICV-only analysis, which used the absolute gate.
- Channel map verified indirectly, by marker co-localisation, since the tiles do
  not retain OME channel names.

## 6. Reproducing

```bash
# 1. merged CellProfiler output -> reduced-set per-sample CSVs
python3 make_cafe_csvs.py OUTDIR run_8-26-26/merged_*
#    --tcell-k 2.5     stricter CD3      --cd8-ratio 3.0  stricter CD8
#    --tcell-abs       pipeline's original absolute CD3 gate

# 2. cluster one cell type
python3 cluster_cells.py OUTDIR Microglia
python3 cluster_cells.py OUTDIR Tcells --resolution 0.4 --neighbors 15
```

Per cell type the outputs are UMAPs (clusters / group / sample / split-by-group /
T-cell subtype), a marker dotplot, a cluster-frequency barplot, count, frequency
and median-expression tables, and `adata_clustered.h5ad`.

`global_x` / `global_y` are carried on every cell, so spatial analyses (e.g.
proximity of reactive astrocytes to T cells) are possible without re-running
segmentation.
