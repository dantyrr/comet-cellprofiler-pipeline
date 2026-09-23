# COMET ICV + IP — corrected analysis, 8 brains, 4 cell types

**Supersedes [`results_8-26-26.md`](results_8-26-26.md), which should not be used.**
That analysis rested on an astrocyte segmentation later shown to be selecting
neuropil rather than cells, and its headline finding does not survive correction.

Run 2026-08-28 · 8 brains · nucleus-seeded segmentation
([`pipeline/COMET_28ch_seeded.cppipe`](../pipeline/COMET_28ch_seeded.cppipe))

---

## 1. What was wrong, and how it was found

The original pipeline identified astrocytes by thresholding GFAP directly
(adaptive Minimum Cross-Entropy, floor 0.02). Checking the resulting objects
against the image showed they were **not GFAP-enriched at all**:

| | object's own marker | tissue baseline | enrichment |
|---|---|---|---|
| Microglia → Iba1 | 11.55 | 1.70 | **6.80×** |
| Astrocytes → GFAP | 6.22 | 7.67 | **0.81× (depleted)** |

Across four brains, only 5.7–16.8% of "astrocyte" objects exceeded the tissue
95th percentile for GFAP, and their median sat at the tissue median.

**Cause.** The threshold floor was set *below* GFAP's own tissue background. In
8-bit units the astrocyte floor is 5.10 while GFAP's tissue median is 6.0 and
mean 7.67; the microglia floor is 2.55 against an Iba1 median of 0.0. The same
design succeeds for a sparse marker and fails for a dense one — GFAP fine
processes tile essentially the whole neuropil, so adaptive thresholding in a
50-px window always finds "foreground". Global thresholding does not help:
Minimum Cross-Entropy on log(GFAP) also lands at ~6, because on this channel the
background *is* GFAP signal.

Background subtraction was already applied upstream in Horizon Viewer, and GFAP
correlates only 0.049 with the Cy5 control channel, so the floor is real
biological signal rather than autofluorescence.

## 2. The fix — nucleus-seeded segmentation

Each nucleus is grown into marker-positive territory
(`IdentifySecondaryObjects`, Propagation) instead of thresholding the marker and
hoping the objects are cells:

- **Astrocytes** = Nuclei → GFAP, threshold 0.047, then filter GFAP ≥ 0.047
- **Microglia** = Nuclei → Iba1, threshold 0.012, then filter Iba1 ≥ 0.012
- **Neurons** = `Cells` (nucleus+3px) with NeuN ≥ 0.03 — **unchanged**, this was
  already sound: NeuN is nuclear, so a tight window is correct, and the gate sits
  in a real valley (all-cells median 2.35 → 23.42 when gated)

**One nucleus per object by construction**, which also eliminates the merged
process networks the old pipeline produced (objects with 39 primary processes,
~10% of the dataset, where one object spanned several cells).

| | old | new |
|---|---|---|
| Astrocyte median GFAP | 0.024 (depleted) | **0.113 – 0.142** (~5× enriched) |
| Astrocytes per brain | 8,426 – 20,524 | 2,529 – 6,809 |
| Microglia per brain | 21,227 – 34,802 | 13,168 – 22,285 |

The microglia drop is a correction, not a loss: 21,227 of 101,812 cells is 20.8%
microglia, which is implausible; the new figure is ~5.7%, within the expected
5–10%. `Cells` (101,812) and `Neurons` (42,644) are **bit-identical** to the old
run, confirming the change was confined to the glia.

Final population sizes: **35,480 astrocytes · 132,271 microglia · 460,410 neurons
· 1,162 T cells** (1,115 CD4 / 47 CD8).

## 3. The main finding — activated microglia sit in a distinct niche

This is the most robust result in the project, because it is a **within-brain**
comparison: each brain is its own control, so density, batch, section quality and
animal-to-animal variation all cancel.

Splitting microglia at each brain's median CD68 and comparing neighbourhood
composition within 25 µm:

| neighbour type | high/low ratio | brains agreeing |
|---|---|---|
| **T cells** | **3.82×** | **8/8** |
| **Astrocytes** | **1.98×** | **8/8** |
| Microglia | 1.32× | 8/8 |
| **Neurons** | **0.74×** | 8/8 |

**The effect decays monotonically with radius** — T cells 4.94× at 15 µm, 3.82×
at 25, 2.61× at 50, 2.38× at 100, 2.31× at 200. That decay is the signature of a
genuine *local* association diluting as more surrounding tissue is sampled; a
regional confound would be roughly flat.

### It tracks activation, not brightness

Splitting on each marker in turn (25 µm, T-cell enrichment):

| split marker | T cells | brains agreeing |
|---|---|---|
| **CD74** (MHC-II invariant chain) | **5.91×** | 8/8 |
| Arg1 | 4.49× | 8/8 |
| CD44 | 4.00× | 8/8 |
| CD68 | 3.82× | 8/8 |
| CD206 | 2.54× | 7/8 |
| APOE | 1.75× | 6/8 |
| **Iba1** | **0.90×** | 5/8 — **none** |

**Iba1 is the control that matters.** It is the microglial identity marker, so
Iba1-high means more microglial mass, not activation — and it shows no T-cell
association whatsoever. If the effect were an artifact of brighter objects
sitting in denser neighbourhoods, Iba1 would show it. The ranking runs
antigen-presentation > phagocytic/activation > alternative activation >> identity.

*(GFAP as a split marker gives 12× for astrocyte neighbours and must be ignored —
a microglial object measured as GFAP-high is one with astrocyte signal bleeding
into its mask, so the result is circular.)*

### Both CD4 and CD8

Asking the question from the T-cell side — for each T cell, what fraction of
nearby microglia are CD68-high, against a 50% null:

| radius | CD4 | CD8 |
|---|---|---|
| 25 µm | **83.0%** (p = 1e-137) | **83.6%** (p = 1e-07) |
| 50 µm | 80.7% | 80.0% |
| 100 µm | 77.6% | 75.1% |

Per-cell (so one T cell in a crowded spot cannot dominate): 557/642 CD4 and
**21/27 CD8** have a majority-CD68-high neighbourhood, p = 0.006 for CD8. The
effect is large enough to reach significance on 29 CD8 cells, and is not a CD4
effect that CD8 is riding on.

## 4. Treated vs control — what survives, and what does not

### Marker-level comparison (clustering-independent)

Per-brain median marker intensity, treated/control ratio within route:

| | astrocytes ICV / IP | microglia ICV / IP | |
|---|---|---|---|
| **CD206** | 0.65 / **0.50** | 0.72 / **0.45** | **down in both, ~2× in IP** |
| **APOE** | 0.93 / 0.76 | 0.90 / 0.76 | down in both |
| Iba1 | — | 1.08 / 1.05 | up in both |
| CD68 | — | 1.43 / 1.05 | up in both |
| GFAP | 1.02 / 1.06 | — | up in both, but trivially |
| CD74 | 2.18 / 0.93 | 2.03 / 1.00 | **ICV only** |
| CD44 | 1.76 / 1.07 | 1.97 / 0.94 | ICV only |

CD206 falling ~2-fold in both routes and both cell types is the clearest
treated-vs-control signal. CD74 and CD44 rising ~2× in ICV with no IP effect is
where the two delivery routes differ most.

### What does NOT survive

**The astrocyte GFAP/reactivity finding from `results_8-26-26.md`.** With
corrected segmentation it is not reproducible, and the per-brain values show why:

| | ICV-C | ICV-T | IP-C | IP-T |
|---|---|---|---|---|
| astrocyte % of neighbours | 10.85, 8.72 | 9.06, **16.79** | 6.53, 8.05 | 9.86, 6.03 |
| astrocyte share of glia | 23.4, 20.0 | 21.2, **34.1** | 14.6, 19.6 | 24.8, 12.9 |

One treated brain per route sits inside the control range, and the IP treated
brains fall on opposite sides of their controls. **Within-group spread exceeds
the between-group difference.** Direct measurement agrees: astrocyte median GFAP
changes by only 2–6%.

**Cluster-frequency comparisons in general.** Three clustering configurations
gave three different answers for the GFAP-high population (−0.22/+2.02,
−1.88/+0.81, +4.26/+0.66). The combined-run result turned out to be tracking the
astrocyte:microglia ratio rather than astrocyte state. At n=2 per group, with
partitions that reshuffle when the feature set changes, cluster frequencies are
not a sound basis for treated-vs-control claims.

**Neighbourhood cell counts, unless normalised.** Total neighbours per microglion
rose 1.208× (ICV) and 1.198× (IP) — a uniform density increase that inflates every
cell type equally. Always compare composition (percentages), not raw counts.

## 5. Limitations

- **n = 2 brains per group.** Everything treated-vs-control here is descriptive.
  No significance testing was done: per-cell tests would be pseudoreplication, and
  brain-level tests at 2 vs 2 have no power. The spatial result is different — it
  is within-brain and replicates in all 8, which is why it carries weight.
- **GFAP-based segmentation under-counts astrocytes.** GFAP strongly labels
  fibrous and reactive astrocytes; protoplasmic astrocytes in gray matter are
  GFAP-dim and fall below any usable threshold. These are "GFAP+ astrocytes".
- **The spatial result is correlational.** It says activated microglia and T cells
  are co-located, not which recruits which.
- **No anatomical annotation.** Coordinates are raw pixels with no region labels,
  so a regional confound cannot be fully excluded (the radius decay argues
  strongly against one).
- **Pixel size is assumed at 0.28 µm/px**, inferred rather than read from OME
  `PhysicalSizeX`. Stated distances scale linearly with it; ratios and p-values
  do not.
- **IP_T3_3 reads as an outlier** on several markers (APOE 0.0132 against
  0.0165–0.0223 elsewhere). With two brains per group, one outlier moves a group
  mean substantially.

## 6. Reproducing

```bash
# 1. merged CellProfiler output -> reduced-set per-sample CSVs
python3 analysis/make_cafe_csvs.py OUTDIR run_*/merged_* --extra-morphology

# 2. cluster one cell type, or several jointly
python3 analysis/cluster_cells.py OUTDIR Microglia
python3 analysis/cluster_cells.py OUTDIR Microglia,Astrocytes \
        --cluster-features "expression,AreaShape_MeanRadius"

# 3. the spatial analysis
python3 analysis/spatial_neighborhood.py OUTDIR --radius-um 25 --split-by CD68
```

See [`RUNBOOK.md`](RUNBOOK.md) for the full pipeline and
[`spatial_analysis.md`](spatial_analysis.md) for the spatial method.
