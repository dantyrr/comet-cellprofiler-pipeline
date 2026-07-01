# Channel map and threshold rationale

This document captures the validated scientific decisions baked into
`pipeline/COMET_track1plus2_28ch.cppipe`, and *why* they were made. It's
meant for a future collaborator (or future you) who knows general
bioimage analysis but hasn't used this specific pipeline before. Values
here should not be changed without re-validating against real data — see
`COWORK_BUILD_BRIEF.md`'s note that these were empirically validated
across multiple test images, including positive/negative controls.

## 1. Channel map (28-channel COMET panel)

The pipeline was originally built against a 30-channel image. Production
data uses a 28-channel panel (2 fewer control channels). The map below was
verified against the OME-XML channel names embedded in the source file,
not assumed by offset arithmetic.

Confirmed marker → T-plane index map (0-based):

| Marker | T-plane index (0-based) |
|---|---|
| DAPI | 0 |
| CD206 | 5 |
| CD4 | 6 |
| NeuN | 7 |
| CD8 | 8 |
| CD3 | 10 |
| GFAP | 14 |
| APOE | 15 |
| Iba1 | 17 |
| Arg1 (Arginase1) | 19 |
| CD74 | 21 |
| CD11b | 23 |
| CD68 | 24 |

Full 28-plane panel also includes (not used by this pipeline, available if
extended): Anti-Rabbit_TRITC (1), Anti-Rat_Cy5 (2), CD11c (3), Foxp3 (4),
CD31 (9), anti-Mouse_TRITC (11), Anti-Rabbit_Cy5 (12), Granzyme K (13),
CD45 (16), CD45R/B220 (18), CD44 (20), Na+K+ATPase (22), CoxIV (25),
Mal I_FITC (26), SNL/EBL_Cy5 (27).

**Important:** if this pipeline is applied to a different COMET panel or
acquisition batch, re-verify the channel map from the OME-XML — do not
assume the offset transfers between batches/panels. Verification method:

```python
import tifffile
with tifffile.TiffFile(path) as tif:
    xml = tif.ome_metadata
# then parse <Channel Name="..."> entries in order
```

## 2. Microglia (Iba1) segmentation

- Threshold strategy: **Adaptive**
- Method: **Minimum Cross-Entropy**
- Log transform before thresholding: **Yes**
- Threshold correction factor: **1.0**
- Lower/upper bounds: **0.01 / 1.0**
- Adaptive window size: **50**
- Typical object diameter: 30–150 px, discard outside range

**Rationale:** Iba1 signal is diffuse/ramified with no clean background
valley on raw intensities — plain Minimum Cross-Entropy without a log
transform collapses to a near-zero threshold on sparse tiles (fuses
background haze into one giant object, which then fails the diameter
filter → zero detected cells). The log transform symmetrizes the
background/foreground distributions enough for MinXE to find a real
valley. The 0.01 floor prevents total collapse on very sparse/edge tiles
while still letting Adaptive raise the threshold appropriately on dense
tissue.

## 3. Astrocyte (GFAP) segmentation

- Threshold strategy: **Adaptive**
- Method: **Minimum Cross-Entropy**
- Log transform before thresholding: **Yes**
- Threshold correction factor: **1.0**
- Lower/upper bounds: **0.02 / 1.0**
- Adaptive window size: **50**
- Typical object diameter: 30–150 px, discard outside range

Same rationale as microglia. The floor value (0.02) is roughly double the
Iba1 floor, calibrated from tissue-masked whole-slide intensity
histograms (background-dominant exponential decay, threshold set near the
onset of the bright-tail signal).

A secondary `FilterObjects` step (`AstrocytesFiltered`) drops the dimmest
decile of detected astrocyte objects (`MeanIntensity_GFAP >= 0.024`) to
exclude marginal/possibly-artifactual detections while preserving the
full unfiltered `Astrocytes` object set (with skeleton/morphology
measurements) for anyone who wants the complete population.

## 4. T cell identification and CD4/CD8 subtyping

**Problem discovered during validation:** CD4 and CD8 (both Cy5-labeled in
different COMET cycles) show elevated baseline intensity relative to the
original pipeline-development image, such that absolute intensity gates
(e.g. `CD4 <= 0.025`) systematically reject real CD8 T cells whose CD4
channel sits slightly above that fixed cutoff due to acquisition-batch
intensity drift — **not** biological CD4 co-expression.

**Solution — ratio-based subtyping instead of absolute dual gating:**

1. `Tcells` = objects passing: `CD3 MeanIntensity >= 0.022`, `NeuN
   MeanIntensity <= 0.050`, `Area` in `[100, 950]` px, `FormFactor >=
   0.70`. This identifies T-cell-like objects independent of CD4/CD8
   subtype.
2. `CalculateMath` computes `CD8_CD4_ratio = MeanIntensity_CD8 /
   MeanIntensity_CD4` on the `Cells` object set (must run **before** the
   Tcells `FilterObjects` step so the ratio measurement is available to
   inherit/reference).
3. `CD8_Tcells` = Tcells-equivalent criteria (re-applied directly to
   `Cells`, since CellProfiler 4.2.8's `FilterObjects` GUI does not expose
   inherited Math measurements from a parent object set in its dropdown —
   a known tool limitation, not a design choice) **plus**
   `CD8_CD4_ratio >= 1.0`.
4. `CD4_Tcells` = same base criteria **plus** `CD8_CD4_ratio <= 1.0`.

**Validation performed (documented as evidence, not re-derived here):**

- *Negative control* — CD8-knockout mouse brain (no adoptive transfer): 23
  total T cells detected across 6 image tile-sets, **0 classified as
  CD8** — correct, since CD8 T cells cannot exist in this genotype.
- *Positive control* — CD8-knockout + adoptive CD8+ T cell transfer: 15
  total T cells detected across 5 tile-sets, **3 classified as CD8**
  (ratio 2.0–3.2, i.e. CD8 signal 2-3x brighter than CD4 in the same
  cell), 12 classified as CD4. The 3 CD8-positive cells were independently
  confirmed to be the single brightest CD8 signal in their respective
  images (CD8 MeanIntensity = image maximum in all 3 cases), NeuN
  near-zero, FormFactor 0.90–1.00 — textbook CD8 T cell morphology and
  marker profile.
- The same absolute-gate approach (pre-ratio-fix) **failed** to detect any
  of these 3 known-real CD8 T cells, rejecting each one on the `CD4 <=
  0.025` gate despite correct CD3/CD8 signal, because CD4 sat at
  ~0.029–0.040 in this acquisition batch (bleedthrough/baseline-shift, not
  biology).

This validation pair is the evidence that the ratio method isn't producing
artifactual CD8 counts — keep it in mind if a reviewer asks why ratio
gating was used instead of absolute thresholds.

## 5. Tile overlap and deduplication

- Tiles generated at 4000×4000 px with 100 px overlap between adjacent
  tiles (avoids truncating objects — especially ramified glia — that sit
  on a naive tile boundary).
- `analysis/merge_brain.py` deduplicates by computing each tile's "inner
  zone" (the sub-region of the tile that does *not* overlap any neighbor —
  50 px stripped from each interior edge, full extent kept on
  image-boundary edges) and keeping only objects whose global centroid
  falls inside their own tile's inner zone.
- Typical dedup drop rate observed: 3–7% of objects per class (larger/more
  numerous objects like astrocytes trend toward the higher end since their
  centroids are more likely to fall in overlap-adjacent regions).
- See `tiling/README.md` for the full manifest schema and the requirement
  that any reimplemented tiling step preserve this overlap.
