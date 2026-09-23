# Spatial analysis

Every cell carries `global_x` / `global_y` — whole-slide pixel coordinates, valid
across tile boundaries because `merge_brain.py` deduplicates on the global
centroid. No extra processing is needed to do spatial work.

Results are in [`results_8-28-26.md`](results_8-28-26.md) §3. This document covers
the method and its pitfalls.

## Running it

```bash
python3 analysis/spatial_neighborhood.py DATA_DIR \
        --radius-um 25 --split-by CD68 --px-um 0.28 --out per_microglion.csv
```

For each microglion it counts cells of each type within the radius, then:

1. compares neighbourhood composition between treated and control within route, and
2. with `--split-by MARKER`, splits microglia at **each brain's median** for that
   marker and compares the two halves.

`--split-by` accepts any marker column. The output CSV is one row per microglion,
so any further analysis can be done from it directly.

## Why microglia-centred

Microglia are the numerous population (13k–22k per brain). The same question asked
from the T-cell side rests on 9–14 CD8 cells per brain and is not powered — except
in the specific form "for each T cell, what fraction of nearby microglia are
CD68-high", which *is* well powered because each T cell contributes a fraction
rather than a count.

## Four pitfalls, all of which we hit

**1. Compare composition, not counts.** Total neighbours per microglion rose
~1.20× in both routes between control and treated. That uniform density increase
inflates every cell type roughly equally and will masquerade as "more astrocytes,
more neurons, more microglia near treated microglia". Normalise to percentages.

**2. Prefer within-brain comparisons.** At n=2 per group, between-group
comparisons are dominated by animal-to-animal variation — for most quantities here
the within-group spread exceeds the between-group difference. A median split
*within* each brain makes each brain its own control and cancels density, batch
and section quality at once. This is why the CD68 result replicates 8/8 while
treated-vs-control comparisons do not.

**3. Sweep the radius.** A genuine local association decays monotonically toward
1.0 as the radius grows (T cells: 4.94× at 15 µm → 2.31× at 200 µm). A regional
confound — both cell types concentrating in the same anatomical region — gives a
roughly flat ratio instead. The decay is the main evidence that the effect is
local rather than regional.

**4. Watch for circularity.** Splitting microglia on **GFAP** gives a 12×
astrocyte-neighbour ratio, which is meaningless: a microglial object measured as
GFAP-high is one with astrocyte signal bleeding into its mask. Any split marker
that is primarily expressed by the *neighbouring* cell type will do this. Iba1 is
the useful control in the other direction — it is the microglial identity marker,
so it should show no association, and it doesn't (0.90×, 5/8 brains).

## Pixel size

Distances assume **0.28 µm/px**, inferred from a COMET export (10.46 mm over
37,341 px) rather than read from the file. Confirm against OME `PhysicalSizeX`
before reporting distances:

```python
import tifffile, re
with tifffile.TiffFile(path) as t: xml = t.ome_metadata
re.search(r'PhysicalSizeX="([^"]+)"', xml).group(1)
```

Ratios and p-values are unaffected; only the stated radii scale.

## What is not available

There is **no anatomical annotation** — coordinates are raw pixels with no
cortex/hippocampus/white-matter labels. To ask "where in the section", either draw
ROIs (QuPath or ImageJ) and do point-in-polygon against these coordinates, or
register to the Allen atlas. The first is straightforward; the second is a
project.

## Natural extensions

- Continuous marker level rather than a median split, to test whether the
  association is graded
- A permutation null (shuffle cell-type labels within a brain) to put a p-value on
  colocalisation directly
- Neighbourhood composition around astrocytes or neurons, not only microglia
- Distance to a manually drawn landmark — ventricle, lesion, midline
