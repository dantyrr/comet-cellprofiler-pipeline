# Container

`cellprofiler.def` is a Singularity/Apptainer definition file that builds
(or you can just pull) the official CellProfiler 4.2.8 Docker image.

```bash
# build from the .def file
singularity build cellprofiler.sif cellprofiler.def

# or skip building and just pull the base image directly
singularity pull cellprofiler.sif docker://cellprofiler/cellprofiler:4.2.8
```

The resulting `cellprofiler.sif` is a multi-GB binary and is intentionally
**not** committed to this repo (see `.gitignore`). Build or pull it once per
machine/cluster and place it at `$COMET_PROJECT_ROOT/cellprofiler.sif`.

## Version pinning — important

This pipeline was validated specifically on **CellProfiler 4.2.8**. Other
versions are **not guaranteed** to reproduce identical segmentation results —
adaptive thresholding behavior (particularly Minimum Cross-Entropy, used
throughout this pipeline) has changed across CellProfiler versions
historically. Do not swap in a newer/older CellProfiler image without
re-validating segmentation output against known test images.
