# ASODesigner

![Python](https://img.shields.io/badge/python-3.9--3.12-blue.svg)
![Status](https://img.shields.io/badge/status-experimental-orange.svg)
![License](https://img.shields.io/badge/license-CC%20BY%204.0-lightgrey.svg)

> Feature extraction and analysis utilities for antisense oligonucleotide (ASO) design, built for the TAU-Israel 2025 iGEM project.

## Features

- MOE and LNA candidate ranking pipeline with optional feature breakdowns and off-target scoring.
- Modular feature calculators covering GC metrics, hybridization, RNA accessibility, folding, and toxicity heuristics.
- One-step asset bootstrap that downloads the human GFF database and Bowtie index structure required by the pipeline.
- Ready to embed in FastAPI backends or standalone discovery notebooks.

## Installation

### From PyPI

```bash
pip install asodesigner
```

## Required Assets

The generator expects the human annotation database and Bowtie index structure to exist under `/tmp/.cache/asodesigner`. Download them once via:

```python
from asodesigner.download_assets import ensure_assets

# Validates and downloads if necessary the needed files
ensure_assets() 
```

The helper skips files that already exist and only downloads missing assets.

## Quick Start

Generate top ASO candidates for a human gene, complete with feature annotations:

```python
from asodesigner.aso_generator import design_asos

# Retrieve the top 3 MOE + LNA designs for MALAT1
candidates = design_asos(
    organismName="human",
    geneName="MALAT1",
    geneData=None,
    top_k=3,
    includeFeatureBreakdown=True,
)

print(candidates[["Sequence", "mod_pattern"]])
```

- Set `geneData` to a custom transcript sequence to work outside the reference genome.
- With `includeFeatureBreakdown=True`, additional columns (e.g., `exp_ps_hybr`, `gc_content`, `at_skew`, `off_target`, `on_target`) are attached to each row.
- For lower-level feature utilities, explore modules under `src/asodesigner/`.


## Development Workflow

1. Update or add functionality under `src/asodesigner/`.
2. Keep imports relative within the package (for example, `from .util import helper`).
3. Run `pytest` (or `python -m pytest`) to execute the available unit tests.
4. Optionally run `python -m compileall src/asodesigner` to double-check importability before packaging.

## Extending the Project

- **Feature metrics** – Implement additional sequence, structural, or accessibility metrics under `src/asodesigner/features/`. Many modules (e.g., `seq_features.py`, `hybridization.py`) expose template-style functions you can mirror. 
- **Pipeline enrichment** – The cross-chemistry ASO pipeline lives in `src/asodesigner/aso_generator.py`. Add new feature columns inside `add_features_for_output` or extend the returned DataFrame schema to expose your metrics downstream.
- **Constants and configuration** – Global paths and dataset references live in `src/asodesigner/consts.py`. Update these when introducing new organism builds or experimental assets so the rest of the codebase can locate them.
- **Utility helpers** – Shared logic (reverse complement, translation tables, etc.) sits under `src/asodesigner/util.py` and related utilities. Enhance these modules when new workflows require additional helpers.
- **Data workflows** – Reference datasets and caches under `src/data/` pair with the code in `src/asodesigner`. When extending to other organisms or assemblies, follow the existing directory layout so asset downloaders and consts remain consistent.

Have improvements to share? Open an issue or PR—we welcome new metrics, pipeline enrichments, and broader organism support.

## License

Released under an MIT-style license tailored for academic and research use. See `LICENSE` for the complete terms and instructions for commercial enquiries.
