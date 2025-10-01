# ASOdesigner

ASOdesigner provides feature extraction utilities that support the design of antisense oligonucleotides (ASOs). The
package bundles folding, accessibility, and sequence analysis helpers that were originally created for the TAU-Israel
2025 iGEM project.

## Installation

You can install the library directly from PyPI once it is published:

```bash
pip install asodesigner
```

For local development and experimentation inside this repository, install the package in editable mode together with
its runtime dependencies:

```bash
pip install -r requirements.txt
```

## Usage

Functions are organised by module. For example, to compute a few sequence-derived features:

```python
from asodesigner.features import seq_features

gc = seq_features.get_gc_content("ATGGCC")
```

Refer to the individual module docstrings for detailed behaviour and expected inputs.

## Development Workflow

1. Update or add functionality under `src/asodesigner/`.
2. Keep imports relative within the package (for example, use `from .util import helper`).
3. Run `python -m compileall src/asodesigner` or your preferred test suite to make sure the code still imports
   correctly.


## Genome Cache

`asodesigner.genome.get_hg38_genome_path()` downloads hg38 into `cache/genomes/hg38`
when the file is first requested and returns the local FASTA path. Set
`ASODESIGNER_HG38_PATH` if you already have a FASTA on disk, or
`ASODESIGNER_HG38_URL` to use a different mirror.

## Docker Compose Workflow

`docker-compose.yml` keeps the genome in a named volume so large files stay out
of the repository. Start the stack with the default curl download:

```bash
docker compose up
```

Switch to the simple Python downloader by setting `GENOME_FETCH_METHOD=python`.
Both services mount the `genome-cache` volume, so subsequent runs reuse the same
FASTA file.

## Building and Publishing to PyPI

The project uses a modern `pyproject.toml` driven build. To prepare a release:

1. **Choose a version** – bump the version field in `pyproject.toml` following semantic versioning.
2. **Build the artifacts**:
   ```bash
   python -m build
   ```
3. **Inspect the build** – verify the generated wheel and source distribution under `dist/`.
4. **Upload to PyPI** (or TestPyPI first) using `twine`:
   ```bash
   python -m twine upload dist/*
   ```

Remember that publishing requires a PyPI account and API token. For TestPyPI, swap the repository URL accordingly.

## License

The code is released under the Creative Commons Attribution 4.0 International License (CC BY 4.0). See `LICENSE` for
full details.
