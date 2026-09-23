# Contributing to spoligotyper

Thanks for your interest! Bug reports, questions, ideas and pull requests are all welcome.
Please follow the [code of conduct](CODE_OF_CONDUCT.md).

## Reporting a bug or asking a question
[Open an issue](https://github.com/duceppemo/spoligotyper/issues/new/choose) and pick the matching template.
For bugs, please include the version (`spoligotyper --version`), the command you ran and the messages it printed
(ideally from a run with `-v`). For a spoligotype that looks wrong, include the report line with its `SpacerCount`
column. Check the [troubleshooting page](https://github.com/duceppemo/spoligotyper/wiki/Troubleshooting) and the
[FAQ](https://github.com/duceppemo/spoligotyper/wiki/FAQ) first.

## Development setup
```
git clone https://github.com/duceppemo/spoligotyper
cd spoligotyper
conda env create -f environment.yml
conda activate spoligotyper
pip install --no-deps -e .
pytest --cov
```
The end-to-end tests run Seal on small synthetic genomes and reads, and are skipped if `seal.sh` is not installed.
Lint with [ruff](https://docs.astral.sh/ruff/): `ruff check .`, or install the git hook with `pre-commit install`.
The tests run on GitHub Actions for Python 3.10, 3.12 and 3.14.

To check the results on real data, run the [example](examples/run_example.sh): `bash examples/run_example.sh`.

## Reference data and validation
* `spoligotyper/data/` holds the reference data; its README gives the source of each file.
* `python scripts/make_reference_data.py` rebuilds the species markers and lineage SNP sequences from public genomes
  (needs BBTools and internet access); the result is deterministic.
* `bash validation/run_validation.sh` checks the results on reference genomes and read sets of known species,
  lineage and spoligotype; update the [Validation](https://github.com/duceppemo/spoligotyper/wiki/Validation) page
  with its report when the results change.
* `integrations/` holds the nf-core module and the Galaxy tool: update their version for each release.

## Pull requests
* Branch from `main` and keep each pull request focused on one change.
* Add or update tests for any change in behaviour.
* Match the style of the surrounding code.
* Update the documentation in `wiki/` if the change affects users (see below), and add a line to
  `wiki/Changelog.md` under "Unreleased".

## Documentation
The [wiki](https://github.com/duceppemo/spoligotyper/wiki) is maintained in the [`wiki/`](wiki/) folder of this
repository and published automatically when changes reach `main`. **Do not edit the wiki on GitHub directly**:
those changes would be overwritten. Edit the files in `wiki/` in your pull request instead.

## Releases (maintainers)
1. Update the version in `spoligotyper/__init__.py` and `CITATION.cff` (`version` and `date-released`).
2. Rename "Unreleased" in `wiki/Changelog.md` to the version and date.
3. Commit, tag (`git tag -a vX.Y.Z`), push the commit and the tag, then create the GitHub release once the tests pass.
4. Publishing the GitHub release automatically uploads the package to
   [PyPI](https://pypi.org/project/spoligotyper/) (`.github/workflows/publish.yml`, trusted publishing: no token
   needed; the release tag must match the package version).
5. Bioconda: update `recipe/meta.yaml` (version and sha256 of the new GitHub tarball). Bioconda's bot usually opens
   the pull request to `bioconda-recipes` on its own for new releases; otherwise copy the recipe and open it by hand.
