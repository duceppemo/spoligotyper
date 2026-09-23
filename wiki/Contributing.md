# Contributing

Bug reports, questions and pull requests are welcome. See
[CONTRIBUTING.md](https://github.com/duceppemo/spoligotyper/blob/main/CONTRIBUTING.md) for how to report a bug,
set up a development environment, open a pull request and make a release, and the
[code of conduct](https://github.com/duceppemo/spoligotyper/blob/main/CODE_OF_CONDUCT.md).

## Continuous integration
* Lint with ruff, and tests on Python 3.10, 3.12 and 3.14, for every push and pull request.
* The end-to-end tests run Seal on synthetic genomes and reads with known spoligotypes.
* Coverage is uploaded to Codecov from the Python 3.12 job with the `CODECOV_TOKEN` secret.
* Dependabot opens one pull request per month to update the GitHub Actions used by the workflows.

## Editing this wiki
This wiki is maintained in the [`wiki/`](https://github.com/duceppemo/spoligotyper/tree/main/wiki) folder of the
main repository and published to the GitHub wiki automatically when changes are pushed to `main`.
**Do not edit the wiki on GitHub directly**: changes would be overwritten. Edit the files in `wiki/` and open a pull
request instead.

* Page file names become page titles: `Output-files.md` → "Output files".
* Link to other pages without the extension: `[Usage](Usage)`.
* `_Sidebar.md` is the navigation shown on every page.
* Images are stored in the repository's `assets/` folder and linked with their
  `https://raw.githubusercontent.com/duceppemo/spoligotyper/main/assets/...` URL.

## Logo
`assets/logo_source.svg` is the original logo. Its text needs the Sora font, so `assets/make_logo.py` converts it to
outlines and writes `assets/logo.svg`, `assets/logo_dark.svg` (for dark backgrounds) and the bitmap used in the PDF
report, `spoligotyper/data/logo.png`. See the script for its requirements.
