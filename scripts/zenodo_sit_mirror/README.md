# Zenodo mirror of the SIT list

`spoligotyper-download-sit` downloads the SITVIT2 pattern list of SpolLineages (GPL-3.0) from GitHub, at the commit
and checksum pinned in `spoligotyper/sitdb.py`. This folder prepares a Zenodo mirror of the same file, used when
GitHub cannot be reached, and gives the list a DOI.

1. `bash prepare.sh sit_mirror` downloads the pinned `Spoligo_list.csv` and its `LICENSE`, checks the checksum and
   writes `README.md` in `sit_mirror/`.
2. On Zenodo, create a new upload with the three files of `sit_mirror/` and the metadata of `zenodo_metadata.json`
   (title, type "Dataset", license GPL-3.0, description, related identifiers), and publish it.
3. Add the file URL of the record (`https://zenodo.org/records/<id>/files/Spoligo_list.csv?download=1`) to `URLS` in
   `spoligotyper/sitdb.py`, after the GitHub URL, and its DOI to the documentation.

The mirror holds the unmodified file under its original license (GPL-3.0), with attribution to its authors. Keep it
unmodified: `spoligotyper-download-sit` checks the same SHA-256 checksum for both sources.
