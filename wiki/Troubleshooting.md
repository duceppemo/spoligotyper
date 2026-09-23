# Troubleshooting

### `"seal.sh" was not found`
`spoligotyper` uses the `seal.sh` found in your `PATH` first, then the one installed next to its own Python
interpreter (the `bin/` folder of the conda environment). Calling `/path/to/envs/spoligotyper/bin/spoligotyper`
without activating the environment therefore works if BBTools is installed in that environment.

Otherwise, install BBTools: `conda install -c bioconda bbmap`.

### `Seal failed`
The message ends with the last lines printed by Seal. Common causes:
* **Java errors** (`java: command not found`, `UnsupportedClassVersionError`): BBTools needs Java. Installing
  BBTools with conda installs a compatible Java.
* **Truncated or corrupted input file**: check it with `gzip -t file.fastq.gz`.

Run with `-v` to see the full Seal command and output.

### Java memory errors
`spoligotyper` gives Seal 1 GB of memory (`--memory 1g`), which is plenty for 43 spacers: Seal's memory does not
grow with the size of the input. If Java cannot reserve it (`Could not reserve enough space for object heap`, e.g.
on a login node or in a small container), lower it: `--memory 500m`.

Before version 0.2.0, Seal sized its memory from the computer's free memory (196 GB on a 512 GB server), which
failed on shared computers and clusters.

### `No spacer found`
None of the 43 spacers is in the file. Check that the sample is from the *M. tuberculosis* complex (e.g. with
[Kraken2](https://github.com/DerrickWood/kraken2) or [mashID](https://github.com/duceppemo/mashID)), and that the
file is not empty or truncated.

### `spacer(s) called absent were seen in fewer than 5 reads`
Some spacers were found, but in too few reads to be called present. This happens with low coverage data, or with
mixed or contaminated samples. See [How it works](How-it-works#minimum-count) to decide whether to lower
`--min-count`.

### `is a fasta file: with --min-count N, spacers are probably missed`
Assemblies contain each spacer once. Do not set `--min-count` for fasta files: it is 1 by default.

### `-r2 is only for paired-end fastq files`
`-r2` must be the R2 reads of a paired-end run. Assemblies are given with `-r1` alone.

### `ModuleNotFoundError: No module named 'pkg_resources'`
That is version 0.1, which does not work with recent versions of setuptools. Update to version 0.2.0 or later
(see [Installation](Installation#updating-from-01)).
