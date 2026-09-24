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

### `Several files give the same sample name`
In batch mode (`-i`), each sample name must come from a single fasta file, a single fastq file, or an R1/R2 pair.
The message lists the conflicting files: rename, move or remove them. See
[Usage](Usage#a-folder-of-samples-batch-mode) for how sample names are derived from file names.

### Some samples `failed`
In batch mode, a sample that cannot be typed (truncated file, not a sequence file, ...) is reported as `failed` with
its error in the table and the PDF report, and the other samples are still typed. spoligotyper exits with code 1.

### `mixed sample?`
The lineage SNPs, a region of difference or the spacer counts suggest several strains. See
[Species and lineage](Species-and-lineage#mixed-samples). Check the sample (e.g. single colony re-culture) before
reporting its spoligotype: the pattern of a mixed sample is the union of those of its strains.

### `too little MTBC DNA to check the species`
Spacers were found but the MTBC control regions have fewer than 3 reads each: the depth is too low, or most reads
come from something else. The species and lineage are not called.

### `No spacer found`
None of the 43 spacers is in the file. Check that the sample is from the *M. tuberculosis* complex (e.g. with
[Kraken2](https://github.com/DerrickWood/kraken2) or [mashID](https://github.com/duceppemo/mashID)), and that the
file is not empty or truncated.

### `spacer(s) called absent were seen in fewer than 5 reads`
Some spacers were found, but in too few reads to be called present. This happens with low coverage data, or with
mixed or contaminated samples. See [How it works](How-it-works#minimum-count) to decide whether to lower
`--min-count`.

### `assembly typed with minimum count N: spacers are probably missed`
Assemblies contain each spacer once. Do not set `--min-count` for assemblies: it is 1 by default.

### `typed as reads in fasta format`
The fasta file has more than 1,000 sequences and more than 13 Mb (3 times the genome): it holds reads (e.g. from
`fasterq-dump --fasta`), not an assembly, and is typed as reads (minimum count 5, depth and MTBC fraction estimated).
If it really is an assembly (e.g. a metagenome assembly), use `--min-count 1`.

### Paths with spaces, commas, "xmx" or "xms"
They are supported. Seal cannot read them (BBTools splits arguments on spaces and commas, and reads any argument
containing "xmx" or "xms" as a Java memory setting), so spoligotyper gives Seal links with neutral names in a
temporary folder. The same links give Seal the right extension for files named without one (e.g. Galaxy's `.dat`
files) or gzipped without `.gz`. Only the system temporary folder itself must be free of these (set `TMPDIR`
otherwise).

### `-r2 is only for paired-end fastq files`
`-r2` must be the R2 reads of a paired-end run. Assemblies are given with `-r1` alone.

### `ModuleNotFoundError: No module named 'pkg_resources'`
That is version 0.1, which does not work with recent versions of setuptools. Update to version 0.2.0 or later
(see [Installation](Installation#updating-from-01)).
