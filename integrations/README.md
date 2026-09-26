# Workflow integrations

Ready-to-submit wrappers for workflow systems. They use the bioconda package (and its BioContainers image), so they
can only be submitted once the matching spoligotyper version is on bioconda.

Submission branches are ready in the maintainer's forks: `spoligotyper` in
[duceppemo/modules](https://github.com/duceppemo/modules/tree/spoligotyper) (nf-core) and
[duceppemo/tools-iuc](https://github.com/duceppemo/tools-iuc/tree/spoligotyper) (Galaxy). Keep them in sync with
this folder.

## nf-core module: `nf-core/modules/nf-core/spoligotyper/`
Files in the layout of the [nf-core/modules](https://github.com/nf-core/modules) repository. To submit, copy the
`spoligotyper/` folder to `modules/nf-core/` of a clone of nf-core/modules, then from its root:
```
nf-core modules lint spoligotyper
nf-core modules test spoligotyper   # Creates the snapshot file tests/main.nf.test.snap
```
and open a pull request. Update the version in `main.nf` (container) and `environment.yml` for each release.

The SIT database is not included: to fill the SIT columns, download it with `spoligotyper-download-sit` and give it
with `--sit-db` in `task.ext.args` (the module would need an extra input for the file).

Outputs (`<prefix>` is the sample id, or `task.ext.prefix`): `<prefix>_spoligotyping.txt` (table),
`<prefix>_spoligotyping.json`, `<prefix>_spoligotyping.pdf` (optional: not written with `--no-pdf` in
`task.ext.args`) and `<prefix>_spoligotyping_mqc.json` for MultiQC. Reads can be single-end, paired-end or an assembly.

## Galaxy tool: `galaxy/spoligotyper.xml`
Checked with `planemo lint` and `planemo shed_lint` (with `.shed.yml`), and its command was tested by rendering the template for single-end, paired-end and
fasta inputs (`planemo test` needs the bioconda package: run it before submitting). To publish it, open a pull request to
[galaxyproject/tools-iuc](https://github.com/galaxyproject/tools-iuc) (folder `tools/spoligotyper/`, with
`test-data/` and a `.shed.yml`), or upload it to the Galaxy Tool Shed. Update `@TOOL_VERSION@` for each release.
