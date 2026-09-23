# spoligotyper

**spoligotyper** performs *in silico* spoligotyping of *Mycobacterium tuberculosis* complex (MTBC) samples from
sequencing reads (fastq) or genome assemblies (fasta). It finds the 43 spacers of the direct repeat (DR) locus with
[Seal](https://sourceforge.net/projects/bbmap/) from BBTools, and reports the spoligotype as binary, octal and
hexadecimal codes, and as an SB number from the [Mbovis.org](https://www.mbovis.org/) database.

It works for all members of the complex (*M. tuberculosis*, *M. bovis*, *M. caprae*, *M. pinnipedii*,
*M. microti*, *M. africanum*, *M. canettii*, ...). SB numbers are only defined for the animal-adapted lineages
(RD9-deleted: *M. bovis*, *M. caprae*, ...); for human-adapted lineages, use the octal code (see the [FAQ](FAQ)).

```
$ spoligotyper -r1 AF2122_97.fastq.gz -o results/
Sample     SpacerCount           Binary                                       Octal            Hexadecimal        Spoligotype
AF2122_97  56:47:0:58:63:0:...   1101101000001110111111111111111111111100000  664073777777600  6D-03-5F-7F-FF-60  SB0140
```

## Getting started
* [Installation](Installation)
* [Tutorial](Tutorial): three public genomes, step by step
* [Usage](Usage): input files, options and examples

## Understanding the results
* [Output files](Output-files): the report, column by column
* [How it works](How-it-works): spoligotyping, spacer detection, and the choice of the minimum count

## Reference
* [Troubleshooting](Troubleshooting)
* [FAQ](FAQ)
* [Changelog](Changelog)
* [Contributing](Contributing)

## Citing
Please cite spoligotyper (GitHub's "Cite this repository" button on the
[repository page](https://github.com/duceppemo/spoligotyper)) and BBTools: Bushnell B. BBTools.
https://sourceforge.net/projects/bbmap/
