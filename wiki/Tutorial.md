# Tutorial

This tutorial types three public samples with known spoligotypes. It takes a few minutes, most of it to download
210 MB of reads.

| Sample | Data | Expected spoligotype |
|---|---|---|
| `AF2122_97` | *M. bovis* AF2122/97, Illumina single-end reads ([ERR1744454](https://www.ebi.ac.uk/ena/browser/view/ERR1744454)) | SB0140 |
| `NC_002945` | *M. bovis* AF2122/97 reference genome ([NC_002945.4](https://www.ncbi.nlm.nih.gov/nuccore/NC_002945.4)) | SB0140 |
| `H37Rv` | *M. tuberculosis* H37Rv reference genome ([NC_000962.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3)) | octal 777777477760771, no SB number |

## Run it all at once
The repository includes a script that downloads the data, types the three samples and checks the results:
```
git clone https://github.com/duceppemo/spoligotyper
bash spoligotyper/examples/run_example.sh 4   # 4 threads
```
It ends with `OK: results match expected_results.tsv`. The rest of this page does the same steps by hand.

## 1. Download the data
```
mkdir -p tutorial/data && cd tutorial
curl -L -o data/AF2122_97.fastq.gz \
    https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR174/004/ERR1744454/ERR1744454.fastq.gz
efetch='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id='
curl -L -o data/NC_002945.fasta "${efetch}NC_002945.4"
curl -L -o data/H37Rv.fasta "${efetch}NC_000962.3"
```

## 2. Reads
```
spoligotyper -r1 data/AF2122_97.fastq.gz -o results/
```
```
12:00:00 INFO    Spoligotyping AF2122_97 (fastq, minimum count 5)
12:00:07 INFO    AF2122_97: SB0140 (octal 664073777777600), M. bovis, lineage BOV
Sample     SpacerCount           Binary                                       Octal            ...  Spoligotype  ...  Species   Lineage  ...
AF2122_97  56:47:0:58:63:0:...   1101101000001110111111111111111111111100000  664073777777600  ...  SB0140       ...  M. bovis  BOV      ...
12:00:07 INFO    Report saved in results/AF2122_97_spoligotyping.txt
12:00:07 INFO    Report saved in results/AF2122_97_spoligotyping.json
12:00:07 INFO    Report saved in results/AF2122_97_spoligotyping.pdf
12:00:07 INFO    Done in 7.8 s
```
The main columns of the table:

| Column | Value |
|---|---|
| SpacerCount | 56:47:0:58:63:0:76:0:0:0:0:0:50:43:38:0:36:29:45:50:47:51:50:54:140:74:73:78:57:56:50:46:39:41:45:45:54:53:0:0:0:0:0 |
| Octal | 664073777777600 |
| Spoligotype | SB0140 |
| Depth | 65 |
| Species | M. bovis |
| RD9, RD4, RD1 | deleted, deleted, present |
| Lineage | BOV |
| MTBCFraction | 1.00 |

The sample is named after the file, and the results are saved as a table and a PDF report. Each spacer is either
found in about 30 to 80 reads (present), or in none (absent): a clear-cut result. The depth, 65x, is estimated from
the number of bases, and all the reads appear to be MTBC (`MTBCFraction` 1.00). The species, *M. bovis*, comes from
the regions of difference RD9 and RD4, both deleted, and the lineage from the SNP barcode: see
[Species and lineage](Species-and-lineage). Spacer 25 has about twice as many reads as the others because it is
present twice in the DR locus of AF2122/97. With paired-end reads, add `-r2 R2.fastq.gz`.

The pattern, 1101101000001110111111111111111111111100000, is SB0140, the spoligotype of AF2122/97.
It shows the classic *M. bovis* signature: spacers 3, 9, 16 and 39 to 43 are absent.

## 3. Assemblies
```
spoligotyper -r1 data/NC_002945.fasta -o results/
spoligotyper -r1 data/H37Rv.fasta -o results/
```
For fasta files, each spacer is found once or not at all, so the minimum count is automatically 1. The AF2122/97
assembly gives the same spoligotype as its reads, SB0140.

H37Rv gives `Spoligo not found`: SB numbers only exist for the animal-adapted lineages. Its octal code,
777777477760771, is the one to use for *M. tuberculosis*, for example to look up its shared international type
(SIT) in SITVIT (see the [FAQ](FAQ#my-sample-is-spoligo-not-found)). The species is *M. tuberculosis* (RD9, RD4 and
RD1 present) and the lineage 4.9, "Euro-American (H37Rv-like)".

## 4. All at once, with a PDF report
Type every sample of the `data/` folder in one run:
```
spoligotyper -i data/ -o batch/ --operator "Your Name"
```
The reads and the two assemblies are detected automatically. `batch/spoligotyping.tsv` has one line per sample, and
`batch/spoligotyping_report.pdf` has a summary page, one section per sample with the reads supporting each spacer,
and the run information (software versions, checksums, parameters) for quality assurance. See
[Output files](Output-files#pdf-report-spoligotyping_reportpdf-or-sample_spoligotypingpdf), and the
[report of this tutorial](https://github.com/duceppemo/spoligotyper/blob/main/assets/example_report.pdf) as an example.

## Next steps
* [Output files](Output-files): what each column and each part of the PDF report means
* [How it works](How-it-works): how spacers are detected, and when to change `--min-count`
* [Usage](Usage): all the options
