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
Sample     SpacerCount                                                                                                                  Binary                                       Octal            Hexadecimal        Spoligotype
AF2122_97  56:47:0:58:63:0:76:0:0:0:0:0:50:43:38:0:36:29:45:50:47:51:50:54:140:74:73:78:57:56:50:46:39:41:45:45:54:53:0:0:0:0:0  1101101000001110111111111111111111111100000  664073777777600  6D-03-5F-7F-FF-60  SB0140
12:00:04 INFO    AF2122_97: SB0140 (octal 664073777777600)
12:00:04 INFO    Report saved in results/AF2122_97_spoligotyping.txt (4.2 s)
```
The sample is named after the file. Each spacer is either found in about 30 to 80 reads (present), or in none
(absent): a clear-cut result. Spacer 25 has about twice as many reads as the others because it is present twice in
the DR locus of AF2122/97. With paired-end reads, add `-r2 R2.fastq.gz`.

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
(SIT) in SITVIT (see the [FAQ](FAQ#my-sample-is-spoligo-not-found)).

## 4. Combine the reports
Each sample has its own report. To make one table, with the header once:
```
awk 'FNR == 1 && NR > 1 {next} 1' results/*_spoligotyping.txt > all_samples.tsv
```

## Next steps
* [Output files](Output-files): what each column means
* [How it works](How-it-works): how spacers are detected, and when to change `--min-count`
* [Usage](Usage): all the options
