# Validation

`run_validation.sh` checks spoligotyper against:
* 16 public reference genomes (`genomes.tsv`): 13 *M. tuberculosis* complex genomes of known lineage (lineages 1 to
  4 and 6, *M. bovis*, BCG, *M. canettii*) and 3 non-tuberculous mycobacteria;
* public Illumina reads of *M. bovis* AF2122/97 (ERR1744454);
* simulated reads: pure samples at 10x and 30x, single-end and paired-end, a mix of two strains (70% H37Rv, 30%
  *M. bovis*), and a culture contaminated with *M. marinum* (equal depth, so 40% of the reads are MTBC: its genome is larger).

```
bash run_validation.sh 8   # 8 threads
```
It downloads about 600 MB into `data/`, writes the reports in `results/`, and prints `results/validation.md`. The
results of the latest release are on the [Validation](https://github.com/duceppemo/spoligotyper/wiki/Validation)
page of the wiki.
