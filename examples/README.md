# Example

`run_example.sh` downloads three public samples, spoligotypes them, and compares the results with
`expected_results.tsv`:

| Sample | Data | Expected spoligotype |
|---|---|---|
| `AF2122_97` | *M. bovis* AF2122/97, Illumina single-end reads (ENA ERR1744454, 210 MB) | SB0140 |
| `NC_002945` | *M. bovis* AF2122/97 reference genome (NCBI NC_002945.4) | SB0140 |
| `H37Rv` | *M. tuberculosis* H37Rv reference genome (NCBI NC_000962.3) | octal 777777477760771, not an SB pattern |

```
bash run_example.sh 4   # 4 threads
```
Downloads are kept in `data/` for later runs, and the reports are saved in `results/`. The PDF report of this
example is in [`assets/example_report.pdf`](../assets/example_report.pdf). The
[tutorial](https://github.com/duceppemo/spoligotyper/wiki/Tutorial) walks through the same steps by hand.
