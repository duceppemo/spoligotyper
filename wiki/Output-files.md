# Output files

## `<sample>_spoligotyping.txt`
A tab-separated table with a header and one line per sample. The same table is printed to the screen.

| Column | Example | Description |
|---|---|---|
| `Sample` | `AF2122_97` | Sample name |
| `SpacerCount` | `56:47:0:58:...` | Number of reads containing each spacer, spacers 1 to 43, separated by `:`. For assemblies, the number of contigs |
| `Binary` | `1101101000001...` | 43 digits: 1 = spacer present (count ≥ `--min-count`), 0 = absent |
| `Octal` | `664073777777600` | 15-digit octal code |
| `Hexadecimal` | `6D-03-5F-7F-FF-60` | Hexadecimal code, 6 blocks |
| `Spoligotype` | `SB0140` | SB number of the pattern in the [Mbovis.org](https://www.mbovis.org/) database, or `Spoligo not found` |

### Octal code
The binary pattern is cut into 14 groups of 3 spacers, and each group is written as one octal digit (000 = 0,
001 = 1, ..., 111 = 7). Spacer 43 is the 15th digit, 0 or 1. It is the standard code used by SITVIT and in
publications, and it can be converted back to the binary pattern.

### Hexadecimal code
The binary pattern is cut into 6 blocks of 7, 7, 7, 7, 8 and 7 spacers, and each block is written as a 2-digit
hexadecimal number.

### Reading the counts
`SpacerCount` shows how confident each call is:
* **Reads**: present spacers usually have tens of reads and absent spacers 0. Counts close to `--min-count`
  deserve a second look, and spoligotyper warns about absent spacers seen in 1 to 4 reads. See
  [How it works](How-it-works#minimum-count).
* **Paired-end reads**: when a spacer is found in one read of a pair, both reads are counted, so counts are about
  twice the number of DNA fragments containing the spacer.
* **Assemblies**: counts are 0 or 1 (occasionally 2, if a spacer is split over two contigs or repeated).
