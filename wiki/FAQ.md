# FAQ

### My sample is `Spoligo not found`
The pattern is not in the [Mbovis.org](https://www.mbovis.org/) database. Either:
* **It is a human-adapted lineage** (*M. tuberculosis*, *M. africanum*, ...): SB numbers only exist for the
  animal-adapted lineages (*M. bovis*, *M. caprae*, *M. pinnipedii*, *M. microti*, ...). Use the octal code, for
  example to find the shared international type (SIT) and lineage in the
  SITVIT database (H37Rv, 777777477760771, is SIT451).
* **It is a new pattern**: new *M. bovis* patterns can be submitted to Mbovis.org to get an SB number.
* **A spacer is miscalled**: look at the `SpacerCount` column for counts close to `--min-count`, and see
  [How it works](How-it-works#minimum-count).

### Can I use another database?
Yes, with `--db my_database.txt`. The file has one pattern per line, with 3 columns separated by spaces or tabs:
the octal code, the name, and the binary pattern. Lines starting with `#` are ignored. For example, with SIT numbers:
```
# octal          name    binary
777777477760771  SIT451  1111111111111111111001111111111100001111111
```
spoligotyper checks that each octal code matches its binary pattern.

### How recent is the included database?
It is a snapshot of the Mbovis.org database with 1,976 SB patterns. To use a newer version, download it from
Mbovis.org, format it as above and use `--db`.

### Can I type nanopore reads?
Yes. Seal allows 1 mismatch per 25 bp spacer, which suits accurate reads (recent flow cells and basecallers).
With less accurate reads, fewer reads match each spacer: check the `SpacerCount` column, and lower `--min-count` if
present spacers have low counts.

### Does spoligotyper work on *M. tuberculosis*?
Yes: spacer detection and the binary, octal and hexadecimal codes work for the whole complex. Only the SB number
is specific to the animal-adapted lineages (see above).

### Why are there 43 spacers?
The 43 spacers are those of the standard spoligotyping membrane
([Kamerbeek *et al.* 1997](https://doi.org/10.1128/jcm.35.4.907-914.1997)), so that *in silico* results can be
compared with laboratory spoligotyping. Other spacers exist in some strains but are not part of the standard
pattern.

### Should I type reads or the assembly?
Reads, when you have them: the DR locus is repetitive and can be broken or collapsed in short-read assemblies.
See [How it works](How-it-works#limitations).

### How do I cite spoligotyper?
See [Home](Home#citing).
