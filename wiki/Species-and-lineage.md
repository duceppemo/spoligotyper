# Species and lineage

Besides the spoligotype, spoligotyper checks every sample for:
* **the species**, from regions of difference (RD) deleted in some members of the complex;
* **the lineage**, from the SNP barcode of [Coll *et al.* 2014](https://doi.org/10.1038/ncomms5812);
* **how much of the sample is MTBC**, to detect contaminated samples;
* **mixed samples**, containing several strains.

It takes a second pass over the reads (a few seconds). `--no-species` skips it.

## Species: regions of difference
spoligotyper measures the read depth of three regions of difference, relative to MTBC-specific control regions:

| Region | Deleted in |
|---|---|
| RD9 | *M. africanum* and all the animal-adapted lineages, including *M. bovis* |
| RD4 | *M. bovis* and BCG only |
| RD1 | BCG (and, with a different deletion, *M. microti*) |

A region is **deleted** when its depth is at most 10% of the control depth, **present** when it is at least 50%, and
**partial** in between (a mixed sample, or a region only partly deleted).

| RD9 | RD4 | RD1 | Species |
|---|---|---|---|
| present | present | present | *M. tuberculosis* (or *M. canettii* when no lineage SNP is found) |
| deleted | present | present | *M. africanum* (with the SNPs of lineage 5 or 6), an animal-adapted lineage that is not *M. bovis* (*M. caprae*, *M. pinnipedii*, ...; BOV_AFRI SNP only), or either when no lineage SNP is found |
| deleted | present | deleted | Animal-adapted lineage with an RD1 deletion, e.g. *M. microti* |
| deleted | deleted | present | *M. bovis* |
| deleted | deleted | deleted | *M. bovis* BCG |

The RD segments and the control regions are 100 bp pieces of the H37Rv genome, chosen by
[`scripts/make_reference_data.py`](https://github.com/duceppemo/spoligotyper/blob/main/scripts/make_reference_data.py)
from public genomes: each is found in every MTBC genome that has the region, and in none of the genomes lacking it
(AF2122/97, BCG Pasteur, *M. africanum* GM041182) or of six non-tuberculous mycobacteria (*M. marinum*,
*M. kansasii*, *M. avium*, *M. abscessus*, *M. ulcerans*, *M. smegmatis*).

If the control regions are not found (fewer than 3 reads each, or half of them missing), the sample is reported as
**MTBC not detected**, and no species or lineage is called.

## Lineage: SNP barcode
The barcode of Coll *et al.* has 62 SNPs, each specific to a lineage or sublineage: lineages 1 to 7 and their
sublineages (e.g. 4.3.4.2, LAM), *M. bovis* (BOV), and the clade of *M. bovis* and lineage 6 (BOV_AFRI). For each SNP,
spoligotyper counts the reads carrying each allele (exact 31-mers). A lineage is called when at least 80% of the reads
covering its SNP (and at least 3 reads, or 1 contig) carry its allele. The most specific lineage is reported, with its
name and its main spoligotype families (e.g. `4.3.4.2`, "Euro-American (LAM)", "LAM11-ZWE, LAM9, LAM1, LAM4").

The spoligotype families (Beijing, EAI, CAS, LAM, T, H, X, ...) are given as the families most often found in each
lineage by Coll *et al.*: spoligotype families are not always monophyletic, while the SNP lineages are.

A few SNP sequences are also found in non-tuberculous mycobacteria. These SNPs are not used to detect mixed samples,
and are ignored when the sample looks contaminated (less than 80% of MTBC reads).

## Contamination: fraction of MTBC reads
For reads, the depth of the MTBC control regions is compared with the depth expected from the number of bases
sequenced. A pure culture gives about 100%; a sample with 50% of reads from another organism gives about 50%.
spoligotyper warns below 60%. The estimate assumes a 4.4 Mb genome; plasmid-rich or very uneven libraries can
lower it slightly.

## Mixed samples
A mix of strains is flagged when:
* both alleles of a lineage SNP are each carried by at least 15% of 10 or more reads: the lineage is then reported as,
  e.g., `mixed: 4 68%, BOV 35%` (percentage of reads with each lineage's allele), and the species as
  "MTBC, mixed sample?";
* SNPs of incompatible lineages are present (e.g. lineage 2 and lineage 4);
* a region of difference is partially deleted;
* some present spacers have much fewer reads than the others (less than 40% of the median, when the median is at
  least 30 reads).

The spoligotype of a mixed sample is the union of the spoligotypes of its strains, and usually matches neither.

## Consistency checks
spoligotyper warns when the evidence disagrees: an RD9 intact but lineage SNPs of an RD9-deleted lineage, lineage SNPs
of *M. bovis* with RD4 present, or an SB number (defined for the RD9-deleted lineages only) for a sample with RD9
intact.

## Validation
See [Validation](Validation): 16 reference genomes of known species and lineage, and pure, mixed and contaminated
read sets.
