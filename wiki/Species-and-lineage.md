# Species and lineage

Besides the spoligotype, spoligotyper checks every sample for:
* **the species**, from regions of difference (RD) deleted in some members of the complex;
* **the lineage**, from the SNP barcode of [Coll *et al.* 2014](https://doi.org/10.1038/ncomms5812);
* **how much of the sample is MTBC**, to detect contaminated samples;
* **mixed samples**, containing several strains.

It takes a second pass over the reads (a few seconds). `--no-species` skips it.

## Species: regions of difference
spoligotyper measures the read depth of the five regions of difference of the classical RD PCR scheme, relative to
MTBC-specific control regions:

| Region | Deleted in |
|---|---|
| RD1 | BCG and the Dassie bacillus |
| RD4 | *M. bovis* and BCG (and some *M. canettii*) |
| RD7 | Lineage 6 (*M. africanum* West African 2) and all the animal-adapted lineages |
| RD9 | *M. africanum* (lineages 5 and 6) and all the animal-adapted lineages |
| RD12 | *M. bovis*, BCG, *M. caprae* and *M. orygis* (and some *M. canettii*) |

A region is **deleted** when its depth is at most 10% of the control depth, **present** when it is at least 50%, and
**partial** in between (a mixed sample, or a region only partly deleted).

| RD1 | RD4 | RD7 | RD9 | RD12 | Species |
|---|---|---|---|---|---|
| + | + | + | + | + | *M. tuberculosis*; *M. canettii* when there is no standard spacer and no specific lineage SNP |
| + | + | + | − | + | *M. africanum* (lineage 5, West African 1) |
| + | + | − | − | + | *M. africanum* (lineage 6, West African 2), *M. microti*, *M. pinnipedii* or *M. mungi*: lineage 6 with its SNP, the others with the BOV_AFRI SNP only |
| + | − or + | + | + or − | − or + | *M. canettii*, when RD7 is present but RD4 or RD12 is deleted |
| + | + | − | − | − | *M. orygis* or *M. caprae* |
| + | − | − | − | − | *M. bovis* |
| − | − | − | − | − | *M. bovis* BCG |
| − | + | − | − | + | Dassie bacillus |

Any other combination is reported as an unusual RD profile. Notes:
* In some RD tables, lineage 5 and lineage 6 are called *M. africanum* "1b" and "1a": lineage 5 keeps RD7, lineage 6
  (like *M. africanum* GM041182) lost it, like the animal lineages.
* *M. microti* has its own deletion near RD1 (RD1<sup>mic</sup>), but it does not include the RD1 segments used here:
  *M. microti* is RD1 +, as in the RD PCR scheme.
* *M. canettii* is diverse: of four *M. canettii* genomes, one lacks RD12, one lacks RD4, and two have the RD profile
  of *M. tuberculosis*. None has any of the 43 standard spacers, which identifies them. Some carry the lineage 4 SNP
  (the only lineage defined by the H37Rv allele): the lineage SNPs are not reliable for *M. canettii*.
* *M. caprae* and *M. orygis* carry the SNP of the *M. bovis* clade (BOV) of the barcode.
* No genome of the Dassie bacillus is publicly available: its row follows the RD PCR scheme and was not validated.

The RD segments and the control regions are 100 bp pieces of the H37Rv genome, chosen by
[`scripts/make_reference_data.py`](https://github.com/duceppemo/spoligotyper/blob/main/scripts/make_reference_data.py)
from public genomes: each is found in every MTBC genome that has the region, and in none of the genomes lacking it
(AF2122/97, BCG Pasteur, *M. africanum* GM041182, *M. canettii* CIPT 140010059) or of six non-tuberculous
mycobacteria (*M. marinum*, *M. kansasii*, *M. avium*, *M. abscessus*, *M. ulcerans*, *M. smegmatis*). They were then
checked on 15 other genomes (*M. caprae*, *M. orygis*, *M. microti*, *M. pinnipedii*, *M. mungi*, *M. africanum*,
*M. canettii*): see [Validation](Validation).

If the control regions are not found (fewer than 3 reads each, or half of them missing), the sample is reported as
**MTBC not detected**, and no species or lineage is called.

## Lineage: SNP barcode
The barcode of Coll *et al.* has 62 SNPs, each specific to a lineage or sublineage: lineages 1 to 7 and their
sublineages (e.g. 4.3.4.2, LAM), the *M. bovis* clade (BOV, also carried by *M. caprae* and *M. orygis*), and the clade
of the animal lineages and lineage 6 (BOV_AFRI). For each SNP,
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
spoligotyper warns when the evidence disagrees:
* RD9 intact but lineage SNPs of an RD9-deleted lineage (5, 6, animal lineages), or the reverse;
* lineage SNPs of the *M. bovis* clade with RD4 present (except with the *M. caprae* / *M. orygis* profile);
* RD4 deleted but lineage SNPs of a lineage other than *M. bovis*;
* RD7 present with lineage 6 or animal-lineage SNPs, or deleted with lineage 5 SNPs;
* lineage SNPs in a sample identified as *M. canettii*;
* an SB number (defined for the RD9-deleted lineages only) for a sample with RD9 intact.

## Validation
See [Validation](Validation): 31 reference genomes of known species and lineage, and pure, mixed and contaminated
read sets.
