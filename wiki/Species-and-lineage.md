# Species and lineage

Besides the spoligotype, spoligotyper checks every sample for:
* **the species**, from regions of difference (RD) deleted in some members of the complex;
* **the lineage**, from the SNP barcode of [Coll *et al.* 2014](https://doi.org/10.1038/ncomms5812);
* **how much of the sample is MTBC**, to detect contaminated samples;
* **mixed samples**, containing several strains.

It takes a second pass over the reads (a few seconds). `--no-species` skips it.

## Species: regions of difference
spoligotyper determines the five regions of difference of the classical RD PCR scheme in silico, from the sequencing
data. Each region is the part of the H37Rv genome missing from *M. bovis* AF2122/97 (RD4, RD7, RD9, RD12) or from
BCG Pasteur (RD1), and is represented by 100 bp segments inside it:

| Region | H37Rv region (NC_000962.3) | Segments |
|---|---|---|
| RD1 | 4,350,251-4,359,740 (9,490 bp) | 20 |
| RD4 | 1,696,001-1,708,740 (12,740 bp) | 20 |
| RD7 | 2,208,001-2,220,740 (12,740 bp) | 20 |
| RD9 | 2,330,051-2,332,140 (2,090 bp) | 8 |
| RD12 | 3,485,101-3,487,540 (2,440 bp) | 20 |

A segment is **found** when Seal finds it in the reads at a depth of at least 5% of the depth of the MTBC control
regions, or in the assembly. A region is:
* **present** when all its segments are found (in reads, up to 10% of them, at least 1, may be missing: low depth);
* **deleted** when at most 10% of its segments (at least 1: stray reads) are found;
* **partial** (partially deleted) in between: the report gives the H37Rv coordinates of the missing segments;
* **reduced** (reads only) when its segments are found, but at less than 50% of the control depth: a mix of strains
  with and without the region?

The + and − of the RD profile below mean that the DNA of the region is present or absent, like a PCR with primers
inside the region (amplification = present). spoligotyper does not measure amplicon sizes or deletion junctions:
assays that distinguish RDs by the size of an amplicon spanning the region may report partial or strain-specific
deletions differently. A partially deleted region counts as + when at least half of its segments are found, − otherwise.

The regions are usually deleted in:

| Region | Deleted in |
|---|---|
| RD1 | BCG and the Dassie bacillus |
| RD4 | *M. bovis* and BCG (and some *M. canettii*) |
| RD7 | Lineage 6 (*M. africanum* West African 2) and all the animal-adapted lineages |
| RD9 | *M. africanum* (lineages 5 and 6) and all the animal-adapted lineages |
| RD12 | *M. bovis*, BCG, *M. caprae* and *M. orygis* (and some *M. canettii*) |

| RD1 | RD4 | RD7 | RD9 | RD12 | Species |
|---|---|---|---|---|---|
| + | + | + | + | + | *M. tuberculosis*; *M. canettii* when there is no standard spacer and no specific lineage SNP |
| + | + | + | − | + | *M. africanum* (lineage 5, West African 1) |
| + | + | − | − | + | *M. africanum* (lineage 6, West African 2), *M. microti*, *M. pinnipedii* or *M. mungi*: lineage 6 with its SNP; *M. microti* when RD1 is partially deleted by RD1<sup>mic</sup>; otherwise the three, with the BOV_AFRI SNP |
| + | − or + | + | + or − | − or + | *M. canettii*, when RD7 is present but RD4 or RD12 is deleted |
| + | + | − | − | − | *M. orygis* or *M. caprae* |
| + | − | − | − | − | *M. bovis* |
| − | − | − | − | − | *M. bovis* BCG |
| − | + | − | − | + | Dassie bacillus |

Any other combination is reported as an unusual RD profile. Notes:
* In some RD tables, lineage 5 and lineage 6 are called *M. africanum* "1b" and "1a": lineage 5 keeps RD7, lineage 6
  (like *M. africanum* GM041182) lost it, like the animal lineages.
* *M. microti* has its own deletion, RD1<sup>mic</sup>, which removes the part of RD1 from Rv3871 to Rv3876: RD1 is
  reported as partially deleted (11 of 20 segments found, H37Rv 4,350,651-4,354,450 missing, in three *M. microti*
  genomes) and counts as + in the RD profile, as *M. microti* is RD1 + in the RD PCR scheme. This RD1<sup>mic</sup>
  pattern identifies *M. microti* among the three species of its row. One error is tolerated on each
  side (a few stray reads on one RD1<sup>mic</sup> segment, or no read on one GC-rich segment outside it), as seen in
  real Illumina reads of *M. microti* (see [Validation](Validation)).
* Other partial deletions seen in the validation genomes: RD4 of *M. canettii* ET1291 (9 of 20 segments), RD1 of the
  *M. mungi* draft genome (17 of 20) and RD7 of *M. africanum* RB30065 (18 of 20). They are reported with their
  coordinates and do not change the species.
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
sequenced. A pure culture gives about 100% with simulated reads, and 69% to 100% with the real reads of pure cultures
of the [Validation](Validation) (duplicate reads and uneven coverage lower it); a sample with 50% of reads from another
organism gives about 50%. spoligotyper warns below 60%. The estimate assumes a 4.4 Mb genome.

## Mixed samples
A mix of strains is flagged when:
* both alleles of a lineage SNP are each carried by at least 15% of 10 or more reads: the lineage is then reported as,
  e.g., `mixed: 4 68%, BOV 35%` (percentage of reads with each lineage's allele), and the species as
  "MTBC, mixed sample?";
* SNPs of incompatible lineages are present (e.g. lineage 2 and lineage 4);
* a region of difference is present at reduced depth (all its segments found, at less than half the depth of the
  control regions);
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
