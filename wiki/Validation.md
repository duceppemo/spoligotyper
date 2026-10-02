# Validation

spoligotyper was checked on public reference genomes and on read sets of known species, lineage and
spoligotype. Everything below is reproducible with
[`validation/run_validation.sh`](https://github.com/duceppemo/spoligotyper/tree/main/validation), which downloads
the data, simulates the read sets (BBTools `randomreads.sh`, fixed seeds) and checks each result.

**Result: 59 of 59 checks passed.**

## Reference genomes
* *M. tuberculosis*, *M. bovis*, BCG, *M. africanum* GM041182: the expected lineages are the predictions of
  [Coll *et al.* 2014](https://doi.org/10.1038/ncomms5812) (Supplementary Table 4) for the same genomes;
  spoligotyper may report a more specific sublineage on the same branch (EAI5: 1.1.2 within 1.1). Documented
  spoligotypes: H37Rv 777777477760771 (SIT451), AF2122/97 SB0140, BCG SB0120, and the Beijing signature (spacers
  1-34 absent, 35-43 present) for CCDC5079.
* Other members of the complex (*M. caprae*, *M. orygis*, *M. microti*, *M. pinnipedii*, *M. mungi*, *M. africanum*,
  *M. canettii*): the expected species is the NCBI taxonomy of the genome, with the RD profile of the RD PCR scheme
  (see [Species and lineage](Species-and-lineage)), and the lineage the clade of the barcode (BOV: *M. bovis* clade;
  BOV_AFRI: animal lineages and lineage 6). These genomes were not used to choose the RD segments: they are an
  independent test. The *M. africanum* RB30001 genome has the spoligotype of GM041182 (AFRI_1, lineage 6); the
  lineage of RB30065 is not documented.
* La lineage: the lineage of the livestock-associated complex
  ([Zwyer *et al.* 2021](https://doi.org/10.12688/openreseurope.14029.2)): AF2122/97 belongs to clonal complex Eu1
  (La1.8.1), BCG to La1.2, *M. caprae* to La2 and *M. orygis* to La3; the other genomes have none.
* The RD profile is shown after the species: RD1, RD4, RD7, RD9 and RD12, + present, - deleted, p partially
  deleted (see [Species and lineage](Species-and-lineage#species-regions-of-difference)), r present at reduced depth.
* SIT and SITVIT2 family: from the SITVIT2 patterns published with SpolLineages (`spoligotyper-download-sit`).
  The documented SITs are checked: H37Rv SIT451, BCG SIT482, and SIT1 for the Beijing strain CCDC5079.

| Sample | Organism | SB | SIT (family) | Octal | Species (RD1 RD4 RD7 RD9 RD12) | Lineage | Expected lineage | La lineage | Result |
|---|---|---|---|---|---|---|---|---|---|
| H37Rv | M. tuberculosis H37Rv | Not in Mbovis.org | SIT451 (T-H37Rv) | 777777477760771 | M. tuberculosis (+ + + + +) | 4.9 | 4.9 | - | OK |
| CDC1551 | M. tuberculosis CDC1551 | Not in Mbovis.org | SIT549 (X3) | 700076757760771 | M. tuberculosis (+ + + + +) | 4.1.1.3 | 4.1.1.3 | - | OK |
| Erdman | M. tuberculosis Erdman ATCC 35801 | Not in Mbovis.org | SIT1230 (H1) | 777757774020771 | M. tuberculosis (+ + + + +) | 4.1.2.1 | 4.1.2.1 | - | OK |
| F11 | M. tuberculosis F11 | Not in Mbovis.org | SIT33 (LAM3) | 776177607760771 | M. tuberculosis (+ + + + +) | 4.3.2.1 | 4.3.2.1 | - | OK |
| KZN1435 | M. tuberculosis KZN 1435 | Not in Mbovis.org | SIT60 (LAM4) | 777777607760731 | M. tuberculosis (+ + + + +) | 4.3.3 | 4.3.3 | - | OK |
| CCDC5079 | M. tuberculosis CCDC5079 (Beijing) | Not in Mbovis.org | SIT1 (Beijing) | 000000000003771 | M. tuberculosis (+ + + + +) | 2.2.1 | 2.2.1 | - | OK |
| CAS_NITR204 | M. tuberculosis CAS/NITR204 | Not in Mbovis.org | Not in SITVIT2 list | 677777441741771 | M. tuberculosis (+ + + + +) | 3 | 3 | - | OK |
| EAI5_NITR206 | M. tuberculosis EAI5/NITR206 | Not in Mbovis.org | Not in SITVIT2 list | 667777467740071 | M. tuberculosis (+ + + + +) | 1.1.2 | 1.1 | - | OK |
| RGTB423 | M. tuberculosis RGTB423 | Not in Mbovis.org | Not in SITVIT2 list | 777736033740711 | M. tuberculosis (+ + + + +) | 1.2.2 | 1.2.2 | - | OK |
| GM041182 | M. africanum GM041182 (lineage 6) | SB0147 | SIT181 (AFRI_1) | 770777777777671 | M. africanum (lineage 6) (+ + - - +) | 6 | 6 | - | OK |
| AF2122_97 | M. bovis AF2122/97 | SB0140 | SIT683 (BOV_2) | 664073777777600 | M. bovis (+ - - - -) | BOV | BOV | La1.8.1 | OK |
| BCG_Pasteur | M. bovis BCG Pasteur 1173P2 | SB0120 | SIT482 (BOV_1) | 676773777777600 | M. bovis BCG (- - - - -) | BOV | BOV | La1.2 | OK |
| M_canettii | M. canettii CIPT 140010059 | SB2277 | SIT2669 (ATYPIC) | 000000000000000 | M. canettii (+ + + + -) | - | - | - | OK |
| M_marinum | M. marinum M | SB2277 | SIT2669 (ATYPIC) | 000000000000000 | MTBC not detected (? ? ? ? ?) | - | - | - | OK |
| M_kansasii | M. kansasii ATCC 12478 | SB2277 | SIT2669 (ATYPIC) | 000000000000000 | MTBC not detected (? ? ? ? ?) | - | - | - | OK |
| M_avium | M. avium 104 | SB2277 | SIT2669 (ATYPIC) | 000000000000000 | MTBC not detected (? ? ? ? ?) | - | - | - | OK |
| M_africanum_RB30001 | M. africanum RB30001 (same spoligotype as GM041182, AFRI_1) | SB0147 | SIT181 (AFRI_1) | 770777777777671 | M. africanum (lineage 6) (+ + - - +) | 6 | 6 | - | OK |
| M_africanum_RB30065 | M. africanum RB30065 (lineage not documented) | Not in Mbovis.org | Orphan (AFRI_2) | 474077607177071 | M. africanum (lineage 5) (+ + p - +) | 5 | not documented | - | OK |
| M_caprae_Allgaeu | M. caprae Allgaeu | SB0418 | SIT647 (BOV_4-CAPRAE) | 200003777377600 | M. caprae (+ + - - -) | BOV | BOV | La2 | OK |
| M_caprae_SY-1 | M. caprae SY-1 | SB0418 | SIT647 (BOV_4-CAPRAE) | 200003777377600 | M. caprae (+ + - - -) | BOV | BOV | La2 | OK |
| M_orygis_51145 | M. orygis 51145 | Not in Mbovis.org | Not in SITVIT2 list | 600000000000271 | M. orygis (+ + - - -) | BOV | BOV | La3 | OK |
| M_orygis_NIAB | M. orygis NIAB_BDWBCSHFL_1 | Not in Mbovis.org | Not in SITVIT2 list | 600740007774671 | M. orygis (+ + - - -) | BOV | BOV | La3 | OK |
| M_microti_OV254 | M. microti OV254 | SB0118 | SIT539 (microti) | 000000000000600 | M. microti (p + - - +) | BOV_AFRI | BOV_AFRI | - | OK |
| M_microti_MausIV | M. microti Maus IV | SB0118 | SIT539 (microti) | 000000000000600 | M. microti (p + - - +) | BOV_AFRI | BOV_AFRI | - | OK |
| M_microti_94-2272 | M. microti 94-2272 | SB0118 | SIT539 (microti) | 000000000000600 | M. microti (p + - - +) | BOV_AFRI | BOV_AFRI | - | OK |
| M_pinnipedii_MP1 | M. pinnipedii MP1 (draft) | SB0155 | SIT593 (PINI1) | 074000037777600 | M. microti, M. pinnipedii or M. mungi (+ + - - +) | BOV_AFRI | BOV_AFRI | - | OK |
| M_pinnipedii_BAA-688 | M. pinnipedii ATCC BAA-688 (draft) | Not in Mbovis.org | Not in SITVIT2 list | 074000033747400 | M. microti, M. pinnipedii or M. mungi (+ + - - +) | BOV_AFRI | BOV_AFRI | - | OK |
| M_mungi_BM22813 | M. mungi BM22813 (draft) | SB1960 | SIT3151 (mungi) | 672600000000671 | M. microti, M. pinnipedii or M. mungi (p + - - +) | BOV_AFRI | BOV_AFRI | - | OK |
| M_canettii_ET1291 | M. canettii ET1291 | SB2277 | SIT2669 (ATYPIC) | 000000000000000 | M. canettii (+ p + + +) | - | - | - | OK |
| M_canettii_CIPT140070010 | M. canettii CIPT 140070010 | SB2277 | SIT2669 (ATYPIC) | 000000000000000 | M. canettii (+ + + + +) | - | - | - | OK |
| M_canettii_CIPT140070017 | M. canettii CIPT 140070017 | SB2277 | SIT2669 (ATYPIC) | 000000000000000 | M. canettii (+ + + + +) | 4 | not documented | - | OK |

The non-tuberculous mycobacteria and the *M. canettii* genomes have none of the 43 standard spacers: their pattern
is SB2277, the pattern with no spacer, flagged with a warning (see the [FAQ](FAQ#a-sample-with-no-spacer-is-sb2277)).
*M. canettii* CIPT 140070017 carries the lineage 4 SNP (the H37Rv allele); it is identified as *M. canettii* from its
missing spacers, with a warning that the lineage SNPs are not reliable for *M. canettii*.

## Reads
Public reads (ENA), of strains of known spoligotype, species and lineage, most of them the strains of the reference
genomes above:
* **ERR1744454**: *M. bovis* AF2122/97, Illumina single-end.
* **SRR12006063**: *M. tuberculosis* H37Rv, Illumina HiSeq 4000 paired-end. This lab stock is not quite the reference
  genome: in most reads, spacer 39 is followed directly by the direct repeat and the sequence that follows spacer 43
  in the reference, so most cells lost spacers 40 to 43, and a minority kept spacers 41 to 43 (recombination between
  direct repeats during passage). Another H37Rv run (ERR15989112) also lacks spacer 40. The resulting pattern,
  777777477760731, is SIT1647 of the T-H37Rv family in SITVIT2. spoligotyper flags the minority as a mixed sample.
* **ERR027297**: *M. microti* Maus IV, Illumina GAII paired-end (2010, first 1.5 million pairs). With the strong GC
  bias of these reads, one RD1<sup>mic</sup> segment gets a few stray reads and one GC-rich RD1 segment none:
  *M. microti* is still recognised from its RD1<sup>mic</sup> deletion, which tolerates one error on each side.
* **SRR16643349**: *M. orygis* 51145, Illumina MiniSeq paired-end. 8 reads (median of the present spacers: 73) carry
  a variant of spacer 3 with one SNP, flanked by direct repeats, absent from the PacBio assembly of the strain: a
  minority population or cross-contamination. Spacer 3 is called present (octal 700000000000271, instead of
  600000000000271 for the assembly), and flagged; the octal code is not checked.
* **ERR2383628**: *M. africanum* RB30001 (lineage 6), Illumina HiSeq 2500 paired-end (first million pairs).
* **SRR18636082**, **SRR23035463**: *M. canettii* ET1291, Illumina NextSeq paired-end, and nanopore (first 30,000
  reads, mean length 4.3 kb): the nanopore reads get the long-read warning.

Simulated reads:
* **sim_...**: 150 bp reads with sequencing errors, simulated from the reference genomes above at 10x or 30x.
* **sim_mixed_H37Rv70_AF2122_30**: 70% H37Rv and 30% *M. bovis* reads.
* **sim_contaminated_H37Rv15x_marinum15x**: H37Rv and *M. marinum* reads at the same depth. The *M. marinum* genome
  is larger (6.6 Mb), so 40% of the reads are MTBC: the estimated MTBC fraction, 0.40, is right.

| Sample | SB | Octal | Species | Lineage | La lineage | MTBC fraction | Warnings | Result |
|---|---|---|---|---|---|---|---|---|
| ERR1744454 | SB0140 | 664073777777600 | M. bovis | BOV | La1.8.1 | 1.00 | - | OK |
| SRR12006063 | Not in Mbovis.org | 777777477760731 | M. tuberculosis | 4.9 | - | 0.69 | spacers with few reads, mixed sample | OK |
| ERR027297 | SB0118 | 000000000000600 | M. microti | BOV_AFRI | - | 0.70 | - | OK |
| SRR16643349 | Not in Mbovis.org | 700000000000271 | M. orygis | BOV | La3 | 0.87 | mixed sample | OK |
| ERR2383628 | SB0147 | 770777777777671 | M. africanum (lineage 6) | 6 | - | 1.00 | - | OK |
| SRR18636082 | SB2277 | 000000000000000 | M. canettii | - | - | 0.86 | RD4 partially deleted, no spacer (M. canettii) | OK |
| SRR23035463 | SB2277 | 000000000000000 | M. canettii | - | - | 0.95 | long reads, RD4 partially deleted, no spacer (M. canettii) | OK |
| sim_H37Rv_30x | Not in Mbovis.org | 777777477760771 | M. tuberculosis | 4.9 | - | 1.00 | - | OK |
| sim_H37Rv_30x_PE | Not in Mbovis.org | 777777477760771 | M. tuberculosis | 4.9 | - | 0.97 | - | OK |
| sim_H37Rv_10x | Not in Mbovis.org | 777677475760761 | M. tuberculosis | 4.9 | - | 1.00 | low depth, spacers with few reads | OK |
| sim_Beijing_30x | Not in Mbovis.org | 000000000003771 | M. tuberculosis | 2.2.1 | - | 1.00 | - | OK |
| sim_AF2122_97_30x | SB0140 | 664073777777600 | M. bovis | BOV | La1.8.1 | 1.00 | - | OK |
| sim_mixed_H37Rv70_AF2122_30 | Not in Mbovis.org | 777777777767771 | MTBC, mixed sample? | mixed: 4 68%, BOV 35%, BOV_AFRI 26% | mixed: La1/La2 26%, La1 28%, La1.8 31%, La1.8.1 28% | 1.00 | mixed sample, mixed sample, spacers with few reads | OK |
| sim_contaminated_H37Rv15x_marinum15x | Not in Mbovis.org | 777777477760771 | M. tuberculosis | 4.9 | - | 0.40 | contamination, low MTBC depth | OK |

* At 10x, three present spacers get fewer than 5 reads and are called absent: the pattern is wrong, but flagged
  ("low depth", "spacers with few reads"). Use `--min-count 2` or `3` for low depth data, see
  [How it works](How-it-works#minimum-count).
* The mixed sample is detected from both alleles of the lineage SNPs; its spoligotype is the union of the two
  strains' patterns.
* The contaminated sample keeps its spoligotype, species and lineage, with a contamination warning.

## Livestock lineages
One run per lineage and sublineage of the livestock-associated complex, from the genomes of
[Zwyer *et al.* 2021](https://doi.org/10.12688/openreseurope.14029.2) (extended data Table 1), with the sublineage
and the SB number given in the paper (first 600,000 read pairs; all the reads for the three runs whose reads are
sorted by position: their first reads cover only part of the genome). The check is the livestock lineage and the
species; "SNPs" gives, for each group called, the SNPs with the derived allele out of those covered.
spoligotyper finds the SB number of the paper for 13 of the 14 runs. The *M. orygis* run (ERR2659154) differs from
the paper's SB0422 at spacer 3 only: 140 of its reads carry a variant of spacer 3 with 2 mismatches
(CCGTGCTTCCAGTGATC**A**CCTT**G**TA), which spoligotyper, allowing 1 mismatch per spacer, does not count, while the
paper's in silico typing (KvarQ) does. The *M. orygis* NIAB genome above also differs from SB0422 at spacer 3 only.
SRR1173570 (from a chimpanzee) has only about 50% of MTBC reads, which spoligotyper flags.

| Run | Origin | Expected | La lineage | SNPs | Species | SB (paper) | SB | Result |
|---|---|---|---|---|---|---|---|---|
| SRR8065072 | Uganda, cattle | La1.1 | La1.1 | La1 4/4, La1.1 4/4 | M. bovis | SB1405 | SB1405 | OK |
| SRR7851302 | France, cattle | La1.2 | La1.2 | La1 4/4, La1.2 4/4 | M. bovis | SB0134 | SB0134 | OK |
| SRR1173570 | Uganda, chimpanzee | La1.3 | La1.3 | La1 4/4, La1.3 4/4 | M. bovis | SB0133 | SB0133 | OK |
| ERR229952 | Russia, - | La1.4 | La1.4 | La1 4/4, La1.4 5/5 | M. bovis | SB1919 | SB1919 | OK |
| ERR2212116 | Germany, red deer | La1.5 | La1.5 | La1 4/4, La1.5 5/5 | M. bovis | SB0989 | SB0989 | OK |
| ERR1203064 | Ghana, human | La1.6 | La1.6 | La1 4/4, La1.6 4/4 | M. bovis | SB0944 | SB0944 | OK |
| ERR550954 | Germany, - | La1.7.1 | La1.7.1 | La1 4/4, La1.7 4/4, La1.7.1 5/5 | M. bovis | SB0339 | SB0339 | OK |
| ERR2815520 | Germany, cattle | La1.7.X | La1.7.X | La1 4/4, La1.7 4/4, La1.7.X-unk4 5/5 | M. bovis | SB0120 | SB0120 | OK |
| SRR8064859 | Zambia, cattle | La1.7.X | La1.7.X | La1 4/4, La1.7 4/4, La1.7.X-unk5 5/5 | M. bovis | SB0120 | SB0120 | OK |
| ERR841810 | United Kingdom, cattle | La1.8.1 | La1.8.1 | La1 4/4, La1.8 4/4, La1.8.1 5/5 | M. bovis | SB0140 | SB0140 | OK |
| SRR7851304 | France, cattle | La1.8.2 | La1.8.2 | La1 4/4, La1.8 4/4, La1.8.2 5/5 | M. bovis | SB0134 | SB0134 | OK |
| SRR1792479 | USA, cervid | La1.8.X | La1.8.X | La1 4/4, La1.8 4/4, La1.8.X-unk6 5/5 | M. bovis | SB1069 | SB1069 | OK |
| SRR13888769 | Spain, goat | La2 | La2 | La2 5/5 | M. caprae | SB0415 | SB0415 | OK |
| ERR2659154 | Australia, human | La3 | La3 | La3 5/5 | M. orygis | SB0422 | Not in Mbovis.org | OK |

## Comparison with SpoTyping
[SpoTyping](https://github.com/xiaeryu/SpoTyping-v2.0) 2.1 (Xia *et al.* 2016, BLAST-based), an independent in silico
spoligotyping tool, was run on the 13 original MTBC genomes (`--seq`): the spoligotypes are identical for 12 of
13 genomes.

| Sample | spoligotyper | SpoTyping | Difference |
|---|---|---|---|
| H37Rv | 777777477760771 | 777777477760771 | identical |
| CDC1551 | 700076757760771 | 700076757760771 | identical |
| Erdman | 777757774020771 | 777757774020771 | identical |
| F11 | 776177607760771 | 776177607760771 | identical |
| KZN1435 | 777777607760731 | 777777607760731 | identical |
| CCDC5079 | 000000000003771 | 000000000003771 | identical |
| CAS_NITR204 | 677777441741771 | 677777441741771 | identical |
| EAI5_NITR206 | 667777467740071 | 667777477740671 | spacers 24, 37, 38 |
| RGTB423 | 777736033740711 | 777736033740711 | identical |
| GM041182 | 770777777777671 | 770777777777671 | identical |
| AF2122_97 | 664073777777600 | 664073777777600 | identical |
| BCG_Pasteur | 676773777777600 | 676773777777600 | identical |
| M_canettii | 000000000000000 | 000000000000000 | identical |

SpoTyping's EAI5/NITR206 pattern has spacers 24, 37 and 38 present. In this assembly, the closest sequences to these
spacers have 6, 4 and 5 mismatches out of 25 bp, so spoligotyper, which allows 1 mismatch, calls them absent.
