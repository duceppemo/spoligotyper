# Validation

spoligotyper was checked on public reference genomes and on read sets of known species, lineage and
spoligotype. Everything below is reproducible with
[`validation/run_validation.sh`](https://github.com/duceppemo/spoligotyper/tree/main/validation), which downloads
the data, simulates the read sets (BBTools `randomreads.sh`, fixed seeds) and checks each result.

**Result: 39 of 39 checks passed.**

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
* The RD profile is shown after the species: RD1, RD4, RD7, RD9 and RD12, + present, - deleted.

| Sample | Organism | Spoligotype | Octal | Species (RD1 RD4 RD7 RD9 RD12) | Lineage | Expected lineage | Result |
|---|---|---|---|---|---|---|---|
| H37Rv | M. tuberculosis H37Rv | Spoligo not found | 777777477760771 | M. tuberculosis (+ + + + +) | 4.9 | 4.9 | OK |
| CDC1551 | M. tuberculosis CDC1551 | Spoligo not found | 700076757760771 | M. tuberculosis (+ + + + +) | 4.1.1.3 | 4.1.1.3 | OK |
| Erdman | M. tuberculosis Erdman ATCC 35801 | Spoligo not found | 777757774020771 | M. tuberculosis (+ + + + +) | 4.1.2.1 | 4.1.2.1 | OK |
| F11 | M. tuberculosis F11 | Spoligo not found | 776177607760771 | M. tuberculosis (+ + + + +) | 4.3.2.1 | 4.3.2.1 | OK |
| KZN1435 | M. tuberculosis KZN 1435 | Spoligo not found | 777777607760731 | M. tuberculosis (+ + + + +) | 4.3.3 | 4.3.3 | OK |
| CCDC5079 | M. tuberculosis CCDC5079 (Beijing) | Spoligo not found | 000000000003771 | M. tuberculosis (+ + + + +) | 2.2.1 | 2.2.1 | OK |
| CAS_NITR204 | M. tuberculosis CAS/NITR204 | Spoligo not found | 677777441741771 | M. tuberculosis (+ + + + +) | 3 | 3 | OK |
| EAI5_NITR206 | M. tuberculosis EAI5/NITR206 | Spoligo not found | 667777467740071 | M. tuberculosis (+ + + + +) | 1.1.2 | 1.1 | OK |
| RGTB423 | M. tuberculosis RGTB423 | Spoligo not found | 777736033740711 | M. tuberculosis (+ + + + +) | 1.2.2 | 1.2.2 | OK |
| GM041182 | M. africanum GM041182 (lineage 6) | SB0147 | 770777777777671 | M. africanum (lineage 6) (+ + - - +) | 6 | 6 | OK |
| AF2122_97 | M. bovis AF2122/97 | SB0140 | 664073777777600 | M. bovis (+ - - - -) | BOV | BOV | OK |
| BCG_Pasteur | M. bovis BCG Pasteur 1173P2 | SB0120 | 676773777777600 | M. bovis BCG (- - - - -) | BOV | BOV | OK |
| M_canettii | M. canettii CIPT 140010059 | SB2277 | 000000000000000 | M. canettii (+ + + + -) | - | - | OK |
| M_marinum | M. marinum M | SB2277 | 000000000000000 | MTBC not detected (? ? ? ? ?) | - | - | OK |
| M_kansasii | M. kansasii ATCC 12478 | SB2277 | 000000000000000 | MTBC not detected (? ? ? ? ?) | - | - | OK |
| M_avium | M. avium 104 | SB2277 | 000000000000000 | MTBC not detected (? ? ? ? ?) | - | - | OK |
| M_africanum_RB30001 | M. africanum RB30001 (same spoligotype as GM041182, AFRI_1) | SB0147 | 770777777777671 | M. africanum (lineage 6) (+ + - - +) | 6 | 6 | OK |
| M_africanum_RB30065 | M. africanum RB30065 (lineage not documented) | Spoligo not found | 474077607177071 | M. africanum (lineage 5) (+ + + - +) | 5 | not documented | OK |
| M_caprae_Allgaeu | M. caprae Allgaeu | SB0418 | 200003777377600 | M. orygis or M. caprae (+ + - - -) | BOV | BOV | OK |
| M_caprae_SY-1 | M. caprae SY-1 | SB0418 | 200003777377600 | M. orygis or M. caprae (+ + - - -) | BOV | BOV | OK |
| M_orygis_51145 | M. orygis 51145 | Spoligo not found | 600000000000271 | M. orygis or M. caprae (+ + - - -) | BOV | BOV | OK |
| M_orygis_NIAB | M. orygis NIAB_BDWBCSHFL_1 | Spoligo not found | 600740007774671 | M. orygis or M. caprae (+ + - - -) | BOV | BOV | OK |
| M_microti_OV254 | M. microti OV254 | SB0118 | 000000000000600 | M. microti, M. pinnipedii or M. mungi (+ + - - +) | BOV_AFRI | BOV_AFRI | OK |
| M_microti_MausIV | M. microti Maus IV | SB0118 | 000000000000600 | M. microti, M. pinnipedii or M. mungi (+ + - - +) | BOV_AFRI | BOV_AFRI | OK |
| M_microti_94-2272 | M. microti 94-2272 | SB0118 | 000000000000600 | M. microti, M. pinnipedii or M. mungi (+ + - - +) | BOV_AFRI | BOV_AFRI | OK |
| M_pinnipedii_MP1 | M. pinnipedii MP1 (draft) | SB0155 | 074000037777600 | M. microti, M. pinnipedii or M. mungi (+ + - - +) | BOV_AFRI | BOV_AFRI | OK |
| M_pinnipedii_BAA-688 | M. pinnipedii ATCC BAA-688 (draft) | Spoligo not found | 074000033747400 | M. microti, M. pinnipedii or M. mungi (+ + - - +) | BOV_AFRI | BOV_AFRI | OK |
| M_mungi_BM22813 | M. mungi BM22813 (draft) | SB1960 | 672600000000671 | M. microti, M. pinnipedii or M. mungi (+ + - - +) | BOV_AFRI | BOV_AFRI | OK |
| M_canettii_ET1291 | M. canettii ET1291 | SB2277 | 000000000000000 | M. canettii (+ - + + +) | - | - | OK |
| M_canettii_CIPT140070010 | M. canettii CIPT 140070010 | SB2277 | 000000000000000 | M. canettii (+ + + + +) | - | - | OK |
| M_canettii_CIPT140070017 | M. canettii CIPT 140070017 | SB2277 | 000000000000000 | M. canettii (+ + + + +) | 4 | not documented | OK |

The non-tuberculous mycobacteria and the *M. canettii* genomes have none of the 43 standard spacers: their pattern
is SB2277, the pattern with no spacer, flagged with a warning (see the [FAQ](FAQ#a-sample-with-no-spacer-is-sb2277)).
*M. canettii* CIPT 140070017 carries the lineage 4 SNP (the H37Rv allele); it is identified as *M. canettii* from its
missing spacers, with a warning that the lineage SNPs are not reliable for *M. canettii*.

## Reads
* **ERR1744454**: Illumina reads of *M. bovis* AF2122/97.
* **sim_...**: 150 bp reads with sequencing errors, simulated from the reference genomes above at 10x or 30x.
* **sim_mixed_H37Rv70_AF2122_30**: 70% H37Rv and 30% *M. bovis* reads.
* **sim_contaminated_H37Rv15x_marinum15x**: H37Rv and *M. marinum* reads at the same depth. The *M. marinum* genome
  is larger (6.6 Mb), so 40% of the reads are MTBC: the estimated MTBC fraction, 0.40, is right.

| Sample | Spoligotype | Octal | Species | Lineage | MTBC fraction | Warnings | Result |
|---|---|---|---|---|---|---|---|
| ERR1744454 | SB0140 | 664073777777600 | M. bovis | BOV | 1.00 | - | OK |
| sim_H37Rv_30x | Spoligo not found | 777777477760771 | M. tuberculosis | 4.9 | 1.00 | - | OK |
| sim_H37Rv_30x_PE | Spoligo not found | 777777477760771 | M. tuberculosis | 4.9 | 0.97 | - | OK |
| sim_H37Rv_10x | Spoligo not found | 777677475760761 | M. tuberculosis | 4.9 | 1.00 | low depth, spacers with few reads | OK |
| sim_Beijing_30x | Spoligo not found | 000000000003771 | M. tuberculosis | 2.2.1 | 1.00 | - | OK |
| sim_AF2122_97_30x | SB0140 | 664073777777600 | M. bovis | BOV | 1.00 | - | OK |
| sim_mixed_H37Rv70_AF2122_30 | Spoligo not found | 777777777767771 | MTBC, mixed sample? | mixed: 4 68%, BOV 35%, BOV_AFRI 26% | 1.00 | mixed sample, spacers with few reads | OK |
| sim_contaminated_H37Rv15x_marinum15x | Spoligo not found | 777777477760771 | M. tuberculosis | 4.9 | 0.40 | contamination, low MTBC depth | OK |

* At 10x, three present spacers get fewer than 5 reads and are called absent: the pattern is wrong, but flagged
  ("low depth", "spacers with few reads"). Use `--min-count 2` or `3` for low depth data, see
  [How it works](How-it-works#minimum-count).
* The mixed sample is detected from both alleles of the lineage SNPs; its spoligotype is the union of the two
  strains' patterns.
* The contaminated sample keeps its spoligotype, species and lineage, with a contamination warning.

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
