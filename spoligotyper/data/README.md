# Reference data

| File | Content | Source |
|---|---|---|
| `spoligo_spacers.fasta` | The 43 spacer sequences (25 bp) of the standard spoligotyping membrane | Kamerbeek J *et al.* J Clin Microbiol 35:907-914 (1997). https://doi.org/10.1128/jcm.35.4.907-914.1997 |
| `spoligotype_db.txt` | 1,976 SB numbers with their octal and binary patterns | Snapshot of the [Mbovis.org](https://www.mbovis.org/) database. Smith NH, Upton P. Infect Genet Evol 12:873-876 (2012). https://doi.org/10.1016/j.meegid.2011.08.002 |
| `markers.fasta` | 100 bp segments of the H37Rv genome (NC_000962.3): MTBC-specific controls, and segments of the regions of difference RD1, RD4, RD7, RD9 and RD12 | Built by `scripts/make_reference_data.py` from public genomes (see the script). Regions of difference: Brosch R *et al.* PNAS 99:3684-3689 (2002). https://doi.org/10.1073/pnas.052548299 |
| `lineage_barcode.tsv` | The 62 SNPs of the lineage barcode, lineage names and main spoligotype families | Coll F *et al.* A robust SNP barcode for typing *Mycobacterium tuberculosis* complex strains. Nat Commun 5:4812 (2014). https://doi.org/10.1038/ncomms5812. Supplementary Table 3 and Supplementary Data 1, reused under the [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/) license |
| `lineage_snps.fasta` | For each barcode SNP, 61 bp of H37Rv centred on the SNP, with each allele | Built by `scripts/make_reference_data.py` from `lineage_barcode.tsv` and H37Rv |
| `logo.png` | Logo, for the PDF report | Built by `assets/make_logo.py` |
