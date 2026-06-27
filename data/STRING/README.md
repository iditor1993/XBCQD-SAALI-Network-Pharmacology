# STRING Database Files

This directory contains STRING database files for human protein-protein interactions (PPI).

## Files

### Included (committed to repository)

- `9606.protein.info.v12.0.txt` (6.0 MB)
  - Protein annotation information for *Homo sapiens*

### Excluded (must download separately)

- `9606.protein.links.v12.0.txt` (~602 MB)
  - Protein-protein interaction network for *Homo sapiens*
  - **Excluded from GitHub due to file size limit (100MB)**

## Download Instructions

1. Visit: https://string-db.org/cgi/download?species_text=Homo+sapiens
2. Download: `9606.protein.links.v12.0.txt` (under "Protein links")
3. Place the downloaded file in this directory: `data/STRING/9606.protein.links.v12.0.txt`

## Version

- STRING v12.0
- Species: Homo sapiens (Taxonomy ID: 9606)
- Last accessed: 2025-04

## Reference

Szklarczyk D, et al. (2023). "The STRING database in 2023: protein-protein association networks and functional characterization of user-uploaded gene/measurement sets." *Nucleic Acids Research*, 51(D1): D638-D646.
