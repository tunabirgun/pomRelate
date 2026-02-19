# pomRelate

**A static, client-side bioinformatics tool for cross-species gene analysis in model organisms, centered on *Schizosaccharomyces pombe*.**

pomRelate enables researchers to map orthologs, explore protein-protein interactions, perform functional enrichment analysis, and visualize phylogenetic gene trees across 7 model organisms — entirely in the browser with no backend server required.


## Features

- **Cross-Species Ortholog Mapping** — Identify orthologs across 7 model organisms using phylogeny-based NOG assignments from STRING/eggNOG at the Eukaryota level, with alias-based name matching as fallback
- **Protein-Protein Interaction (PPI) Tables** — Browse interactions from STRING v12.0 with configurable score thresholds (400–999)
- **Interactive PPI Network** — Force-directed network visualization with zoom, pan, and drag. Hub genes identified by degree centrality
- **GO Annotations** — Per-gene Gene Ontology terms (Biological Process, Molecular Function, Cellular Component)
- **KEGG Pathway Annotations** — Per-gene KEGG pathway mappings
- **GO Enrichment Analysis** — Over-representation analysis using Fisher's Exact Test (hypergeometric) with Benjamini-Hochberg FDR correction
- **KEGG Enrichment Analysis** — Pathway enrichment with the same statistical framework
- **Publication-Quality Plots** — Bar charts, dot plots, and hierarchical clustering dendrograms with 7 color palettes
- **Hierarchical Clustering Tree** — UPGMA dendrogram clustering enriched terms by gene set overlap (Jaccard distance)
- **Interactive Gene Tooltips** — Hover over any gene to see annotation, GO terms, KEGG pathways, and links to UniProt/STRING/AmiGO/KEGG/eggNOG
- **Phylogeny Analysis** — Per-gene phylogenetic trees from eggNOG v7, mapped via STRING orthologous groups at the Eukaryota level
- **Phylogenetic Export Formats** — Newick, NEXUS, and PhyloXML
- **Multiple Export Formats** — CSV, PNG, SVG, and PDF
- **Dark / Light Theme** — Persistent theme preference

## Species Coverage

| Species | Taxon ID | KEGG Code | Notes |
|---|---|---|---|
| *Schizosaccharomyces pombe* | 4896 | spo | Default source species |
| *Saccharomyces cerevisiae* | 4932 | sce | |
| *Drosophila melanogaster* | 7227 | dme | |
| *Caenorhabditis elegans* | 6239 | cel | |
| *Arabidopsis thaliana* | 3702 | ath | |
| *Mus musculus* | 10090 | mmu | |
| *Homo sapiens* | 9606 | hsa | |

## Methods

### Orthology Mapping

Cross-species ortholog identification uses a two-tier approach:

1. **Primary — NOG-based orthology**: Each protein's eggNOG orthologous group (NOG) assignment is obtained from STRING v12.0 at the Eukaryota level (taxid 2759). Proteins sharing the same NOG ID are orthologs.
2. **Fallback — Alias matching**: For genes without NOG assignments, cross-species mapping falls back to name/alias matching using STRING alias data.

### PPI Networks

Interaction data is sourced from STRING v12.0. The network visualization uses a synchronous force-directed layout (300 iterations). Hub genes are identified based on degree centrality (top 20%, minimum degree 3).

### Enrichment Analysis

GO and KEGG enrichment is performed using a Fisher's Exact Test (hypergeometric test) with Benjamini-Hochberg FDR correction. Background sets are species-specific genome-wide annotations.

### Phylogeny Analysis

Gene trees are derived from eggNOG v7 pre-computed protein family phylogenies, pruned to model organism species. Orthologous group (NOG) assignments from STRING v12.0 at the Eukaryota level (taxid 2759) link genes to trees.

## Data Sources

| Database | Version | URL |
|---|---|---|
| STRING | v12.0 | https://string-db.org |
| KEGG | Current | https://www.kegg.jp |
| Gene Ontology | via STRING v12.0 | https://geneontology.org |
| eggNOG | v7 | https://eggnogdb.org |

## Usage

1. Select a **source species** (defaults to *S. pombe*)
2. Enter **gene names** (e.g., `cdc2`, `cdc13`, `wee1`, `rad21`)
3. Optionally select **target species** for cross-species ortholog lookup
4. Adjust the **PPI score threshold** (default: 700)
5. Click **Analyze**

## Building Data

The recommended approach rebuilds all species data from scratch and automatically splits any files over 95 MB into chunks for GitHub compatibility:

```bash
python scripts/rebuild_and_chunk.py
```

This downloads per-species data from STRING v12.0 and KEGG, then produces chunked JSON files with a `manifest.json` in each species directory listing the chunk counts. The web app reads the manifest at runtime and reassembles chunks transparently via `loadChunkedJSON()`.

Individual build scripts are also available:

```bash
# Per-species data (STRING + KEGG)
python scripts/build_species_data.py

# Phylogeny data (STRING orthology + eggNOG v7 trees)
python scripts/build_phylogeny_data.py
```

Requirements: Python 3.8+ (standard library only).

## Local Development

```bash
cd pomRelate
python -m http.server 8000
```

## Limitations

- **Orthology coverage** — Genes without STRING/eggNOG Eukaryota-level NOG assignments fall back to alias-based matching
- **Static data** — Pre-downloaded data does not auto-update
- **Large proteomes** — Mouse and human PPI/alias files can be very large; they are automatically chunked into <100 MB pieces for GitHub hosting

## References

- Szklarczyk, D., et al. (2023). The STRING database in 2023. *Nucleic Acids Research*, *51*(D1), D483–D489.
- Kanehisa, M., et al. (2023). KEGG for taxonomy-based analysis. *Nucleic Acids Research*, *51*(D1), D587–D592.
- The Gene Ontology Consortium. (2023). The GO knowledgebase in 2023. *Genetics*, *224*(1), iyad031.
- Hernández-Plaza, A., et al. (2026). eggNOG v7. *Nucleic Acids Research*, *54*(D1), D402.
