# pomRelate

Browser-based gene exploration across seven model organisms, centered on *Schizosaccharomyces pombe*. The static interface combines NOG-based cross-species matching, STRING interaction networks, GO/KEGG enrichment and precomputed eggNOG gene trees.

## Use

Serve the repository locally:

```bash
python -m http.server 8000
```

Open `http://localhost:8000`, select the source species, enter gene names or protein IDs, optionally choose target species, set the PPI threshold and select **Analyze**. The default PPI threshold is 700. Tables and figures export as CSV, PNG, SVG or PDF; gene trees also export as Newick, NEXUS or PhyloXML.

## Data and interpretation

Supports *S. pombe*, *S. cerevisiae*, *D. melanogaster*, *C. elegans*, *A. thaliana*, mouse and human. NOG assignments use the Eukaryota level.

The application uses bundled STRING v12.0 data, GO/KEGG annotations and eggNOG v7 trees. NOG-based matching falls back to aliases where assignments are missing; alias matches can be incomplete or misleading and require review. Shared group membership should be interpreted in its orthology context.

GO/KEGG over-representation uses hypergeometric tests and Benjamini-Hochberg adjustment against species-specific annotation backgrounds. Enriched-term clustering uses Jaccard overlap and UPGMA; it is distinct from the precomputed gene phylogenies. Review identifiers, background coverage and tree-match warnings before biological interpretation.

The data are static and require explicit updates. KEGG coverage varies across species. Large bundled files can make the clone and browser loading substantial; a local static server is required for the documented workflow.

## References

- Szklarczyk, D., et al. (2023). The STRING database in 2023. *Nucleic Acids Research*, *51*(D1), D483–D489.
- Kanehisa, M., et al. (2023). KEGG for taxonomy-based analysis. *Nucleic Acids Research*, *51*(D1), D587–D592.
- The Gene Ontology Consortium. (2023). The GO knowledgebase in 2023. *Genetics*, *224*(1), iyad031.
- Hernández-Plaza, A., et al. (2026). eggNOG v7. *Nucleic Acids Research*, *54*(D1), D402.
