/**
 * pomRelate — Enrichment Analysis Engine
 * Client-side hypergeometric test for GO and KEGG over-representation analysis.
 * No external dependencies.
 */

// ===== Math Utilities =====

// Log-factorial using Stirling's approximation for large n, exact for small n
const _lfCache = [0]; // lfact(0) = 0
function lfact(n) {
    if (n < 0) return 0;
    if (_lfCache[n] !== undefined) return _lfCache[n];
    if (n <= 1) { _lfCache[n] = 0; return 0; }
    // Build up cache iteratively
    let last = _lfCache.length - 1;
    let val = _lfCache[last];
    for (let i = last + 1; i <= n; i++) {
        val += Math.log(i);
        _lfCache[i] = val;
    }
    return _lfCache[n];
}

/**
 * Hypergeometric p-value (upper tail):
 * P(X >= k) where X ~ Hypergeometric(N, K, n)
 *
 * N = total background genes
 * K = genes in background annotated with this term
 * n = size of user gene list (that have annotations)
 * k = genes in user list annotated with this term
 */
function hypergeomPValue(k, n, K, N) {
    if (k <= 0 || n <= 0 || K <= 0 || N <= 0) return 1;
    if (n > N || K > N) return 1;
    if (k > Math.min(n, K)) return 0;

    let pval = 0;
    const minI = Math.max(k, n - (N - K), 0);
    const maxI = Math.min(n, K);
    for (let i = minI; i <= maxI; i++) {
        const lp = lfact(K) - lfact(i) - lfact(K - i)
            + lfact(N - K) - lfact(n - i) - lfact(N - K - n + i)
            - lfact(N) + lfact(n) + lfact(N - n);
        pval += Math.exp(lp);
    }
    return Math.min(pval, 1);
}

/**
 * Benjamini-Hochberg FDR correction.
 * Takes array of {pValue, ...} objects, adds .fdr field in-place.
 * Returns the same array sorted by pValue ascending.
 */
function bhFDR(results) {
    results.sort((a, b) => a.pValue - b.pValue);
    const m = results.length;
    for (let i = m - 1; i >= 0; i--) {
        const rank = i + 1;
        const raw = (results[i].pValue * m) / rank;
        if (i < m - 1) {
            results[i].fdr = Math.min(raw, results[i + 1].fdr);
        } else {
            results[i].fdr = Math.min(raw, 1);
        }
    }
    return results;
}

// ===== GO Enrichment =====

/**
 * Run GO enrichment analysis.
 * @param {string[]} queryProteinIds - resolved protein IDs from user input
 * @param {Object} goData - protein->terms map from species GO data
 * @param {string} [categoryFilter] - optional: "Biological Process", "Molecular Function", "Cellular Component"
 * @returns {Object} { results: [...], stats: { mapped, total, termsTotal } }
 */
function runGOEnrichment(queryProteinIds, goData, categoryFilter) {
    if (!goData || queryProteinIds.length === 0) {
        return { results: [], stats: { mapped: 0, total: 0, termsTotal: 0 } };
    }

    // Background: all proteins with GO annotations
    const bgProteins = Object.keys(goData);
    const N = bgProteins.length;

    // Build term -> set of background proteins
    const termBg = {}; // term -> { desc, category, proteins: Set }
    for (const [pid, terms] of Object.entries(goData)) {
        if (!Array.isArray(terms)) continue;
        for (const t of terms) {
            if (!t || !t.term) continue;
            const category = t.category || 'Unknown';
            if (categoryFilter && !category.includes(categoryFilter)) continue;
            if (!termBg[t.term]) {
                termBg[t.term] = { desc: t.description || '', category, proteins: new Set() };
            }
            termBg[t.term].proteins.add(pid);
        }
    }

    // Query: filter to proteins present in GO background
    const queryInBg = queryProteinIds.filter(pid => goData[pid]);
    const n = queryInBg.length;
    const querySet = new Set(queryInBg);

    if (n === 0) {
        return { results: [], stats: { mapped: 0, total: queryProteinIds.length, termsTotal: Object.keys(termBg).length } };
    }

    // Test each term
    const results = [];
    for (const [term, info] of Object.entries(termBg)) {
        const K = info.proteins.size;
        // Count overlap
        let k = 0;
        const geneHits = [];
        for (const pid of queryInBg) {
            if (info.proteins.has(pid)) {
                k++;
                geneHits.push(pid);
            }
        }
        if (k === 0) continue;

        const pValue = hypergeomPValue(k, n, K, N);
        const expectedK = (K / N) * n;
        const fold = expectedK > 0 ? k / expectedK : 0;

        results.push({
            term,
            description: info.desc,
            category: info.category,
            pValue,
            fdr: 1,
            fold: Math.round(fold * 100) / 100,
            geneCount: k,
            bgCount: K,
            totalGenes: n,
            totalBg: N,
            genes: geneHits,
        });
    }

    bhFDR(results);

    return {
        results,
        stats: { mapped: n, total: queryProteinIds.length, termsTotal: Object.keys(termBg).length },
    };
}

// ===== KEGG Enrichment =====

/**
 * Run KEGG pathway enrichment analysis.
 * Uses KEGG pathway data (from REST API).
 * Resolves STRING protein IDs to KEGG gene names via aliases.
 * @param {string[]} queryProteinIds
 * @param {Object} keggPathwayData - { pathways: {id->name}, gene_pathways: {gene->[pathways]} }
 * @param {Object} aliasData - protein -> [alias strings] from species aliases
 * @param {Object} infoData - protein info for name resolution
 * @returns {Object} { results, stats }
 */
function runKEGGEnrichment(queryProteinIds, keggPathwayData, aliasData, infoData) {
    if (!keggPathwayData || !keggPathwayData.gene_pathways) {
        return { results: [], stats: { mapped: 0, total: queryProteinIds.length, termsTotal: 0 } };
    }

    const genePathways = keggPathwayData.gene_pathways;
    const pathwayNames = keggPathwayData.pathways || {};

    // Build a set of all known KEGG gene names for fast lookup
    const keggGeneSet = new Set(Object.keys(genePathways));

    // Background: all genes with at least one pathway
    const N = keggGeneSet.size;

    // Build pathway -> set of background genes
    const pathwayBg = {};
    for (const [gene, pathways] of Object.entries(genePathways)) {
        if (!Array.isArray(pathways)) continue;
        for (const pw of pathways) {
            if (!pathwayBg[pw]) {
                const lookup = pw.replace(/^path:/, '');
                pathwayBg[pw] = { name: pathwayNames[lookup] || pw, genes: new Set() };
            }
            pathwayBg[pw].genes.add(gene);
        }
    }

    // Resolve query protein IDs to KEGG gene names via aliases
    const queryGeneNames = new Set();
    const pidToKegg = {}; // for reporting back
    for (const pid of queryProteinIds) {
        // Try direct protein ID
        if (keggGeneSet.has(pid)) {
            queryGeneNames.add(pid);
            pidToKegg[pid] = pid;
            continue;
        }
        // Try info preferred name
        if (infoData && infoData[pid]) {
            const name = infoData[pid].name;
            if (name && keggGeneSet.has(name)) {
                queryGeneNames.add(name);
                pidToKegg[pid] = name;
                continue;
            }
        }
        // Search through aliases for a match
        if (aliasData && aliasData[pid]) {
            let found = false;
            for (const alias of aliasData[pid]) {
                if (keggGeneSet.has(alias)) {
                    queryGeneNames.add(alias);
                    pidToKegg[pid] = alias;
                    found = true;
                    break;
                }
                // Also try uppercase version
                const upper = alias.toUpperCase();
                if (keggGeneSet.has(upper)) {
                    queryGeneNames.add(upper);
                    pidToKegg[pid] = upper;
                    found = true;
                    break;
                }
                // Also try lowercase version
                const lower = alias.toLowerCase();
                if (keggGeneSet.has(lower)) {
                    queryGeneNames.add(lower);
                    pidToKegg[pid] = lower;
                    found = true;
                    break;
                }
            }
            if (found) continue;
        }
    }

    const n = queryGeneNames.size;
    if (n === 0) {
        return { results: [], stats: { mapped: 0, total: queryProteinIds.length, termsTotal: Object.keys(pathwayBg).length } };
    }

    const results = [];
    for (const [pw, info] of Object.entries(pathwayBg)) {
        const K = info.genes.size;
        let k = 0;
        const geneHits = [];
        for (const gene of queryGeneNames) {
            if (info.genes.has(gene)) {
                k++;
                geneHits.push(gene);
            }
        }
        if (k === 0) continue;

        const pValue = hypergeomPValue(k, n, K, N);
        const expectedK = (K / N) * n;
        const fold = expectedK > 0 ? k / expectedK : 0;

        results.push({
            term: pw,
            description: info.name,
            category: 'KEGG Pathway',
            pValue,
            fdr: 1,
            fold: Math.round(fold * 100) / 100,
            geneCount: k,
            bgCount: K,
            totalGenes: n,
            totalBg: N,
            genes: geneHits,
        });
    }

    bhFDR(results);

    return {
        results,
        stats: { mapped: n, total: queryProteinIds.length, termsTotal: Object.keys(pathwayBg).length },
    };
}

// Export for use in app.js
window.Enrichment = { runGOEnrichment, runKEGGEnrichment };
