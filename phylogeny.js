/**
 * pomRelate — Phylogeny Analysis Module
 * Newick parser, tree layout, SVG renderer, and tab builder.
 * Renders precomputed gene trees per orthogroup with highlighted query/target genes.
 * No external dependencies.
 */

// ===== Newick Parser =====

/**
 * Parse a Newick-format string into a tree object.
 * Returns { name, branchLength, children[] } or null on failure.
 */
function parseNewick(str) {
    if (!str || typeof str !== 'string') return null;
    str = str.trim().replace(/;\s*$/, '');
    if (!str) return null;

    let pos = 0;

    function parseNode() {
        const node = { name: '', branchLength: null, children: [] };

        if (str[pos] === '(') {
            pos++; // skip '('
            node.children.push(parseNode());
            while (pos < str.length && str[pos] === ',') {
                pos++; // skip ','
                node.children.push(parseNode());
            }
            if (pos < str.length && str[pos] === ')') {
                pos++; // skip ')'
            }
        }

        // Read name (until ':', ',', ')', ';', or end)
        let name = '';
        while (pos < str.length && str[pos] !== ':' && str[pos] !== ',' && str[pos] !== ')' && str[pos] !== ';') {
            name += str[pos];
            pos++;
        }
        node.name = name.trim();

        // Read branch length
        if (pos < str.length && str[pos] === ':') {
            pos++; // skip ':'
            let lenStr = '';
            while (pos < str.length && str[pos] !== ',' && str[pos] !== ')' && str[pos] !== ';') {
                lenStr += str[pos];
                pos++;
            }
            const len = parseFloat(lenStr);
            if (!isNaN(len)) node.branchLength = len;
        }

        return node;
    }

    try {
        const tree = parseNode();
        return tree;
    } catch (e) {
        console.warn('Failed to parse Newick string:', e);
        return null;
    }
}

// ===== Tree Layout =====

/**
 * Count leaf nodes in a subtree.
 */
function countLeaves(node) {
    if (!node.children || node.children.length === 0) return 1;
    let count = 0;
    for (const child of node.children) count += countLeaves(child);
    return count;
}

/**
 * Find the maximum root-to-tip distance.
 */
function maxDepth(node, currentDepth) {
    const d = currentDepth + (node.branchLength || 0);
    if (!node.children || node.children.length === 0) return d;
    let max = 0;
    for (const child of node.children) {
        max = Math.max(max, maxDepth(child, d));
    }
    return max;
}

/**
 * Check if tree has meaningful branch lengths.
 */
function hasBranchLengths(node) {
    if (node.branchLength !== null && node.branchLength > 0) return true;
    if (node.children) {
        for (const child of node.children) {
            if (hasBranchLengths(child)) return true;
        }
    }
    return false;
}

/**
 * Assign x,y coordinates to each node for a rectangular tree layout.
 * Root at left, tips at right. Horizontal = depth, vertical = leaf index.
 *
 * @param {Object} root - parsed Newick tree
 * @param {number} plotW - available width for tree (excluding labels)
 * @param {number} plotH - available height
 * @param {boolean} useBranchLengths - scale x by branch length
 * @returns {Object[]} flat array of { node, x, y, parentX, parentY }
 */
function layoutTree(root, plotW, plotH, useBranchLengths) {
    const leafCount = countLeaves(root);
    const leafSpacing = plotH / Math.max(leafCount, 1);
    const totalDepth = useBranchLengths ? maxDepth(root, 0) : 0;
    const depthScale = useBranchLengths && totalDepth > 0 ? plotW / totalDepth : 0;

    const entries = [];
    let leafIndex = 0;

    // Compute max cladogram depth for equal spacing mode
    function maxCladoDepth(node, level) {
        if (!node.children || node.children.length === 0) return level;
        let max = 0;
        for (const child of node.children) {
            max = Math.max(max, maxCladoDepth(child, level + 1));
        }
        return max;
    }
    const maxLevel = useBranchLengths ? 0 : maxCladoDepth(root, 0);
    const levelStep = maxLevel > 0 ? plotW / maxLevel : plotW;

    function traverse(node, parentX, parentY, cumulativeLen, level) {
        let x, y;

        if (useBranchLengths) {
            x = (cumulativeLen + (node.branchLength || 0)) * depthScale;
        } else {
            x = level * levelStep;
        }

        if (!node.children || node.children.length === 0) {
            // Leaf node
            y = leafIndex * leafSpacing + leafSpacing / 2;
            leafIndex++;
            entries.push({ node, x, y, parentX, parentY: parentY });
        } else {
            // Internal node: children first, then compute y as midpoint
            const childEntries = [];
            for (const child of node.children) {
                const childEntry = traverse(
                    child, x, null,
                    useBranchLengths ? cumulativeLen + (node.branchLength || 0) : 0,
                    level + 1
                );
                childEntries.push(childEntry);
            }
            y = (childEntries[0].y + childEntries[childEntries.length - 1].y) / 2;

            // Update parentY for children
            for (const ce of childEntries) {
                const entry = entries.find(e => e.node === ce.node);
                if (entry) entry.parentY = y;
            }

            entries.push({ node, x, y, parentX, parentY: parentY });
        }

        return { node, x, y };
    }

    traverse(root, 0, 0, 0, 0);

    // Set root's parentX/parentY to its own position
    const rootEntry = entries.find(e => e.node === root);
    if (rootEntry) {
        rootEntry.parentX = rootEntry.x;
        rootEntry.parentY = rootEntry.y;
    }

    return entries;
}

// ===== SVG Helpers (local, matching plots.js pattern) =====

function _makeSVG(w, h) {
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    svg.setAttribute('width', w);
    svg.setAttribute('height', h);
    svg.setAttribute('viewBox', `0 0 ${w} ${h}`);
    return svg;
}

function _addGroup(parent, tx, ty) {
    const g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
    g.setAttribute('transform', `translate(${tx},${ty})`);
    parent.appendChild(g);
    return g;
}

function _addRect(parent, x, y, w, h, fill, className) {
    const r = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    r.setAttribute('x', x); r.setAttribute('y', y);
    r.setAttribute('width', w); r.setAttribute('height', h);
    r.setAttribute('fill', fill);
    if (className) r.classList.add(className);
    parent.appendChild(r);
    return r;
}

function _addLine(parent, x1, y1, x2, y2, stroke, width) {
    const l = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    l.setAttribute('x1', x1); l.setAttribute('y1', y1);
    l.setAttribute('x2', x2); l.setAttribute('y2', y2);
    l.setAttribute('stroke', stroke);
    l.setAttribute('stroke-width', width || 1);
    parent.appendChild(l);
    return l;
}

function _addCircle(parent, cx, cy, r, fill) {
    const c = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
    c.setAttribute('cx', cx); c.setAttribute('cy', cy);
    c.setAttribute('r', r); c.setAttribute('fill', fill);
    parent.appendChild(c);
    return c;
}

function _addText(parent, x, y, content, opts) {
    opts = opts || {};
    const t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    t.setAttribute('x', x); t.setAttribute('y', y);
    t.setAttribute('fill', opts.fill || '#000');
    t.setAttribute('font-size', opts.size || '10px');
    t.setAttribute('font-family', opts.family || "'EB Garamond', Georgia, serif");
    t.setAttribute('text-anchor', opts.anchor || 'start');
    if (opts.weight) t.setAttribute('font-weight', opts.weight);
    if (opts.baseline) t.setAttribute('dominant-baseline', opts.baseline);
    if (opts.style) t.setAttribute('font-style', opts.style);
    t.textContent = content;
    parent.appendChild(t);
    return t;
}

// ===== Phylogeny Tooltip =====

/**
 * Show a tooltip for a phylogeny tree leaf node.
 * Uses the existing gene-tooltip div from app.js.
 */
function _showPhyloTooltip(tipData, event) {
    var tooltipEl = document.getElementById('gene-tooltip');
    if (!tooltipEl) return;

    // Clear any pending hide from the main tooltip system
    if (window.cancelTooltipTimer) window.cancelTooltipTimer();
    clearTimeout(window._phyloTooltipTimer);

    window._phyloTooltipTimer = setTimeout(function() {
        var role = tipData.isQuery ? 'Query gene' : tipData.isTarget ? 'Target species' : 'Other';
        var roleColor = tipData.isQuery ? '#c92a2a' : tipData.isTarget ? '#1864ab' : '#666';
        var theme = document.documentElement.getAttribute('data-theme');
        if (theme === 'dark') {
            roleColor = tipData.isQuery ? '#ff8787' : tipData.isTarget ? '#4dabf7' : '#999';
        }

        var html = '<div class="gene-tooltip-header">';
        if (tipData.prefName) {
            html += '<span class="gene-tooltip-name">' + _esc(tipData.prefName) + '</span>';
        }
        html += '<span class="gene-tooltip-pid">' + _esc(tipData.shortId) + '</span>';
        html += '</div>';

        // Species
        html += '<div class="gene-tooltip-annotation" style="margin-bottom:4px;">';
        html += '<em>' + _esc(tipData.species) + '</em>';
        html += '</div>';

        // Role badge
        html += '<div style="margin-bottom:6px;">';
        html += '<span style="display:inline-block;padding:1px 8px;border-radius:3px;font-size:0.75rem;';
        html += 'background:' + roleColor + '22;color:' + roleColor + ';border:1px solid ' + roleColor + '44;">';
        html += role + '</span>';
        html += '</div>';

        // Details section
        html += '<div class="gene-tooltip-section">';
        html += '<div class="gene-tooltip-section-title">Details</div>';
        html += '<ul class="gene-tooltip-terms">';
        html += '<li><strong>Full ID:</strong> ' + _esc(tipData.geneName) + '</li>';
        if (tipData.ogId) {
            html += '<li><strong>Orthogroup:</strong> ' + _esc(tipData.ogId) + '</li>';
        }
        if (tipData.branchLength !== null && tipData.branchLength !== undefined) {
            html += '<li><strong>Branch length:</strong> ' + tipData.branchLength.toFixed(5) + '</li>';
        }
        html += '</ul></div>';

        // Reference links — use string_taxid for STRING URLs if available
        html += '<div class="gene-tooltip-refs">';
        if (tipData.speciesTaxid && tipData.shortId) {
            var stateObj = window.state || {};
            var spList = stateObj.speciesList || [];
            var spEntry = spList.find(function(s) { return s.taxid === tipData.speciesTaxid; });
            var strTaxid = (spEntry && spEntry.string_taxid) || tipData.speciesTaxid;
            html += '<a href="https://string-db.org/network/' + encodeURIComponent(strTaxid + '.' + tipData.shortId) + '" target="_blank" rel="noopener">STRING</a>';
        }
        if (tipData.shortId) {
            html += '<a href="https://www.uniprot.org/uniprotkb?query=' + encodeURIComponent(tipData.shortId) + '" target="_blank" rel="noopener">UniProt</a>';
        }
        html += '</div>';

        tooltipEl.innerHTML = html;
        tooltipEl.hidden = false;

        // Position near cursor
        var pad = 12;
        var x = event.clientX + pad;
        var y = event.clientY + pad;
        var vw = window.innerWidth;
        var vh = window.innerHeight;

        tooltipEl.style.left = x + 'px';
        tooltipEl.style.top = y + 'px';

        // Adjust after rendering to avoid overflow
        requestAnimationFrame(function() {
            var rect = tooltipEl.getBoundingClientRect();
            if (rect.right > vw - 8) {
                tooltipEl.style.left = Math.max(8, event.clientX - rect.width - pad) + 'px';
            }
            if (rect.bottom > vh - 8) {
                tooltipEl.style.top = Math.max(8, event.clientY - rect.height - pad) + 'px';
            }
            tooltipEl.classList.add('visible');
        });
    }, 200);
}

// ===== SVG Tree Renderer =====

/**
 * Render a phylogenetic tree as SVG.
 * @param {Object} root - parsed Newick tree
 * @param {Object} opts - { queryGenes: Set, isTargetSpecies: fn, getSpeciesName: fn, title: string, members: [], ogId: string }
 * @returns {SVGElement}
 */
function renderTreeSVG(root, opts) {
    opts = opts || {};
    const queryGenes = opts.queryGenes || new Set();
    const isTargetSpecies = opts.isTargetSpecies || (() => false);
    const getSpeciesName = opts.getSpeciesName || (() => '');
    const memberNames = opts.memberNames || {};
    const members = opts.members || [];
    const ogId = opts.ogId || '';

    const theme = document.documentElement.getAttribute('data-theme');
    const isDark = theme === 'dark';
    const textColor = isDark ? '#d4d4d4' : '#1a1a1a';
    const textMuted = isDark ? '#888888' : '#666666';
    const supportColor = isDark ? '#7c7c7c' : '#999999';
    const branchColor = isDark ? '#666666' : '#888888';
    const bgColor = isDark ? '#1a1a1a' : '#ffffff';
    const queryColor = isDark ? '#ff8787' : '#c92a2a';
    const targetColor = isDark ? '#4dabf7' : '#1864ab';
    const defaultNodeColor = isDark ? '#adb5bd' : '#868e96';

    const useBL = hasBranchLengths(root);
    const leafCount = countLeaves(root);

    // Adaptive spacing: denser for large trees
    const leafH = leafCount > 80 ? 16 : leafCount > 40 ? 18 : 22;
    const fontSize = leafCount > 80 ? '8px' : leafCount > 40 ? '9px' : '10px';

    const margin = { top: 40, right: 220, bottom: 50, left: 20 };
    const plotW = Math.max(200, Math.min(500, leafCount * 25));
    const plotH = leafCount * leafH;
    const width = margin.left + plotW + margin.right;
    const height = margin.top + plotH + margin.bottom;

    const entries = layoutTree(root, plotW, plotH, useBL);

    // Pre-scan support values to determine scale (0-1 vs 0-100)
    // If ANY internal node numeric value > 1, all are treated as 0-100 scale
    let supportScale01 = true;
    for (const entry of entries) {
        const nd = entry.node;
        if (nd.children && nd.children.length > 0 && nd.name) {
            const v = parseFloat(nd.name);
            if (!isNaN(v) && v > 1) { supportScale01 = false; break; }
        }
    }

    const svg = _makeSVG(width, height);
    _addRect(svg, 0, 0, width, height, bgColor, 'plot-bg');

    // Title
    if (opts.title) {
        _addText(svg, width / 2, 22, opts.title, {
            size: '13px', weight: '600', fill: textColor, anchor: 'middle'
        });
    }

    const g = _addGroup(svg, margin.left, margin.top);

    // Draw branches (rectangular: horizontal + vertical connector)
    for (const entry of entries) {
        if (entry.node === root) continue;
        const { x, y, parentX, parentY } = entry;
        if (parentX === null || parentX === undefined) continue;
        _addLine(g, parentX, y, x, y, branchColor, 1.2);
        _addLine(g, parentX, parentY, parentX, y, branchColor, 1.2);
    }

    // Draw nodes
    for (const entry of entries) {
        const { node, x, y } = entry;
        const isLeaf = !node.children || node.children.length === 0;
        const geneName = node.name || '';
        const speciesTaxid = geneName.split('.')[0];

        // Determine node color
        let nodeColor = defaultNodeColor;
        let nodeR = 3;
        if (queryGenes.has(geneName)) {
            nodeColor = queryColor;
            nodeR = 5;
        } else if (isTargetSpecies(speciesTaxid)) {
            nodeColor = targetColor;
            nodeR = 4;
        }

        if (isLeaf) {
            // Wrap leaf in interactive group for tooltip
            const leafG = _addGroup(g, 0, 0);
            leafG.style.cursor = 'pointer';

            _addCircle(leafG, x, y, nodeR, nodeColor);

            // Label: preferred name + species
            const spName = getSpeciesName(speciesTaxid);
            const shortId = geneName.includes('.') ? geneName.split('.').slice(1).join('.') : geneName;
            const prefName = memberNames[geneName];
            // Show preferred name before the ID when available
            const label = prefName
                ? (spName ? `${prefName} | ${shortId} (${spName})` : `${prefName} | ${shortId}`)
                : (spName ? `${shortId} (${spName})` : shortId);

            const textEl = _addText(leafG, x + 8, y, label, {
                size: fontSize, fill: textColor, baseline: 'middle'
            });

            if (queryGenes.has(geneName)) {
                textEl.setAttribute('font-weight', '700');
                textEl.setAttribute('fill', queryColor);
            } else if (isTargetSpecies(speciesTaxid)) {
                textEl.setAttribute('fill', targetColor);
            }

            // Invisible hit area for easier hovering
            const hitRect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
            hitRect.setAttribute('x', x - nodeR);
            hitRect.setAttribute('y', y - leafH / 2);
            hitRect.setAttribute('width', plotW + margin.right - x + nodeR);
            hitRect.setAttribute('height', leafH);
            hitRect.setAttribute('fill', 'transparent');
            leafG.insertBefore(hitRect, leafG.firstChild);

            const tipData = {
                geneName: geneName,
                shortId: shortId,
                prefName: prefName || '',
                species: spName || speciesTaxid,
                speciesTaxid: speciesTaxid,
                isQuery: queryGenes.has(geneName),
                isTarget: isTargetSpecies(speciesTaxid),
                branchLength: node.branchLength,
                ogId: ogId
            };

            leafG.addEventListener('mouseenter', function(e) {
                _showPhyloTooltip(tipData, e);
            });
            leafG.addEventListener('mouseleave', function(e) {
                var related = e.relatedTarget;
                var tooltip = document.getElementById('gene-tooltip');
                if (related && tooltip && tooltip.contains(related)) return;
                clearTimeout(window._phyloTooltipTimer);
                if (window.hideGeneTooltip) window.hideGeneTooltip();
            });
        } else {
            // Internal node: small dot
            _addCircle(g, x, y, 2, branchColor);

            // Show bootstrap/support values if node name is numeric
            if (geneName) {
                const support = parseFloat(geneName);
                if (!isNaN(support)) {
                    // Use pre-scanned scale to display as percentage
                    const displayVal = supportScale01
                        ? (support * 100).toFixed(0)
                        : Math.round(support).toString();
                    const threshold = supportScale01 ? 0.5 : 50;
                    if (support >= threshold) {
                        _addText(g, x - 4, y - 8, displayVal, {
                            size: '7px', fill: supportColor, anchor: 'end'
                        });
                    }
                }
            }
        }
    }

    // Scale bar (if branch lengths present)
    if (useBL) {
        const totalD = maxDepth(root, 0);
        if (totalD > 0) {
            // Compute a "nice" round number for the scale bar (~10-30% of tree depth)
            const target = totalD * 0.2;
            const mag = Math.pow(10, Math.floor(Math.log10(target)));
            const niceSteps = [1, 2, 5, 10];
            var barLen = niceSteps.reduce(function(best, s) {
                var v = s * mag;
                return Math.abs(v - target) < Math.abs(best - target) ? v : best;
            }, mag);
            const barPx = (barLen / totalD) * plotW;
            const barY = plotH + 20;
            const scaleColor = isDark ? '#b0b0b0' : '#333333';
            _addLine(g, 0, barY, barPx, barY, scaleColor, 2);
            _addLine(g, 0, barY - 4, 0, barY + 4, scaleColor, 1.5);
            _addLine(g, barPx, barY - 4, barPx, barY + 4, scaleColor, 1.5);
            // Format: remove trailing zeros
            var barLabel = barLen < 0.001 ? barLen.toExponential(1) : parseFloat(barLen.toPrecision(3)).toString();
            _addText(g, barPx / 2, barY + 16, barLabel, {
                size: '10px', fill: scaleColor, anchor: 'middle', weight: '500'
            });
        }
    }

    return svg;
}

// ===== Query Gene Matching =====

/**
 * Find tree tip names that correspond to the query gene.
 * eggNOG v7 trees use taxid.UniProtID format — same as STRING protein IDs,
 * so direct matching (S1) should work in most cases.
 * Fallback strategies are kept for edge cases.
 * @returns {Object} { tips: Set, strategy: 'exact'|'id'|'species'|'fallback' }
 */
function _findQueryTipNames(matchedId, members, newickStr, sourceTaxid) {
    // Extract tip names from Newick
    var tipMatches = newickStr.match(/[\(,]([A-Za-z0-9_.]+):/g);
    if (!tipMatches) return { tips: new Set([matchedId]), strategy: 'fallback' };
    var tipNames = tipMatches.map(function(t) { return t.replace(/^[\(,]/, '').replace(/:$/, ''); });

    // S1: direct match (primary — works for eggNOG v7 trees)
    if (tipNames.indexOf(matchedId) >= 0) return { tips: new Set([matchedId]), strategy: 'exact' };

    // S2: protein ID portion match (fallback)
    var proteinId = matchedId.indexOf('.') >= 0 ? matchedId.split('.').pop() : matchedId;
    for (var ti = 0; ti < tipNames.length; ti++) {
        var tipId = tipNames[ti].indexOf('.') >= 0 ? tipNames[ti].split('.').pop() : tipNames[ti];
        if (tipId === proteinId) return { tips: new Set([tipNames[ti]]), strategy: 'id' };
    }

    // S3: match all tips from the source species (same taxid prefix)
    var taxPrefix = sourceTaxid + '.';
    var sourceTips = tipNames.filter(function(t) { return t.indexOf(taxPrefix) === 0; });
    if (sourceTips.length > 0) return { tips: new Set(sourceTips), strategy: 'species' };

    return { tips: new Set([matchedId]), strategy: 'fallback' };
}

// ===== Taxid Name Resolution =====

/**
 * Static NCBI taxid → scientific name mapping.
 */
var _taxidNames = null;

// ===== Tab Builder =====

/**
 * Look up orthogroup for a resolved gene.
 * Tries multiple ID forms: taxid.proteinId, proteinId, preferred name.
 */
function findOrthogroup(proteinId, sourceTaxid, phyloData, aliasData) {
    const geneToOg = phyloData.orthogroups.gene_to_og;
    if (!geneToOg) return null;

    // Try taxid.proteinId
    const fullId = `${sourceTaxid}.${proteinId}`;
    if (geneToOg[fullId]) return { ogId: geneToOg[fullId], matchedId: fullId };

    // Try raw proteinId
    if (geneToOg[proteinId]) return { ogId: geneToOg[proteinId], matchedId: proteinId };

    // Try all keys that end with the proteinId
    for (const key of Object.keys(geneToOg)) {
        if (key.endsWith('.' + proteinId)) return { ogId: geneToOg[key], matchedId: key };
    }

    return null;
}

/**
 * Build the Phylogeny tab content.
 * Called from app.js after analysis completes.
 */
function buildPhylogenyTab(resolvedGenes, sourceTaxid, targetTaxids, phyloData) {
    const container = document.querySelector('#tab-phylogeny');
    if (!container) return;

    if (!phyloData || !phyloData.orthogroups || !phyloData.trees) {
        container.innerHTML = '<p class="no-data">No phylogeny data available. Precomputed orthogroup and gene tree files are required.</p>';
        return;
    }

    const found = resolvedGenes.filter(function(g) { return g.proteinId; });
    if (found.length === 0) {
        container.innerHTML = '<p class="no-data">No resolved genes to display phylogeny for.</p>';
        return;
    }

    // Species name lookup
    // eggNOG v7 trees use strain-level taxids matching STRING, so direct lookup works
    const taxidNames = (phyloData && phyloData.taxidNames) || {};
    const getSpeciesName = function(taxid) {
        const sp = (window.state || {}).speciesList;
        if (sp) {
            var match = sp.find(function(s) { return s.taxid === taxid; });
            if (match) return match.compact_name;
        }
        // Fallback to static NCBI taxid→name mapping
        if (taxidNames[taxid]) return taxidNames[taxid];
        return taxid;
    };

    const targetSet = new Set(targetTaxids || []);
    const isTargetSpecies = function(taxid) {
        return targetSet.has(taxid);
    };

    let html = '';
    const treeElements = []; // { containerId, svg, newick, ogId }

    for (const gene of found) {
        const proteinId = gene.proteinId;
        const query = gene.query;

        const result = findOrthogroup(proteinId, sourceTaxid, phyloData);
        if (!result) {
            html += '<div class="result-section">';
            html += '<div class="result-section-title"><span class="result-gene-badge">' + _esc(query) + '</span></div>';
            html += '<p class="gene-not-found">No orthogroup found for this gene.</p>';
            html += '</div>';
            continue;
        }

        const ogId = result.ogId;
        const matchedId = result.matchedId;
        const members = phyloData.orthogroups.og_members[ogId] || [];
        const newick = phyloData.trees[ogId];

        html += '<div class="result-section">';
        html += '<div class="result-section-title">';
        html += '<span class="result-gene-badge">' + _esc(query) + '</span>';
        html += ' <span class="tag tag-phylo">' + _esc(ogId) + '</span>';
        html += '</div>';
        html += '<div style="font-size:0.78rem;color:var(--text-muted);margin:-0.3rem 0 0.4rem 0;">Pre-computed gene tree from eggNOG v7 (Hern\u00e1ndez-Plaza et al., 2026)</div>';

        // Export bar
        const treeContainerId = 'phylo-tree-' + ogId + '-' + proteinId.replace(/[^a-zA-Z0-9]/g, '_');
        html += '<div class="download-bar">';
        html += '<span class="download-label">Export Tree:</span>';
        html += '<button class="download-btn" onclick="window.Phylogeny.exportTree(\'' + treeContainerId + '\', \'png\')">PNG</button>';
        html += '<button class="download-btn" onclick="window.Phylogeny.exportTree(\'' + treeContainerId + '\', \'svg\')">SVG</button>';
        html += '<button class="download-btn" onclick="window.Phylogeny.exportTree(\'' + treeContainerId + '\', \'pdf\')">PDF</button>';
        if (newick) {
            html += '<button class="download-btn" onclick="window.Phylogeny.exportNewick(\'' + _esc(ogId) + '\')">Newick</button>';
            html += '<button class="download-btn" onclick="window.Phylogeny.exportNexus(\'' + _esc(ogId) + '\')">NEXUS</button>';
            html += '<button class="download-btn" onclick="window.Phylogeny.exportPhyloXML(\'' + _esc(ogId) + '\')">PhyloXML</button>';
        }
        html += '</div>';

        // Tree container
        html += '<div id="' + treeContainerId + '" class="phylo-tree-container plot-container"></div>';

        // Summary table
        const representedSpecies = new Set(members.map(function(m) { return m.species; }));
        const missingTargets = targetTaxids.filter(function(t) { return !representedSpecies.has(t); });

        html += '<div class="table-responsive"><table class="result-table"><thead><tr>';
        html += '<th>Orthogroup</th><th>Members</th><th>Species Represented</th><th>Missing Targets</th>';
        html += '</tr></thead><tbody>';
        html += '<tr>';
        html += '<td><code>' + _esc(ogId) + '</code></td>';
        html += '<td>' + members.length + '</td>';
        html += '<td>' + representedSpecies.size + ' / ' + (targetTaxids.length + 1) + '</td>';
        html += '<td>';
        if (missingTargets.length === 0) {
            html += '<span style="color: var(--text-muted); font-style: italic;">None</span>';
        } else {
            html += missingTargets.map(function(t) { return _escItalicSpecies(getSpeciesName(t)); }).join(', ');
        }
        html += '</td>';
        html += '</tr></tbody></table></div>';

        // Member gene list (compact)
        if (members.length > 0) {
            html += '<details style="margin-bottom: 0.5rem;"><summary style="cursor: pointer; font-size: 0.85rem; color: var(--text-secondary);">Show ' + members.length + ' member gene(s)</summary>';
            html += '<div class="table-responsive"><table class="result-table"><thead><tr>';
            html += '<th>Gene</th><th>Species</th><th>Name</th>';
            html += '</tr></thead><tbody>';
            for (const m of members) {
                const isQuery = m.gene === matchedId;
                const rowStyle = isQuery ? ' style="background: rgba(255, 153, 153, 0.1)"' : '';
                html += '<tr' + rowStyle + '>';
                html += '<td><code>' + _esc(m.gene) + '</code></td>';
                html += '<td>' + _escItalicSpecies(getSpeciesName(m.species)) + '</td>';
                html += '<td>' + _esc(m.name || '') + '</td>';
                html += '</tr>';
            }
            html += '</tbody></table></div></details>';
        }

        html += '</div>';

        // Build member name mapping for better tip labels
        var memberNameMap = {};
        for (var mi = 0; mi < members.length; mi++) {
            if (members[mi].name) memberNameMap[members[mi].gene] = members[mi].name;
        }

        // Queue tree rendering
        if (newick) {
            var queryMatch = _findQueryTipNames(matchedId, members, newick, sourceTaxid);
            treeElements.push({
                containerId: treeContainerId,
                newick: newick,
                ogId: ogId,
                matchedId: matchedId,
                queryGenes: queryMatch.tips,
                matchStrategy: queryMatch.strategy,
                isTargetSpecies: isTargetSpecies,
                getSpeciesName: getSpeciesName,
                memberNames: memberNameMap,
                members: members,
                title: ogId + ' Gene Tree'
            });
        }
    }

    const missing = resolvedGenes.filter(function(g) { return !g.proteinId; }).map(function(g) { return g.query; });
    if (missing.length > 0) {
        html += '<div class="not-found-summary">Not found in source species: ' + _esc(missing.join(', ')) + '</div>';
    }

    container.innerHTML = html || '<p class="no-data">No phylogeny results.</p>';

    for (const te of treeElements) {
        const treeContainer = document.getElementById(te.containerId);
        if (!treeContainer) continue;

        const tree = parseNewick(te.newick);
        if (!tree) {
            treeContainer.innerHTML = '<p class="no-data">Failed to parse gene tree.</p>';
            continue;
        }

        const tipCount = countLeaves(tree);

        if (tipCount > 200) {
            treeContainer.innerHTML = '<div class="phylo-large-tree-notice">'
                + '<p><strong>' + tipCount + ' tips</strong> — this tree is too large for in-browser rendering.</p>'
                + '<p>Download the tree in Newick, NEXUS, or PhyloXML format and open it in a dedicated viewer '
                + 'such as <a href="https://itol.embl.de" target="_blank" rel="noopener">iTOL</a>, '
                + '<a href="http://tree.bio.ed.ac.uk/software/figtree/" target="_blank" rel="noopener">FigTree</a>, '
                + 'or <a href="https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/dendroscope/" target="_blank" rel="noopener">Dendroscope</a>.</p>'
                + '</div>';
            continue;
        }

        const svg = renderTreeSVG(tree, {
            queryGenes: te.queryGenes,
            isTargetSpecies: te.isTargetSpecies,
            getSpeciesName: te.getSpeciesName,
            memberNames: te.memberNames,
            members: te.members,
            ogId: te.ogId,
            title: te.title
        });

        if (svg) {
            // Show warning if species-level fallback was used
            if (te.matchStrategy === 'species') {
                var warn = document.createElement('div');
                warn.style.cssText = 'font-size:0.78rem;color:var(--text-muted);padding:0.25rem 0.5rem;font-style:italic;';
                warn.textContent = 'Exact gene match not found in tree — highlighting all co-orthologs from the query species.';
                treeContainer.appendChild(warn);
            }
            // Show notice when tree has fewer species than the orthogroup
            var treeSpecies = new Set();
            var tipRe = /[\(,](\d+)\./g;
            var tipMatch;
            while ((tipMatch = tipRe.exec(te.newick)) !== null) {
                treeSpecies.add(tipMatch[1]);
            }
            var ogSpecies = new Set(te.members.map(function(m) { return m.species; }));
            if (treeSpecies.size < ogSpecies.size) {
                var missingFromTree = [];
                ogSpecies.forEach(function(sp) {
                    if (!treeSpecies.has(sp)) missingFromTree.push(te.getSpeciesName(sp));
                });
                var notice = document.createElement('div');
                notice.style.cssText = 'font-size:0.78rem;color:var(--text-muted);padding:0.3rem 0.5rem;margin-bottom:0.3rem;border-left:3px solid var(--border);background:var(--bg-panel);border-radius:2px;';
                notice.innerHTML = '<strong>Tree covers ' + treeSpecies.size + '/' + ogSpecies.size + ' species.</strong> '
                    + missingFromTree.length + ' species have orthologs in this group but were not included in the '
                    + 'eggNOG v7 phylogenetic analysis for this protein family: '
                    + '<em>' + missingFromTree.join(', ') + '</em>.';
                treeContainer.appendChild(notice);
            }
            treeContainer.appendChild(svg);
        } else {
            treeContainer.innerHTML = '<p class="no-data">Could not render tree.</p>';
        }
    }

    container.querySelectorAll('table').forEach(function(table) {
        if (typeof makeTableSortable === 'function') makeTableSortable(table);
    });
}

// ===== Export Functions =====

function exportTree(containerId, format) {
    const container = document.getElementById(containerId);
    if (!container) return alert('Tree not found.');
    const svg = container.querySelector('svg');
    if (!svg) return alert('No tree to export.');

    const name = containerId.replace('phylo-tree-', 'phylogeny_');
    if (format === 'png') window.Export.downloadPNG(svg, name + '.png');
    else if (format === 'svg') window.Export.downloadSVG(svg, name + '.svg');
    else if (format === 'pdf') window.Export.downloadPDF(svg, name + '.pdf');
}

function exportNewick(ogId) {
    var phyloData = (window.state || {}).phylogenyData;
    if (!phyloData || !phyloData.trees || !phyloData.trees[ogId]) {
        return alert('Newick data not found for ' + ogId);
    }

    var newickStr = phyloData.trees[ogId];
    if (!newickStr.trim().endsWith(';')) newickStr += ';';

    var blob = new Blob([newickStr], { type: 'text/plain;charset=utf-8;' });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = ogId + '_tree.nwk';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(function() { URL.revokeObjectURL(a.href); }, 1000);
}

function exportPhylogenyCSV(resolvedGenes, sourceTaxid, targetTaxids, phyloData) {
    if (!phyloData || !phyloData.orthogroups) return alert('No phylogeny data to export.');

    var taxidNames = (phyloData && phyloData.taxidNames) || {};
    var getSpeciesName = function(taxid) {
        var sp = (window.state || {}).speciesList;
        if (sp) {
            var match = sp.find(function(s) { return s.taxid === taxid; });
            if (match) return match.compact_name;
        }
        if (taxidNames[taxid]) return taxidNames[taxid];
        return taxid;
    };

    var headers = ['Gene', 'Protein ID', 'Orthogroup', 'Members', 'Species Represented', 'Missing Targets'];
    var rows = [];

    var found = resolvedGenes.filter(function(g) { return g.proteinId; });
    for (var i = 0; i < found.length; i++) {
        var gene = found[i];
        var result = findOrthogroup(gene.proteinId, sourceTaxid, phyloData);
        var ogId = result ? result.ogId : '';
        var members = result ? (phyloData.orthogroups.og_members[ogId] || []) : [];
        var representedSpecies = new Set(members.map(function(m) { return m.species; }));
        var missingTargets = targetTaxids.filter(function(t) { return !representedSpecies.has(t); });

        rows.push([
            gene.query,
            gene.proteinId,
            ogId,
            members.length,
            representedSpecies.size,
            '"' + missingTargets.map(function(t) { return getSpeciesName(t); }).join(', ') + '"'
        ]);
    }

    var csv = [headers.join(',')].concat(rows.map(function(r) { return r.join(','); })).join('\n');
    var blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'phylogeny_summary.csv';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(function() { URL.revokeObjectURL(a.href); }, 1000);
}

// ===== Helpers =====

function _esc(str) {
    if (!str) return '';
    var div = document.createElement('div');
    div.textContent = str;
    return div.innerHTML;
}

function _escItalicSpecies(name) {
    if (!name) return '';
    var parts = name.split(' ');
    if (parts.length >= 2) {
        return '<em>' + _esc(parts[0]) + '</em> ' + _esc(parts.slice(1).join(' '));
    }
    return '<em>' + _esc(name) + '</em>';
}

// ===== Additional Export Formats =====

/**
 * Export tree in NEXUS format (compatible with PAUP*, MrBayes, FigTree, etc.)
 */
function exportNexus(ogId) {
    var phyloData = (window.state || {}).phylogenyData;
    if (!phyloData || !phyloData.trees || !phyloData.trees[ogId]) {
        return alert('Tree data not found for ' + ogId);
    }

    var newickStr = phyloData.trees[ogId].trim();
    if (!newickStr.endsWith(';')) newickStr += ';';

    var tree = parseNewick(newickStr);
    var tips = [];
    function collectTips(node) {
        if (!node) return;
        if (!node.children || node.children.length === 0) {
            if (node.name) tips.push(node.name);
        } else {
            for (var i = 0; i < node.children.length; i++) collectTips(node.children[i]);
        }
    }
    collectTips(tree);

    var nexus = '#NEXUS\n\n';
    nexus += 'BEGIN TAXA;\n';
    nexus += '  DIMENSIONS NTAX=' + tips.length + ';\n';
    nexus += '  TAXLABELS\n';
    for (var i = 0; i < tips.length; i++) {
        nexus += "    '" + tips[i].replace(/'/g, "''") + "'\n";
    }
    nexus += '  ;\n';
    nexus += 'END;\n\n';
    nexus += 'BEGIN TREES;\n';
    nexus += '  TREE ' + ogId + ' = [&U] ' + newickStr + '\n';
    nexus += 'END;\n';

    _downloadText(nexus, ogId + '_tree.nex', 'text/plain;charset=utf-8;');
}

/**
 * Export tree in PhyloXML format (compatible with Archaeopteryx, ETE, forester, etc.)
 */
function exportPhyloXML(ogId) {
    var phyloData = (window.state || {}).phylogenyData;
    if (!phyloData || !phyloData.trees || !phyloData.trees[ogId]) {
        return alert('Tree data not found for ' + ogId);
    }

    var tree = parseNewick(phyloData.trees[ogId]);
    if (!tree) return alert('Failed to parse tree for ' + ogId);

    function nodeToXML(node, indent) {
        var pad = '  '.repeat(indent);
        var xml = pad + '<clade>\n';
        if (node.name) {
            xml += pad + '  <name>' + _escXML(node.name) + '</name>\n';
        }
        if (node.branchLength !== null && node.branchLength !== undefined) {
            xml += pad + '  <branch_length>' + node.branchLength + '</branch_length>\n';
        }
        if (node.children) {
            for (var i = 0; i < node.children.length; i++) {
                xml += nodeToXML(node.children[i], indent + 1);
            }
        }
        xml += pad + '</clade>\n';
        return xml;
    }

    var xml = '<?xml version="1.0" encoding="UTF-8"?>\n';
    xml += '<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ';
    xml += 'xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.20/phyloxml.xsd" ';
    xml += 'xmlns="http://www.phyloxml.org">\n';
    xml += '  <phylogeny rooted="false">\n';
    xml += '    <name>' + _escXML(ogId) + '</name>\n';
    xml += '    <description>Gene tree from eggNOG/STRING orthology</description>\n';
    xml += nodeToXML(tree, 2);
    xml += '  </phylogeny>\n';
    xml += '</phyloxml>\n';

    _downloadText(xml, ogId + '_tree.xml', 'application/xml;charset=utf-8;');
}

function _escXML(str) {
    if (!str) return '';
    return str.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}

function _downloadText(content, filename, mimeType) {
    var blob = new Blob([content], { type: mimeType });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(function() { URL.revokeObjectURL(a.href); }, 1000);
}

// ===== Export Module =====

window.Phylogeny = {
    parseNewick: parseNewick,
    buildPhylogenyTab: buildPhylogenyTab,
    findOrthogroup: findOrthogroup,
    exportTree: exportTree,
    exportNewick: exportNewick,
    exportNexus: exportNexus,
    exportPhyloXML: exportPhyloXML,
    exportPhylogenyCSV: exportPhylogenyCSV,
};
