/**
 * pomRelate — Publication-Quality Plot Generation
 * SVG-based enrichment plots: bar chart, dot plot, and hierarchical clustering dendrogram.
 * Designed for academic publication standards.
 */

const PALETTES = {
    'Default': (t, theme) => {
        // Sequential palette: light-to-dark grays with slight blue tint for publication
        t = Math.max(0, Math.min(1, t));
        if (theme === 'dark') {
            const r = Math.round(60 + t * 150);
            const g = Math.round(60 + t * 155);
            const b = Math.round(70 + t * 160);
            return `rgb(${r},${g},${b})`;
        } else {
            const r = Math.round(210 - t * 170);
            const g = Math.round(215 - t * 175);
            const b = Math.round(220 - t * 170);
            return `rgb(${r},${g},${b})`;
        }
    },
    'Viridis': (t) => interpolate(t, [[68, 1, 84], [72, 35, 116], [64, 67, 135], [53, 94, 141], [42, 118, 142], [33, 144, 141], [31, 161, 135], [53, 183, 121], [109, 205, 89], [180, 222, 44], [253, 231, 37]]),
    'Magma': (t) => interpolate(t, [[0, 0, 4], [27, 11, 63], [80, 18, 107], [141, 35, 105], [201, 68, 88], [243, 115, 70], [254, 173, 84], [251, 252, 191]]),
    'Plasma': (t) => interpolate(t, [[13, 8, 135], [83, 2, 163], [156, 23, 158], [204, 71, 120], [237, 121, 83], [253, 180, 47], [240, 249, 33]]),
    'Blues': (t) => interpolate(t, [[247, 251, 255], [222, 235, 247], [198, 219, 239], [158, 202, 225], [107, 174, 214], [66, 146, 198], [33, 113, 181], [8, 81, 156], [8, 48, 107]]),
    'Reds': (t) => interpolate(t, [[255, 245, 240], [254, 224, 210], [252, 187, 161], [252, 146, 114], [251, 106, 74], [239, 59, 44], [203, 24, 29], [165, 15, 21], [103, 0, 13]]),
    'Greys': (t) => interpolate(t, [[255, 255, 255], [240, 240, 240], [217, 217, 217], [189, 189, 189], [150, 150, 150], [115, 115, 115], [82, 82, 82], [37, 37, 37], [0, 0, 0]])
};

function interpolate(t, colors) {
    t = Math.max(0, Math.min(1, t));
    const i = t * (colors.length - 1);
    const i0 = Math.floor(i);
    const i1 = Math.min(i0 + 1, colors.length - 1);
    const f = i - i0;
    const c0 = colors[i0];
    const c1 = colors[i1];
    const r = Math.round(c0[0] + f * (c1[0] - c0[0]));
    const g = Math.round(c0[1] + f * (c1[1] - c0[1]));
    const b = Math.round(c0[2] + f * (c1[2] - c0[2]));
    return `rgb(${r},${g},${b})`;
}

/**
 * Create a horizontal bar chart of top enriched terms.
 * Publication-quality: proper axes, tick marks, legends, no overlapping.
 */
function createBarChart(results, topN = 20, palette = 'Default', title = 'Enrichment Analysis') {
    const data = results.filter(r => r.fdr <= 1).slice(0, topN).reverse(); // reverse for bottom-to-top
    if (data.length === 0) return null;

    const cs = getComputedStyle(document.documentElement);
    const theme = document.documentElement.getAttribute('data-theme');
    const textColor = theme === 'dark' ? '#d4d4d4' : '#1a1a1a';
    const textMuted = theme === 'dark' ? '#888888' : '#666666';
    const axisColor = theme === 'dark' ? '#555555' : '#333333';
    const gridColor = theme === 'dark' ? '#333333' : '#e0e0e0';
    const bgColor = theme === 'dark' ? '#1a1a1a' : '#ffffff';

    const margin = { top: 50, right: 110, bottom: 58, left: 280 };
    const barH = 20;
    const barGap = 7;
    const plotH = data.length * (barH + barGap) - barGap;
    const plotW = 440;
    const width = margin.left + plotW + margin.right;
    const height = margin.top + plotH + margin.bottom;

    const maxVal = Math.max(...data.map(d => -Math.log10(Math.max(d.fdr, 1e-300))));
    const xMax = niceMax(maxVal);
    const xScale = (v) => (v / xMax) * plotW;

    const svg = makeSVG(width, height);

    // Background
    addRect(svg, 0, 0, width, height, bgColor, 'plot-bg');

    // Title
    addText(svg, width / 2, 22, title, {
        size: '14px', weight: '700', fill: textColor, anchor: 'middle',
        family: "'EB Garamond', Georgia, serif"
    });

    // Subtitle: FDR threshold
    const sigCount = data.filter(d => d.fdr < 0.05).length;
    addText(svg, width / 2, 38, `Top ${data.length} terms · ${sigCount} significant (FDR < 0.05)`, {
        size: '10px', fill: textMuted, anchor: 'middle',
        family: "'EB Garamond', Georgia, serif"
    });

    const g = addGroup(svg, margin.left, margin.top);

    // X-axis gridlines and ticks
    const ticks = niceTicksFor(0, xMax, 5);
    for (const t of ticks) {
        const x = xScale(t);
        addLine(g, x, -3, x, plotH, gridColor, 0.5, '2,3');
        addLine(g, x, plotH, x, plotH + 5, axisColor, 1);
        addText(g, x, plotH + 18, formatTick(t), {
            size: '10px', fill: textColor, anchor: 'middle',
            family: "'EB Garamond', Georgia, serif"
        });
    }

    // X-axis label
    addText(svg, margin.left + plotW / 2, height - 10, '−log\u2081\u2080(FDR)', {
        size: '11px', fill: textColor, anchor: 'middle', weight: '500',
        family: "'EB Garamond', Georgia, serif"
    });

    // Color scale: fold enrichment → shades
    const maxFold = Math.max(...data.map(d => d.fold));
    const minFold = Math.min(...data.map(d => d.fold));

    // Bars
    for (let i = 0; i < data.length; i++) {
        const d = data[i];
        const y = i * (barH + barGap);
        const val = -Math.log10(Math.max(d.fdr, 1e-300));
        const w = Math.max(xScale(val), 2);

        const t = maxFold > minFold ? (d.fold - minFold) / (maxFold - minFold) : 0.5;
        const colorFn = PALETTES[palette] || PALETTES['Default'];
        const barColor = colorFn(t, theme);

        const r = addRect(g, 0, y, w, barH, barColor);
        r.setAttribute('rx', '2');

        // Significance marker
        if (d.fdr < 0.001) {
            addText(g, w + 4, y + barH / 2 + 1, '***', { size: '9px', fill: textMuted, anchor: 'start', baseline: 'middle' });
        } else if (d.fdr < 0.01) {
            addText(g, w + 4, y + barH / 2 + 1, '**', { size: '9px', fill: textMuted, anchor: 'start', baseline: 'middle' });
        } else if (d.fdr < 0.05) {
            addText(g, w + 4, y + barH / 2 + 1, '*', { size: '9px', fill: textMuted, anchor: 'start', baseline: 'middle' });
        }

        // Gene count right of bar
        addText(g, w + (d.fdr < 0.05 ? 22 : 4), y + barH / 2 + 1, `n=${d.geneCount}`, {
            size: '9px', fill: textMuted, anchor: 'start', baseline: 'middle',
            family: "'EB Garamond', Georgia, serif"
        });

        // Term label (left of bars)
        const label = truncLabel(d.description, 42);
        addText(g, -8, y + barH / 2 + 1, label, {
            size: '10px', fill: textColor, anchor: 'end', baseline: 'middle',
            family: "'EB Garamond', Georgia, serif"
        });
    }

    // Axes
    addLine(g, 0, -3, 0, plotH, axisColor, 1.2);
    addLine(g, 0, plotH, plotW, plotH, axisColor, 1.2);

    // Legend: fold enrichment gradient
    const legendX = plotW + 45;
    const legendY = 0;
    const legendH = Math.min(plotH * 0.5, 100);
    const legendW = 12;

    addText(g, legendX + legendW / 2, legendY - 6, 'Fold Enrichment', {
        size: '9px', fill: textColor, anchor: 'middle', weight: '600',
        family: "'EB Garamond', Georgia, serif"
    });

    const gradSteps = 20;
    for (let i = 0; i < gradSteps; i++) {
        const frac = i / (gradSteps - 1);
        const cy = legendY + frac * legendH;
        const colorFn = PALETTES[palette] || PALETTES['Default'];
        const color = colorFn(1 - frac, theme);
        addRect(g, legendX, cy, legendW, legendH / gradSteps + 1, color);
    }

    // Legend tick labels
    addText(g, legendX + legendW + 4, legendY + 4, maxFold.toFixed(1), {
        size: '8px', fill: textMuted, anchor: 'start', baseline: 'middle',
        family: "'EB Garamond', Georgia, serif"
    });
    addText(g, legendX + legendW + 4, legendY + legendH, minFold.toFixed(1), {
        size: '8px', fill: textMuted, anchor: 'start', baseline: 'middle',
        family: "'EB Garamond', Georgia, serif"
    });

    // Legend border
    addRect(g, legendX, legendY, legendW, legendH, 'none').setAttribute('stroke', axisColor);

    return svg;
}

/**
 * Create a dot plot of enrichment results.
 * Publication-quality with proper size and color legends.
 */
function createDotPlot(results, topN = 20, palette = 'Default', title = 'Enrichment Dot Plot') {
    const data = results.filter(r => r.fdr <= 1).slice(0, topN).reverse();
    if (data.length === 0) return null;

    const theme = document.documentElement.getAttribute('data-theme');
    const textColor = theme === 'dark' ? '#d4d4d4' : '#1a1a1a';
    const textMuted = theme === 'dark' ? '#888888' : '#666666';
    const axisColor = theme === 'dark' ? '#555555' : '#333333';
    const gridColor = theme === 'dark' ? '#333333' : '#e0e0e0';
    const bgColor = theme === 'dark' ? '#1a1a1a' : '#ffffff';

    const margin = { top: 50, right: 130, bottom: 58, left: 280 };
    const rowH = 26;
    const plotH = data.length * rowH;
    const plotW = 380;
    const width = margin.left + plotW + margin.right;
    const height = margin.top + plotH + margin.bottom;

    const maxFold = Math.max(...data.map(d => d.fold));
    const xMax = niceMax(maxFold);
    const maxGeneCount = Math.max(...data.map(d => d.geneCount));
    const maxLogFDR = Math.max(...data.map(d => -Math.log10(Math.max(d.fdr, 1e-300))));

    const xScale = (v) => (v / xMax) * plotW;
    const rScale = (v) => 3 + (v / Math.max(maxGeneCount, 1)) * 9;

    const svg = makeSVG(width, height);

    addRect(svg, 0, 0, width, height, bgColor, 'plot-bg');

    addText(svg, width / 2, 22, title, {
        size: '14px', weight: '700', fill: textColor, anchor: 'middle',
        family: "'EB Garamond', Georgia, serif"
    });

    const sigCount = data.filter(d => d.fdr < 0.05).length;
    addText(svg, width / 2, 38, `Top ${data.length} terms · ${sigCount} significant (FDR < 0.05)`, {
        size: '10px', fill: textMuted, anchor: 'middle',
        family: "'EB Garamond', Georgia, serif"
    });

    const g = addGroup(svg, margin.left, margin.top);

    // Gridlines
    const ticks = niceTicksFor(0, xMax, 5);
    for (const t of ticks) {
        const x = xScale(t);
        addLine(g, x, -3, x, plotH, gridColor, 0.5, '2,3');
        addLine(g, x, plotH, x, plotH + 5, axisColor, 1);
        addText(g, x, plotH + 18, formatTick(t), {
            size: '10px', fill: textColor, anchor: 'middle',
            family: "'EB Garamond', Georgia, serif"
        });
    }

    addText(svg, margin.left + plotW / 2, height - 10, 'Fold Enrichment', {
        size: '11px', fill: textColor, anchor: 'middle', weight: '500',
        family: "'EB Garamond', Georgia, serif"
    });

    // Horizontal gridlines for each term
    for (let i = 0; i < data.length; i++) {
        const y = i * rowH + rowH / 2;
        addLine(g, 0, y, plotW, y, gridColor, 0.3);
    }

    // Dots
    for (let i = 0; i < data.length; i++) {
        const d = data[i];
        const y = i * rowH + rowH / 2;
        const x = xScale(d.fold);
        const r = rScale(d.geneCount);

        const logFDR = -Math.log10(Math.max(d.fdr, 1e-300));
        const intensity = maxLogFDR > 0 ? Math.min(logFDR / maxLogFDR, 1) : 0.5;
        const colorFn = PALETTES[palette] || PALETTES['Default'];
        const dotColor = colorFn(intensity, theme);

        const c = addCircle(g, x, y, r, dotColor);
        c.setAttribute('stroke', axisColor);
        c.setAttribute('stroke-width', '0.5');

        // Term label
        const label = truncLabel(d.description, 42);
        addText(g, -8, y + 1, label, {
            size: '10px', fill: textColor, anchor: 'end', baseline: 'middle',
            family: "'EB Garamond', Georgia, serif"
        });
    }

    // Axes
    addLine(g, 0, -3, 0, plotH, axisColor, 1.2);
    addLine(g, 0, plotH, plotW, plotH, axisColor, 1.2);

    // ===== Legends =====
    const legX = plotW + 18;

    // Size legend
    addText(g, legX, 0, 'Gene Count', {
        size: '9px', fill: textColor, anchor: 'start', weight: '600',
        family: "'EB Garamond', Georgia, serif"
    });

    const sizeLevels = maxGeneCount <= 3
        ? Array.from({ length: maxGeneCount }, (_, i) => i + 1)
        : [1, Math.ceil(maxGeneCount / 2), maxGeneCount];

    for (let i = 0; i < sizeLevels.length; i++) {
        const ly = 18 + i * 24;
        const lr = rScale(sizeLevels[i]);
        addCircle(g, legX + 10, ly, lr, 'none').setAttribute('stroke', axisColor);
        addText(g, legX + 25, ly + 1, String(sizeLevels[i]), {
            size: '9px', fill: textMuted, anchor: 'start', baseline: 'middle',
            family: "'EB Garamond', Georgia, serif"
        });
    }

    // Color legend: −log₁₀(FDR) gradient
    const colorLegY = 18 + sizeLevels.length * 24 + 14;
    addText(g, legX, colorLegY, '−log\u2081\u2080(FDR)', {
        size: '9px', fill: textColor, anchor: 'start', weight: '600',
        family: "'EB Garamond', Georgia, serif"
    });

    const gradH = 60;
    const gradW = 12;
    const gradY = colorLegY + 10;
    const gradSteps = 15;
    for (let i = 0; i < gradSteps; i++) {
        const frac = i / (gradSteps - 1);
        const cy = gradY + frac * gradH;
        const colorFn = PALETTES[palette] || PALETTES['Default'];
        const color = colorFn(1 - frac, theme);
        addRect(g, legX, cy, gradW, gradH / gradSteps + 1, color);
    }
    addRect(g, legX, gradY, gradW, gradH, 'none').setAttribute('stroke', axisColor);

    addText(g, legX + gradW + 4, gradY + 4, maxLogFDR.toFixed(1), {
        size: '8px', fill: textMuted, anchor: 'start', baseline: 'middle',
        family: "'EB Garamond', Georgia, serif"
    });
    addText(g, legX + gradW + 4, gradY + gradH, '0.0', {
        size: '8px', fill: textMuted, anchor: 'start', baseline: 'middle',
        family: "'EB Garamond', Georgia, serif"
    });

    return svg;
}

// ===== SVG Helpers =====

function makeSVG(w, h) {
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    svg.setAttribute('width', w);
    svg.setAttribute('height', h);
    svg.setAttribute('viewBox', `0 0 ${w} ${h}`);
    return svg;
}

function addGroup(parent, tx, ty) {
    const g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
    g.setAttribute('transform', `translate(${tx},${ty})`);
    parent.appendChild(g);
    return g;
}

function addRect(parent, x, y, w, h, fill, className) {
    const r = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    r.setAttribute('x', x); r.setAttribute('y', y);
    r.setAttribute('width', w); r.setAttribute('height', h);
    r.setAttribute('fill', fill);
    if (className) r.classList.add(className);
    parent.appendChild(r);
    return r;
}

function addLine(parent, x1, y1, x2, y2, stroke, width, dash) {
    const l = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    l.setAttribute('x1', x1); l.setAttribute('y1', y1);
    l.setAttribute('x2', x2); l.setAttribute('y2', y2);
    l.setAttribute('stroke', stroke);
    l.setAttribute('stroke-width', width || 1);
    if (dash) l.setAttribute('stroke-dasharray', dash);
    parent.appendChild(l);
    return l;
}

function addCircle(parent, cx, cy, r, fill) {
    const c = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
    c.setAttribute('cx', cx); c.setAttribute('cy', cy);
    c.setAttribute('r', r); c.setAttribute('fill', fill);
    parent.appendChild(c);
    return c;
}

function addText(parent, x, y, content, opts = {}) {
    const t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    t.setAttribute('x', x); t.setAttribute('y', y);
    t.setAttribute('fill', opts.fill || '#000');
    t.setAttribute('font-size', opts.size || '12px');
    t.setAttribute('font-family', opts.family || "'EB Garamond', Georgia, serif");
    t.setAttribute('text-anchor', opts.anchor || 'start');
    if (opts.weight) t.setAttribute('font-weight', opts.weight);
    if (opts.baseline) t.setAttribute('dominant-baseline', opts.baseline);
    t.textContent = content;
    parent.appendChild(t);
    return t;
}

// ===== Axis Utilities =====

function niceMax(val) {
    if (val <= 0) return 1;
    const mag = Math.pow(10, Math.floor(Math.log10(val)));
    const norm = val / mag;
    if (norm <= 1) return mag;
    if (norm <= 2) return 2 * mag;
    if (norm <= 5) return 5 * mag;
    return 10 * mag;
}

function niceTicksFor(min, max, count) {
    if (max <= min) return [0];
    const step = (max - min) / count;
    const ticks = [];
    for (let i = 0; i <= count; i++) {
        ticks.push(Math.round((min + step * i) * 1000) / 1000);
    }
    return ticks;
}

function formatTick(v) {
    if (v === 0) return '0';
    if (Number.isInteger(v)) return String(v);
    return v.toFixed(1);
}

function truncLabel(str, max) {
    if (!str) return '';
    str = str.charAt(0).toUpperCase() + str.slice(1);
    if (str.length <= max) return str;
    return str.substring(0, max - 1) + '…';
}

// ===== Hierarchical Clustering =====

/**
 * Compute pairwise Jaccard distance matrix between enrichment terms based on gene overlap.
 * Jaccard distance = 1 - |A ∩ B| / |A ∪ B|
 */
function computeJaccardDistanceMatrix(data) {
    const n = data.length;
    const geneSets = data.map(d => new Set(d.genes || []));
    const dist = Array.from({ length: n }, () => new Float64Array(n));

    for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
            const a = geneSets[i];
            const b = geneSets[j];
            let intersection = 0;
            for (const g of a) {
                if (b.has(g)) intersection++;
            }
            const union = a.size + b.size - intersection;
            const d = union === 0 ? 1 : 1 - intersection / union;
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }
    return dist;
}

/**
 * UPGMA (Unweighted Pair Group Method with Arithmetic Mean) agglomerative clustering.
 * Returns a tree: internal nodes have { left, right, height }, leaves have { index, height: 0 }.
 */
function upgmaClustering(distMatrix) {
    const n = distMatrix.length;
    if (n === 0) return null;
    if (n === 1) return { index: 0, height: 0 };

    // Working copy of distances + cluster sizes
    const dist = distMatrix.map(row => Array.from(row));
    const sizes = new Array(n).fill(1);
    const active = new Set(Array.from({ length: n }, (_, i) => i));
    const nodes = Array.from({ length: n }, (_, i) => ({ index: i, height: 0 }));

    for (let step = 0; step < n - 1; step++) {
        // Find closest pair
        let minDist = Infinity, mi = -1, mj = -1;
        const activeArr = [...active];
        for (let ai = 0; ai < activeArr.length; ai++) {
            for (let aj = ai + 1; aj < activeArr.length; aj++) {
                const i = activeArr[ai], j = activeArr[aj];
                if (dist[i][j] < minDist) {
                    minDist = dist[i][j];
                    mi = i; mj = j;
                }
            }
        }

        // Create merged node
        const mergedNode = {
            left: nodes[mi],
            right: nodes[mj],
            height: minDist / 2
        };

        // Update distances (UPGMA: weighted average by cluster size)
        const newIdx = mi;
        const si = sizes[mi], sj = sizes[mj];
        for (const k of active) {
            if (k === mi || k === mj) continue;
            dist[newIdx][k] = (dist[mi][k] * si + dist[mj][k] * sj) / (si + sj);
            dist[k][newIdx] = dist[newIdx][k];
        }

        sizes[newIdx] = si + sj;
        nodes[newIdx] = mergedNode;
        active.delete(mj);
    }

    return nodes[[...active][0]];
}

/**
 * Assign y-positions to leaves (in-order traversal) and compute x from height.
 * Returns { leaves: [{index, x, y}], internals: [{x, yTop, yBottom}] }
 */
function layoutDendrogram(root, maxHeight, plotW, rowH) {
    const leaves = [];
    const branches = [];
    let leafIndex = 0;

    // x-scale: height → horizontal position (0 = leaves on right, maxHeight = root on left)
    const xScale = (h) => plotW * (1 - h / Math.max(maxHeight, 1e-10));

    function traverse(node) {
        if (node.index !== undefined && !node.left) {
            // Leaf
            const y = leafIndex * rowH + rowH / 2;
            leafIndex++;
            leaves.push({ index: node.index, x: xScale(0), y });
            return y;
        }

        const yLeft = traverse(node.left);
        const yRight = traverse(node.right);
        const x = xScale(node.height);

        // Horizontal lines from children to this node's x, then vertical connector
        branches.push({
            x,
            xLeft: xScale(node.left.height !== undefined ? node.left.height : 0),
            xRight: xScale(node.right.height !== undefined ? node.right.height : 0),
            yTop: Math.min(yLeft, yRight),
            yBottom: Math.max(yLeft, yRight),
            yLeft,
            yRight,
            height: node.height
        });

        return (yLeft + yRight) / 2;
    }

    traverse(root);
    return { leaves, branches };
}

/**
 * Create a hierarchical clustering dendrogram of enriched terms.
 * Terms are clustered by Jaccard similarity of their gene sets.
 */
function createClusterTree(results, topN = 20, palette = 'Default', title = 'Enrichment Clustering') {
    // Filter to terms that have at least one gene (required for Jaccard distance)
    const data = results.filter(r => r.fdr <= 1 && r.genes && r.genes.length > 0).slice(0, topN);
    if (data.length < 2) return null;

    const theme = document.documentElement.getAttribute('data-theme');
    const textColor = theme === 'dark' ? '#d4d4d4' : '#1a1a1a';
    const textMuted = theme === 'dark' ? '#888888' : '#666666';
    const axisColor = theme === 'dark' ? '#555555' : '#333333';
    const gridColor = theme === 'dark' ? '#333333' : '#e0e0e0';
    const bgColor = theme === 'dark' ? '#1a1a1a' : '#ffffff';

    // Clustering using Jaccard distance on gene sets
    const distMatrix = computeJaccardDistanceMatrix(data);
    const tree = upgmaClustering(distMatrix);
    if (!tree) return null;

    // Find max tree height for scaling
    function getMaxHeight(node) {
        if (!node.left) return 0;
        return Math.max(node.height || 0, getMaxHeight(node.left), getMaxHeight(node.right));
    }
    const maxTreeHeight = getMaxHeight(tree);

    // Layout
    const margin = { top: 50, right: 20, bottom: 40, left: 280, treeRight: 300 };
    const rowH = 22;
    const plotW = 300; // dendrogram width
    const plotH = data.length * rowH;
    const width = margin.left + plotW + margin.treeRight + margin.right;
    const height = margin.top + plotH + margin.bottom;

    const layout = layoutDendrogram(tree, maxTreeHeight, plotW, rowH);

    // Color scale: −log₁₀(FDR) → intensity
    const maxLogFDR = Math.max(...data.map(d => -Math.log10(Math.max(d.fdr, 1e-300))));
    const colorFn = PALETTES[palette] || PALETTES['Default'];
    const fdrColor = (fdr) => {
        const logFDR = -Math.log10(Math.max(fdr, 1e-300));
        const t = maxLogFDR > 0 ? Math.min(logFDR / maxLogFDR, 1) : 0.5;
        return colorFn(t, theme);
    };

    const svg = makeSVG(width, height);
    addRect(svg, 0, 0, width, height, bgColor, 'plot-bg');

    // Title
    addText(svg, width / 2, 22, title, {
        size: '14px', weight: '700', fill: textColor, anchor: 'middle',
        family: "'EB Garamond', Georgia, serif"
    });
    const sigCount = data.filter(d => d.fdr < 0.05).length;
    addText(svg, width / 2, 38, `Top ${data.length} terms · ${sigCount} significant (FDR < 0.05) · Clustered by gene overlap (Jaccard)`, {
        size: '10px', fill: textMuted, anchor: 'middle',
        family: "'EB Garamond', Georgia, serif"
    });

    const g = addGroup(svg, margin.left, margin.top);

    // Draw branches
    for (const b of layout.branches) {
        // Vertical connector
        addLine(g, b.x, b.yTop, b.x, b.yBottom, axisColor, 1);
        // Horizontal line to left child
        addLine(g, b.x, b.yLeft, b.xLeft, b.yLeft, axisColor, 1);
        // Horizontal line to right child
        addLine(g, b.x, b.yRight, b.xRight, b.yRight, axisColor, 1);
    }

    // Draw leaves: colored dot + term label
    for (const leaf of layout.leaves) {
        const d = data[leaf.index];
        const color = fdrColor(d.fdr);

        // Colored dot at leaf
        const c = addCircle(g, leaf.x + 4, leaf.y, 5, color);
        c.setAttribute('stroke', axisColor);
        c.setAttribute('stroke-width', '0.5');

        // Gene count badge
        addText(g, leaf.x + 14, leaf.y + 1, `${d.geneCount}`, {
            size: '8px', fill: textMuted, anchor: 'start', baseline: 'middle',
            family: "'EB Garamond', Georgia, serif"
        });

        // Term description
        const label = truncLabel(d.description, 50);
        addText(g, leaf.x + 30, leaf.y + 1, label, {
            size: '10px', fill: textColor, anchor: 'start', baseline: 'middle',
            family: "'EB Garamond', Georgia, serif"
        });

        // Significance marker
        if (d.fdr < 0.001) {
            addText(g, leaf.x + 30 + label.length * 5.5 + 4, leaf.y + 1, '***', {
                size: '8px', fill: textMuted, anchor: 'start', baseline: 'middle'
            });
        } else if (d.fdr < 0.01) {
            addText(g, leaf.x + 30 + label.length * 5.5 + 4, leaf.y + 1, '**', {
                size: '8px', fill: textMuted, anchor: 'start', baseline: 'middle'
            });
        } else if (d.fdr < 0.05) {
            addText(g, leaf.x + 30 + label.length * 5.5 + 4, leaf.y + 1, '*', {
                size: '8px', fill: textMuted, anchor: 'start', baseline: 'middle'
            });
        }
    }

    // Distance axis at bottom
    const distTicks = niceTicksFor(0, maxTreeHeight, 4);
    for (const t of distTicks) {
        const x = plotW * (1 - t / Math.max(maxTreeHeight, 1e-10));
        addLine(g, x, plotH + 3, x, plotH + 8, axisColor, 1);
        addText(g, x, plotH + 20, t.toFixed(2), {
            size: '9px', fill: textColor, anchor: 'middle',
            family: "'EB Garamond', Georgia, serif"
        });
    }
    addLine(g, 0, plotH + 3, plotW, plotH + 3, axisColor, 1);
    addText(svg, margin.left + plotW / 2, height - 6, 'Jaccard Distance', {
        size: '10px', fill: textColor, anchor: 'middle', weight: '500',
        family: "'EB Garamond', Georgia, serif"
    });

    // FDR color legend (top-right area)
    const legX = plotW + 240;
    const legY = 0;
    const gradH = Math.min(plotH * 0.4, 80);
    const gradW = 12;

    addText(g, legX, legY - 6, '−log\u2081\u2080(FDR)', {
        size: '9px', fill: textColor, anchor: 'start', weight: '600',
        family: "'EB Garamond', Georgia, serif"
    });

    const gradSteps = 15;
    for (let i = 0; i < gradSteps; i++) {
        const frac = i / (gradSteps - 1);
        const cy = legY + frac * gradH;
        const color = colorFn(1 - frac, theme);
        addRect(g, legX, cy, gradW, gradH / gradSteps + 1, color);
    }
    addRect(g, legX, legY, gradW, gradH, 'none').setAttribute('stroke', axisColor);

    addText(g, legX + gradW + 4, legY + 4, maxLogFDR.toFixed(1), {
        size: '8px', fill: textMuted, anchor: 'start', baseline: 'middle',
        family: "'EB Garamond', Georgia, serif"
    });
    addText(g, legX + gradW + 4, legY + gradH, '0.0', {
        size: '8px', fill: textMuted, anchor: 'start', baseline: 'middle',
        family: "'EB Garamond', Georgia, serif"
    });

    return svg;
}

window.Plots = { createBarChart, createDotPlot, createClusterTree, PALETTES };
