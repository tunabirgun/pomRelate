/**
 * pomRelate — PPI Network Visualization
 * Force-directed layout with zoom, pan, and drag.
 */

/**
 * Build and render a PPI network for the given genes.
 * @returns {Object|null} { svg, nodes }
 */
function buildPPINetworkSVG(resolvedGenes, ppiData, infoData, scoreThreshold, getNameFn, taxid) {
    if (!ppiData) return null;

    const queryIds = new Set(resolvedGenes.filter(g => g.proteinId).map(g => g.proteinId));
    if (queryIds.size === 0) return null;

    // Collect nodes and edges
    const nodeMap = {};
    const edges = [];
    const edgeSet = new Set();

    for (const pid of queryIds) {
        nodeMap[pid] = { id: pid, name: getNameFn(pid), isQuery: true, degree: 0 };
    }

    for (const pid of queryIds) {
        const interactions = ppiData[pid];
        if (!interactions) continue;

        const top = interactions
            .filter(i => i.s >= scoreThreshold)
            .sort((a, b) => b.s - a.s)
            .slice(0, 10);

        for (const { p, s } of top) {
            if (!nodeMap[p]) {
                nodeMap[p] = { id: p, name: getNameFn(p), isQuery: queryIds.has(p), degree: 0 };
            }
            const ek = [pid, p].sort().join('|');
            if (!edgeSet.has(ek)) {
                edgeSet.add(ek);
                edges.push({ source: pid, target: p, score: s });
                nodeMap[pid].degree++;
                nodeMap[p].degree++;
            }
        }
    }

    // Cross-query edges (check both directions since PPI data may be asymmetric)
    const qa = [...queryIds];
    for (let i = 0; i < qa.length; i++) {
        for (let j = i + 1; j < qa.length; j++) {
            const ek = [qa[i], qa[j]].sort().join('|');
            if (edgeSet.has(ek)) continue;
            // Check A→B
            let link = null;
            const intsA = ppiData[qa[i]];
            if (intsA) link = intsA.find(x => x.p === qa[j] && x.s >= scoreThreshold);
            // Check B→A if not found
            if (!link) {
                const intsB = ppiData[qa[j]];
                if (intsB) link = intsB.find(x => x.p === qa[i] && x.s >= scoreThreshold);
            }
            if (link) {
                edgeSet.add(ek);
                edges.push({ source: qa[i], target: qa[j], score: link.s });
                nodeMap[qa[i]].degree++;
                nodeMap[qa[j]].degree++;
            }
        }
    }

    const nodes = Object.values(nodeMap);
    if (nodes.length === 0) return null;

    // Hub detection
    const sorted = [...nodes].sort((a, b) => b.degree - a.degree);
    const hubTh = Math.max(sorted[Math.floor(nodes.length * 0.2)]?.degree || 1, 3);
    for (const n of nodes) n.isHub = n.degree >= hubTh;

    // Initialize random positions
    const width = 800, height = 600;
    const cx = width / 2, cy = height / 2;
    for (const n of nodes) {
        n.x = cx + (Math.random() - 0.5) * 50;
        n.y = cy + (Math.random() - 0.5) * 50;
        n.vx = 0; n.vy = 0;
        n.r = n.isQuery || n.isHub ? 14 : 7;
    }

    runSimulationSync(nodes, edges, width, height);

    // Calculate average interaction confidence per node
    nodes.forEach(n => {
        const incident = edges.filter(e => e.source === n.id || e.target === n.id);
        n.scoreSum = incident.length > 0
            ? Math.round(incident.reduce((sum, e) => sum + e.score, 0) / incident.length)
            : 0;
    });

    const svg = renderNetworkViewer(nodes, edges, width, height, taxid);
    return { svg, nodes };
}

/** Run force-directed simulation synchronously. */
function runSimulationSync(nodes, edges, width, height) {
    const k = Math.sqrt((width * height) / (nodes.length || 1)) * 1.5;
    let alpha = 1.0;
    const idxMap = {};
    nodes.forEach((n, i) => idxMap[n.id] = i);
    const cx = width / 2, cy = height / 2;

    for (let iter = 0; iter < 300; iter++) {
        alpha *= 0.99;
        if (alpha < 0.001) break;

        // Reset forces
        nodes.forEach(n => { n.vx = 0; n.vy = 0; });

        // Repulsion
        for (let i = 0; i < nodes.length; i++) {
            for (let j = i + 1; j < nodes.length; j++) {
                const dx = nodes[i].x - nodes[j].x;
                const dy = nodes[i].y - nodes[j].y;
                let d2 = dx * dx + dy * dy;
                if (d2 === 0) d2 = 0.1;
                const dist = Math.sqrt(d2);

                const f = (k * k) / dist;
                const fx = (dx / dist) * f * alpha;
                const fy = (dy / dist) * f * alpha;

                nodes[i].vx += fx; nodes[i].vy += fy;
                nodes[j].vx -= fx; nodes[j].vy -= fy;

                // Collision
                const rSum = nodes[i].r + nodes[j].r + 4;
                if (dist < rSum) {
                    const overlap = rSum - dist;
                    const pushX = (dx / dist) * overlap * 0.8 * alpha;
                    const pushY = (dy / dist) * overlap * 0.8 * alpha;
                    nodes[i].x += pushX; nodes[i].y += pushY;
                    nodes[j].x -= pushX; nodes[j].y -= pushY;
                }
            }
        }

        // Attraction
        edges.forEach(e => {
            const s = nodes[idxMap[e.source]];
            const t = nodes[idxMap[e.target]];
            if (!s || !t) return;

            const dx = s.x - t.x;
            const dy = s.y - t.y;
            const dist = Math.sqrt(dx * dx + dy * dy) || 1;

            const f = (dist * dist) / k;
            const strength = 0.5 + (e.score / 1000);

            const fx = (dx / dist) * f * strength * alpha;
            const fy = (dy / dist) * f * strength * alpha;

            s.vx -= fx; s.vy -= fy;
            t.vx += fx; t.vy += fy;
        });

        // Center Gravity
        nodes.forEach(n => {
            n.vx += (cx - n.x) * 0.02 * alpha;
            n.vy += (cy - n.y) * 0.02 * alpha;
        });

        // Update positions
        nodes.forEach(n => {
            const v = Math.sqrt(n.vx * n.vx + n.vy * n.vy) || 1;
            const maxV = 15 * alpha;
            if (v > maxV) { n.vx *= maxV / v; n.vy *= maxV / v; }
            n.x += n.vx;
            n.y += n.vy;
        });
    }
}


/**
 * Render the viewer with Zoom/Pan capabilities.
 */
function renderNetworkViewer(nodes, edges, width, height, taxid) {
    const theme = document.documentElement.getAttribute('data-theme');
    const isDark = theme === 'dark';

    const uid = Math.random().toString(36).substr(2, 6);
    const gid = (name) => `${name}-${uid}`;

    // SVG Setup
    const svgW = width;
    const svgH = height;

    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('width', '100%');
    svg.setAttribute('height', '100%');
    svg.setAttribute('viewBox', `0 0 ${svgW} ${svgH}`);
    svg.style.border = isDark ? '1px solid #333' : '1px solid #ddd';
    svg.style.borderRadius = '4px';
    svg.style.background = isDark ? '#1a1a1a' : '#ffffff';
    svg.style.cursor = 'default';
    svg.style.userSelect = 'none';

    // Defs & Gradients
    const defs = document.createElementNS('http://www.w3.org/2000/svg', 'defs');
    function createGrad(id, colorCenter, colorEdge) {
        const rg = document.createElementNS('http://www.w3.org/2000/svg', 'radialGradient');
        rg.setAttribute('id', id);
        rg.setAttribute('cx', '35%'); rg.setAttribute('cy', '35%'); rg.setAttribute('r', '65%');
        const s1 = document.createElementNS('http://www.w3.org/2000/svg', 'stop');
        s1.setAttribute('offset', '0%'); s1.setAttribute('stop-color', colorCenter);
        const s2 = document.createElementNS('http://www.w3.org/2000/svg', 'stop');
        s2.setAttribute('offset', '100%'); s2.setAttribute('stop-color', colorEdge);
        rg.appendChild(s1); rg.appendChild(s2);
        return rg;
    }
    defs.appendChild(createGrad(gid('grad-query'), '#ff9999', '#c92a2a'));
    defs.appendChild(createGrad(gid('grad-hub'), '#74c0fc', '#1864ab'));
    defs.appendChild(createGrad(gid('grad-node'), '#f8f9fa', '#868e96'));
    defs.appendChild(createGrad(gid('grad-node-dark'), '#495057', '#212529'));
    svg.appendChild(defs);

    // Container Group
    const container = document.createElementNS('http://www.w3.org/2000/svg', 'g');
    svg.appendChild(container);

    // Auto-Fit Logic
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    nodes.forEach(n => {
        if (n.x < minX) minX = n.x;
        if (n.x > maxX) maxX = n.x;
        if (n.y < minY) minY = n.y;
        if (n.y > maxY) maxY = n.y;
    });

    // Add padding (node radius ~14)
    minX -= 40; maxX += 40; minY -= 40; maxY += 40;

    const contentW = maxX - minX;
    const contentH = maxY - minY;
    let scale = Math.min(svgW / contentW, svgH / contentH);
    if (scale > 2) scale = 2;

    const midX = (minX + maxX) / 2;
    const midY = (minY + maxY) / 2;
    let translateX = (svgW / 2) - (midX * scale);
    let translateY = (svgH / 2) - (midY * scale);

    container.setAttribute('transform', `translate(${translateX},${translateY}) scale(${scale})`);

    // Draw Elements
    const sortedNodes = [...nodes].sort((a, b) => (a.isQuery ? 1 : 0) - (b.isQuery ? 1 : 0));
    const nodeElements = {};
    const edgeElements = [];

    // Edges — normalize visual weight to absolute STRING score range [0, 1000]
    for (const e of edges) {
        const scoreNorm = e.score / 1000;
        const width = 0.5 + scoreNorm * 2.5;
        const opacity = 0.2 + scoreNorm * 0.5;

        const l = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        l.setAttribute('stroke', isDark ? '#888' : '#666');
        l.setAttribute('stroke-width', width);
        l.setAttribute('stroke-opacity', opacity);

        const s = nodes.find(n => n.id === e.source);
        const t = nodes.find(n => n.id === e.target);
        if (s && t) {
            l.setAttribute('x1', s.x); l.setAttribute('y1', s.y);
            l.setAttribute('x2', t.x); l.setAttribute('y2', t.y);
        }
        container.appendChild(l);
        edgeElements.push({ source: e.source, target: e.target, el: l });
    }

    // Nodes
    for (const n of sortedNodes) {
        const g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        g.setAttribute('transform', `translate(${n.x},${n.y})`);
        g.style.cursor = 'grab';

        const c = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
        c.setAttribute('r', n.r);

        if (n.isQuery) c.setAttribute('fill', `url(#${gid('grad-query')})`);
        else if (n.isHub) c.setAttribute('fill', `url(#${gid('grad-hub')})`);
        else c.setAttribute('fill', isDark ? `url(#${gid('grad-node-dark')})` : `url(#${gid('grad-node')})`);

        c.setAttribute('stroke', isDark ? '#fff' : '#333');
        c.setAttribute('stroke-width', n.isQuery ? 1.5 : 0.5);
        c.setAttribute('stroke-opacity', 0.8);
        g.appendChild(c);

        if (n.isQuery || n.isHub) {
            const t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            t.textContent = n.name;
            t.setAttribute('y', n.r + 12);
            t.setAttribute('text-anchor', 'middle');
            t.setAttribute('fill', isDark ? '#eee' : '#333');
            t.setAttribute('font-size', n.isQuery ? '12px' : '10px');
            t.setAttribute('font-family', 'sans-serif');
            t.setAttribute('font-weight', n.isQuery ? 'bold' : 'normal');
            t.style.pointerEvents = 'none';
            t.style.textShadow = isDark ? '0 1px 2px #000' : '0 1px 2px #fff';
            g.appendChild(t);
        }

        // Tooltip data attributes
        if (taxid) {
            g.dataset.pid = n.id;
            g.dataset.taxid = taxid;
            g.addEventListener('mouseenter', (e) => {
                if (window.showGeneTooltip) window.showGeneTooltip(n.id, taxid, e);
            });
            g.addEventListener('mouseleave', (e) => {
                // Don't hide if mouse is moving to the tooltip
                const related = e.relatedTarget;
                const tooltip = document.getElementById('gene-tooltip');
                if (related && tooltip && tooltip.contains(related)) return;
                if (window.hideGeneTooltip) window.hideGeneTooltip();
            });
        }

        container.appendChild(g);
        nodeElements[n.id] = g;
    }

    // Setup Interaction
    setupInteraction(svg, container, translateX, translateY, scale, nodes, nodeElements, edgeElements);

    return svg;
}

/**
 * Zoom, Pan, Drag
 */
function setupInteraction(svg, container, initX, initY, initScale, nodes, nodeElements, edgeElements) {
    let scale = initScale;
    let translateX = initX;
    let translateY = initY;

    function updateTransform() {
        container.setAttribute('transform', `translate(${translateX},${translateY}) scale(${scale})`);
    }

    // Zoom
    svg.addEventListener('wheel', (e) => {
        e.preventDefault();

        const rect = svg.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;

        const px = (mx - translateX) / scale;
        const py = (my - translateY) / scale;

        const delta = -e.deltaY;
        const factor = Math.pow(1.001, delta);

        const newScale = scale * factor;
        if (newScale < 0.1 || newScale > 10) return;
        scale = newScale;

        translateX = mx - px * scale;
        translateY = my - py * scale;

        updateTransform();
    });

    // Pan
    let isPanning = false;
    let startX = 0, startY = 0;

    svg.addEventListener('mousedown', (e) => {
        if (e.target !== svg) return;
        isPanning = true;
        startX = e.clientX - translateX;
        startY = e.clientY - translateY;
        svg.style.cursor = 'grabbing';
    });

    // Drag Node
    let draggedNode = null;
    let dragInteractOffsetX = 0;
    let dragInteractOffsetY = 0;

    Object.entries(nodeElements).forEach(([id, el]) => {
        el.addEventListener('mousedown', (e) => {
            e.stopPropagation();
            const n = nodes.find(x => x.id === id);
            if (!n) return;
            draggedNode = n;
            el.style.cursor = 'grabbing';

            const rect = svg.getBoundingClientRect();
            const mx = e.clientX - rect.left;
            const my = e.clientY - rect.top;

            const wx = (mx - translateX) / scale;
            const wy = (my - translateY) / scale;

            dragInteractOffsetX = n.x - wx;
            dragInteractOffsetY = n.y - wy;
        });
    });

    function onMove(e) {
        if (!svg.isConnected) {
            window.removeEventListener('mousemove', onMove);
            window.removeEventListener('mouseup', onUp);
            return;
        }

        const rect = svg.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;

        if (isPanning) {
            e.preventDefault();
            translateX = e.clientX - startX;
            translateY = e.clientY - startY;
            updateTransform();
        }
        else if (draggedNode) {
            e.preventDefault();

            // Mouse in world space
            const wx = (mx - translateX) / scale;
            const wy = (my - translateY) / scale;

            draggedNode.x = wx + dragInteractOffsetX;
            draggedNode.y = wy + dragInteractOffsetY;

            const el = nodeElements[draggedNode.id];
            if (el) el.setAttribute('transform', `translate(${draggedNode.x},${draggedNode.y})`);

            edgeElements.forEach(edge => {
                if (edge.source === draggedNode.id || edge.target === draggedNode.id) {
                    const s = nodes.find(n => n.id === edge.source);
                    const t = nodes.find(n => n.id === edge.target);
                    if (s && t) {
                        edge.el.setAttribute('x1', s.x); edge.el.setAttribute('y1', s.y);
                        edge.el.setAttribute('x2', t.x); edge.el.setAttribute('y2', t.y);
                    }
                }
            });
        }
    }

    function onUp() {
        isPanning = false;
        svg.style.cursor = 'default';

        if (draggedNode) {
            const el = nodeElements[draggedNode.id];
            if (el) el.style.cursor = 'grab';
            draggedNode = null;
        }

        if (!svg.isConnected) {
            window.removeEventListener('mousemove', onMove);
            window.removeEventListener('mouseup', onUp);
        }
    }

    window.addEventListener('mousemove', onMove);
    window.addEventListener('mouseup', onUp);
}

window.PPINetwork = { buildPPINetworkSVG };
