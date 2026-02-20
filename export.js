/**
 * pomRelate — Export Utilities
 * CSV download + SVG/PNG/PDF export for enrichment plots.
 */

// ===== CSV Export =====

/**
 * Download enrichment results as CSV.
 * @param {Array} results - enrichment result objects
 * @param {string} filename
 * @param {Function} getNameFn - function to get preferred name from protein ID
 */
function downloadCSV(results, filename, getNameFn) {
    const headers = ['Term', 'Description', 'Category', 'P-Value', 'FDR', 'Fold Enrichment', 'Gene Count', 'Background Count', 'Genes'];
    const rows = results.map(r => [
        r.term,
        `"${(r.description || '').replace(/"/g, '""')}"`,
        `"${(r.category || '').replace(/"/g, '""')}"`,
        r.pValue.toExponential(4),
        r.fdr.toExponential(4),
        r.fold,
        r.geneCount,
        r.bgCount,
        `"${r.genes.map(g => getNameFn ? getNameFn(g) : g).join(', ')}"`,
    ]);

    const csv = [headers.join(','), ...rows.map(r => r.join(','))].join('\n');
    triggerDownload(csv, filename, 'text/csv;charset=utf-8;');
}

// ===== SVG Export =====

function downloadSVG(svgElement, filename) {
    const serializer = new XMLSerializer();
    const clone = svgElement.cloneNode(true);
    injectComputedStyles(svgElement, clone);
    const svgStr = serializer.serializeToString(clone);
    const blob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
    triggerDownloadBlob(blob, filename);
}

// ===== PNG Export =====

function downloadPNG(svgElement, filename, scale = 4) {
    const serializer = new XMLSerializer();
    const clone = svgElement.cloneNode(true);
    injectComputedStyles(svgElement, clone);

    const rect = svgElement.getBoundingClientRect();
    let w = rect.width;
    let h = rect.height;

    const attrW = svgElement.getAttribute('width');
    const attrH = svgElement.getAttribute('height');
    if ((w === 0 || h === 0) && attrW && !attrW.includes('%')) w = parseFloat(attrW);
    if ((w === 0 || h === 0) && attrH && !attrH.includes('%')) h = parseFloat(attrH);

    if (w === 0) w = 800;
    if (h === 0) h = 600;

    clone.setAttribute('width', w);
    clone.setAttribute('height', h);

    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');

    const svgStr = serializer.serializeToString(clone);

    const canvas = document.createElement('canvas');
    canvas.width = w * scale;
    canvas.height = h * scale;

    const ctx = canvas.getContext('2d');
    ctx.scale(scale, scale);

    const bgColor = getComputedStyle(document.documentElement).getPropertyValue('--bg-panel').trim() || '#fff';
    ctx.fillStyle = bgColor;
    ctx.fillRect(0, 0, w, h);

    const img = new Image();
    const svgBlob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(svgBlob);

    img.onload = () => {
        try {
            ctx.drawImage(img, 0, 0, w, h);
            URL.revokeObjectURL(url);
            canvas.toBlob(blob => {
                if (blob) {
                    triggerDownloadBlob(blob, filename);
                } else {
                    console.error('Canvas to Blob failed');
                }
            }, 'image/png');
        } catch (e) {
            console.error('Error drawing image to canvas:', e);
        }
    };
    img.onerror = (e) => {
        console.error('Error loading SVG image for PNG export:', e);
        URL.revokeObjectURL(url);
    };
    img.src = url;
}

// ===== PDF Export =====

/**
 * Simple PDF export: renders SVG as a PNG image embedded in a minimal PDF.
 */
function downloadPDF(svgElement, filename, scale = 2) {
    const serializer = new XMLSerializer();
    const clone = svgElement.cloneNode(true);
    injectComputedStyles(svgElement, clone);
    const svgStr = serializer.serializeToString(clone);

    const canvas = document.createElement('canvas');
    // Use same robust dimension logic as PNG export
    const rect = svgElement.getBoundingClientRect();
    let w = rect.width;
    let h = rect.height;
    const attrW = svgElement.getAttribute('width');
    const attrH = svgElement.getAttribute('height');
    if ((w === 0 || h === 0) && attrW && !attrW.includes('%')) w = parseFloat(attrW);
    if ((w === 0 || h === 0) && attrH && !attrH.includes('%')) h = parseFloat(attrH);
    if (w === 0) w = 800;
    if (h === 0) h = 600;
    canvas.width = w * scale;
    canvas.height = h * scale;

    const ctx = canvas.getContext('2d');
    ctx.scale(scale, scale);

    const bgColor = getComputedStyle(document.documentElement).getPropertyValue('--bg-panel').trim() || '#fff';
    ctx.fillStyle = bgColor;
    ctx.fillRect(0, 0, w, h);

    const img = new Image();
    const svgBlob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
    const url = URL.createObjectURL(svgBlob);

    img.onload = () => {
        ctx.drawImage(img, 0, 0, w, h);
        URL.revokeObjectURL(url);

        const dataUrl = canvas.toDataURL('image/jpeg', 0.95);
        const imgData = atob(dataUrl.split(',')[1]);
        const imgBytes = new Uint8Array(imgData.length);
        for (let i = 0; i < imgData.length; i++) imgBytes[i] = imgData.charCodeAt(i);

        const pdfBytes = buildPDF(imgBytes, w * scale, h * scale, w, h);
        const blob = new Blob([pdfBytes], { type: 'application/pdf' });
        triggerDownloadBlob(blob, filename);
    };
    img.src = url;
}

/**
 * Build a minimal single-page PDF with an embedded JPEG image.
 */
function buildPDF(jpegBytes, imgW, imgH, pageW, pageH) {
    const margin = 36;
    const maxW = 595 - 2 * margin;
    const maxH = 842 - 2 * margin;
    const s = Math.min(maxW / pageW, maxH / pageH, 1);
    const dw = pageW * s;
    const dh = pageH * s;

    const offsets = [];
    let content = '';
    let pos = 0;

    function addObj(str) {
        offsets.push(pos);
        content += str;
        pos = content.length;
    }

    content += '%PDF-1.4\n';
    pos = content.length;

    addObj('1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n');
    addObj('2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n');
    addObj(`3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] /Contents 4 0 R /Resources << /XObject << /Img 5 0 R >> >> >>\nendobj\n`);
    const x = margin;
    const y = 842 - margin - dh;
    const stream = `q ${dw.toFixed(2)} 0 0 ${dh.toFixed(2)} ${x.toFixed(2)} ${y.toFixed(2)} cm /Img Do Q\n`;
    addObj(`4 0 obj\n<< /Length ${stream.length} >>\nstream\n${stream}endstream\nendobj\n`);

    const imgHeader = `5 0 obj\n<< /Type /XObject /Subtype /Image /Width ${imgW} /Height ${imgH} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /DCTDecode /Length ${jpegBytes.length} >>\nstream\n`;
    const imgFooter = `\nendstream\nendobj\n`;

    const encoder = new TextEncoder();
    const headerPart = encoder.encode(content + imgHeader);
    const footerPart = encoder.encode(imgFooter);

    // xref
    const imgObjOffset = pos;
    offsets.push(imgObjOffset);

    const xrefOffset = headerPart.length + jpegBytes.length + footerPart.length;

    let xref = `xref\n0 ${offsets.length + 1}\n`;
    xref += `0000000000 65535 f \n`;
    for (const off of offsets) {
        xref += `${String(off).padStart(10, '0')} 00000 n \n`;
    }
    // Adjust offset for image obj (it's in the binary part)
    // Actually we need to recalculate... let me simplify

    // Simpler approach: combine everything as text except the image stream
    const part1 = content + imgHeader;
    const part1Bytes = encoder.encode(part1);

    const afterImg = imgFooter;
    const afterImgBytes = encoder.encode(afterImg);

    // Recalculate offsets for image object
    const realImgOffset = part1Bytes.length - encoder.encode(imgHeader).length;

    // Build xref
    const totalXrefOffset = part1Bytes.length + jpegBytes.length + afterImgBytes.length;
    let xrefStr = `xref\n0 6\n`;
    xrefStr += `0000000000 65535 f \n`;
    for (let i = 0; i < offsets.length; i++) {
        let off = offsets[i];
        if (i === offsets.length - 1) {
            off = part1Bytes.length - encoder.encode(imgHeader).length;
        }
        xrefStr += `${String(off).padStart(10, '0')} 00000 n \n`;
    }
    xrefStr += `trailer\n<< /Size 6 /Root 1 0 R >>\nstartxref\n${totalXrefOffset}\n%%EOF\n`;

    const xrefBytes = encoder.encode(xrefStr);

    const total = new Uint8Array(part1Bytes.length + jpegBytes.length + afterImgBytes.length + xrefBytes.length);
    total.set(part1Bytes, 0);
    total.set(jpegBytes, part1Bytes.length);
    total.set(afterImgBytes, part1Bytes.length + jpegBytes.length);
    total.set(xrefBytes, part1Bytes.length + jpegBytes.length + afterImgBytes.length);

    return total;
}

// ===== Helpers =====

function triggerDownload(content, filename, mimeType) {
    const blob = new Blob([content], { type: mimeType });
    triggerDownloadBlob(blob, filename);
}

function triggerDownloadBlob(blob, filename) {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

/**
 * Inject computed styles from original SVG into cloned SVG for export.
 */
function injectComputedStyles(original, clone) {
    const origChildren = original.querySelectorAll('*');
    const cloneChildren = clone.querySelectorAll('*');

    for (let i = 0; i < origChildren.length && i < cloneChildren.length; i++) {
        const computed = getComputedStyle(origChildren[i]);
        const important = ['fill', 'stroke', 'stroke-width', 'stroke-dasharray', 'stroke-opacity', 'fill-opacity', 'font-family', 'font-size', 'font-weight', 'text-anchor', 'dominant-baseline', 'opacity', 'visibility'];
        for (const prop of important) {
            const val = computed.getPropertyValue(prop);
            if (val) cloneChildren[i].style.setProperty(prop, val);
        }
    }

    // Set SVG background
    const bg = getComputedStyle(document.documentElement).getPropertyValue('--bg-panel').trim() || '#fff';
    const rect = clone.querySelector('.plot-bg');
    if (rect) rect.setAttribute('fill', bg);
}

window.Export = { downloadCSV, downloadSVG, downloadPNG, downloadPDF };
