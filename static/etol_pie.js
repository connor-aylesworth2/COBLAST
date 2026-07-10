/*
 * eToL / eToL-V domain pie chart.
 *
 * A coarser companion to the heatmap: what share of a job's non-human
 * probe-matched reads falls into each Tree-of-Life domain (Archaea, Bacteria,
 * Fungi, ...). Reuses the same /batch-results/<id>/etol-matrix.json the heatmap
 * loads -- no extra endpoint -- and draws an inline SVG pie so the standalone
 * build stays lean and works offline.
 *
 * The matrix's hit counts are already the non-human reads: the human PGK1/hNSE
 * control probes are not panel records (so no "Human control" row exists), and
 * when the secondary human filter ran, human-derived reads were dropped from the
 * hit list upstream.
 */
(function () {
  "use strict";

  var root = document.getElementById("etol-pie");
  if (!root) {
    return;
  }

  var MATRIX_URL = root.dataset.matrixUrl;

  // Stable colours for the known ToL domains; anything else draws from a palette.
  var DOMAIN_COLORS = {
    "Archaea": "#8c564b",
    "Bacteria": "#1f77b4",
    "Chloroplastida": "#2ca02c",
    "Amoebozoa": "#9467bd",
    "Basal Eukaryota": "#17becf",
    "Fungi": "#ff7f0e",
    "Holozoa/Metazoa": "#d62728",
    "Viruses": "#e377c2",
    "Other": "#7f7f7f",
  };
  var PALETTE = ["#bcbd22", "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173"];
  var extraColors = {};

  function colorFor(domain) {
    if (domain in DOMAIN_COLORS) return DOMAIN_COLORS[domain];
    if (!(domain in extraColors)) {
      extraColors[domain] = PALETTE[Object.keys(extraColors).length % PALETTE.length];
    }
    return extraColors[domain];
  }

  function esc(text) {
    return String(text)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  }

  function setStatus(message) {
    root.querySelector("[data-role=status]").textContent = message || "";
  }

  var matrix = null;

  // Sum matched reads per domain over one sample column (colIndex >= 0) or the
  // whole job (colIndex < 0). Zero-total domains are omitted.
  function domainTotals(colIndex) {
    var totals = {};
    matrix.rows.forEach(function (row, r) {
      var sum = 0;
      if (colIndex < 0) {
        for (var c = 0; c < matrix.cols.length; c += 1) sum += matrix.hits[r][c] || 0;
      } else {
        sum = matrix.hits[r][colIndex] || 0;
      }
      if (sum > 0) totals[row.domain] = (totals[row.domain] || 0) + sum;
    });
    return totals;
  }

  function polar(cx, cy, radius, angle) {
    return [cx + radius * Math.cos(angle), cy + radius * Math.sin(angle)];
  }

  function buildSvg(colIndex) {
    var totals = domainTotals(colIndex);
    var domains = Object.keys(totals).sort(function (a, b) { return totals[b] - totals[a]; });
    var grand = domains.reduce(function (s, d) { return s + totals[d]; }, 0);

    var W = 540, H = 340, cx = 160, cy = 178, R = 130;
    var parts = [
      '<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H +
        '" viewBox="0 0 ' + W + " " + H + '" font-family="system-ui, sans-serif">',
      '<rect width="' + W + '" height="' + H + '" fill="#ffffff"/>',
    ];
    var scope = colIndex < 0 ? "all samples (whole job)" : matrix.cols[colIndex].sample;
    var title = (matrix.preset_label || "eToL") + " — non-human reads by domain";
    parts.push('<text x="16" y="24" font-size="14" font-weight="600">' + esc(title) + "</text>");
    parts.push('<text x="16" y="42" font-size="11" fill="#555">' + esc(scope) + "</text>");

    if (grand <= 0) {
      parts.push('<text x="16" y="80" font-size="12" fill="#666">No non-human probe-matched reads to plot.</text></svg>');
      return parts.join("");
    }

    var angle = -Math.PI / 2; // start at 12 o'clock
    domains.forEach(function (d) {
      var frac = totals[d] / grand;
      var end = angle + frac * 2 * Math.PI;
      var fill = colorFor(d);
      var tip = esc(d + ": " + totals[d] + " reads (" + (frac * 100).toFixed(1) + "%)");
      if (domains.length === 1) {
        // A single slice is a full circle -- one arc can't close a 360° sweep.
        parts.push('<circle cx="' + cx + '" cy="' + cy + '" r="' + R + '" fill="' + fill +
          '"><title>' + tip + "</title></circle>");
      } else {
        var p0 = polar(cx, cy, R, angle);
        var p1 = polar(cx, cy, R, end);
        var large = frac > 0.5 ? 1 : 0;
        parts.push(
          '<path d="M' + cx + " " + cy + " L" + p0[0].toFixed(2) + " " + p0[1].toFixed(2) +
            " A" + R + " " + R + " 0 " + large + " 1 " + p1[0].toFixed(2) + " " + p1[1].toFixed(2) +
            ' Z" fill="' + fill + '" stroke="#ffffff" stroke-width="1"><title>' + tip + "</title></path>"
        );
      }
      angle = end;
    });

    // Legend on the right: domain — percent (read count), largest first.
    var lx = 320, ly = 66;
    domains.forEach(function (d, i) {
      var y = ly + i * 22;
      var pct = (totals[d] / grand * 100).toFixed(1);
      parts.push('<rect x="' + lx + '" y="' + (y - 11) + '" width="13" height="13" fill="' + colorFor(d) + '"/>');
      parts.push('<text x="' + (lx + 19) + '" y="' + y + '" font-size="11" fill="#333">' +
        esc(d + " — " + pct + "% (" + totals[d] + ")") + "</text>");
    });
    parts.push('<text x="16" y="' + (H - 12) + '" font-size="10" fill="#777">Total non-human reads: ' + grand + "</text>");
    parts.push("</svg>");
    return parts.join("");
  }

  function render() {
    var sel = root.querySelector("[data-role=sample]");
    var colIndex = sel ? parseInt(sel.value, 10) : -1;
    if (isNaN(colIndex)) colIndex = -1;
    root.querySelector("[data-role=canvas]").innerHTML = buildSvg(colIndex);
  }

  function populateSamples() {
    var sel = root.querySelector("[data-role=sample]");
    if (!sel) return;
    var opts = ['<option value="-1">All samples (whole job)</option>'];
    matrix.cols.forEach(function (col, i) {
      opts.push('<option value="' + i + '">' + esc(col.sample) + "</option>");
    });
    sel.innerHTML = opts.join("");
  }

  function currentSvg() {
    var svg = root.querySelector("[data-role=canvas] svg");
    return svg ? svg.outerHTML : null;
  }

  function download(name, mime, data) {
    var url = URL.createObjectURL(new Blob([data], { type: mime }));
    var a = document.createElement("a");
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  function exportPng() {
    var svgEl = root.querySelector("[data-role=canvas] svg");
    if (!svgEl) return;
    var w = parseInt(svgEl.getAttribute("width"), 10);
    var h = parseInt(svgEl.getAttribute("height"), 10);
    var img = new Image();
    var url = URL.createObjectURL(new Blob([svgEl.outerHTML], { type: "image/svg+xml" }));
    img.onload = function () {
      var canvas = document.createElement("canvas");
      canvas.width = w * 2;
      canvas.height = h * 2;
      var ctx = canvas.getContext("2d");
      ctx.scale(2, 2);
      ctx.drawImage(img, 0, 0);
      URL.revokeObjectURL(url);
      canvas.toBlob(function (blob) {
        var a = document.createElement("a");
        a.href = URL.createObjectURL(blob);
        a.download = "etol_domain_pie.png";
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(a.href);
      });
    };
    img.src = url;
  }

  root.querySelector("[data-role=png]").addEventListener("click", exportPng);
  root.querySelector("[data-role=svg]").addEventListener("click", function () {
    var svgText = currentSvg();
    if (svgText) download("etol_domain_pie.svg", "image/svg+xml", svgText);
  });
  var sampleSel = root.querySelector("[data-role=sample]");
  if (sampleSel) sampleSel.addEventListener("change", render);

  setStatus("Loading…");
  fetch(MATRIX_URL + "?level=species")
    .then(function (resp) {
      if (!resp.ok) throw new Error("HTTP " + resp.status);
      return resp.json();
    })
    .then(function (data) {
      matrix = data;
      populateSamples();
      setStatus("");
      render();
    })
    .catch(function (err) {
      setStatus("Failed to load pie chart: " + err.message);
    });
})();
