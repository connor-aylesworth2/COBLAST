/*
 * eToL / eToL-V domain pie chart.
 *
 * A coarser companion to the heatmap: what share of a job's non-human
 * probe-matched reads falls into each Tree-of-Life domain (Archaea, Bacteria,
 * Fungi, ...). Reuses the same /batch-results/<id>/etol-matrix.json the heatmap
 * loads -- no extra endpoint -- and draws an inline SVG pie so the standalone
 * build stays lean and works offline.
 *
 * Scopes: the whole job, one design-matrix condition (every sample in that group
 * pooled), one sample, or a grid of one labelled pie per condition -- the summary
 * of each sample group's microbiome. The condition labels are matrix.cols[i].condition,
 * already served for the heatmap's condition strip, so nothing new is fetched.
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
  var GRID = "grid"; // scope value for the by-condition contact sheet
  var TILE_W = 540, TILE_H = 340;

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

  // Quotes are escaped too: condition labels are free text from the design matrix
  // and land in option value="..." attributes as well as in text nodes.
  function esc(text) {
    return String(text)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function csvCell(value) {
    var s = String(value);
    return /[",\n\r]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  }

  function setStatus(message) {
    root.querySelector("[data-role=status]").textContent = message || "";
  }

  var matrix = null;

  // Sum matched reads per domain over a set of sample columns. One column, every
  // column, or a condition's columns all take this same path. Zero-total domains
  // are omitted.
  function domainTotals(colIdxs) {
    var totals = {};
    matrix.rows.forEach(function (row, r) {
      var sum = 0;
      colIdxs.forEach(function (c) {
        sum += matrix.hits[r][c] || 0;
      });
      if (sum > 0) totals[row.domain] = (totals[row.domain] || 0) + sum;
    });
    return totals;
  }

  // Distinct condition labels in column order. The server already grouped the
  // columns by condition (etol_summary.sort_results_by_condition), so first-seen
  // order here is the heatmap's left-to-right order.
  function conditionNames() {
    var names = [];
    matrix.cols.forEach(function (col) {
      var name = String(col.condition || "").trim();
      if (name && names.indexOf(name) < 0) names.push(name);
    });
    return names;
  }

  function conditionScope(name) {
    var cols = [];
    matrix.cols.forEach(function (col, i) {
      if (String(col.condition || "").trim() === name) cols.push(i);
    });
    return {
      value: "cond:" + name,
      label: name,
      subtitle: name + " — " + cols.length + (cols.length === 1 ? " sample" : " samples"),
      csvLabel: "Condition: " + name,
      cols: cols,
    };
  }

  // Every scope the pie can draw. Single source for the drop-down, render(), and
  // the CSV export, so the three cannot drift.
  function scopes() {
    var everyCol = matrix.cols.map(function (col, i) { return i; });
    var list = [{
      value: "all",
      label: "All samples (whole job)",
      subtitle: "all samples (whole job)",
      csvLabel: "All samples (whole job)",
      cols: everyCol,
    }];
    var conditions = conditionNames();
    if (conditions.length) {
      list.push({ value: GRID, label: "All conditions (grid)", subtitle: "", csvLabel: "", cols: everyCol });
      conditions.forEach(function (name) {
        list.push(conditionScope(name));
      });
    }
    matrix.cols.forEach(function (col, i) {
      list.push({
        value: "s:" + i,
        label: col.sample,
        subtitle: col.sample,
        csvLabel: col.sample,
        cols: [i],
      });
    });
    return list;
  }

  function currentScope() {
    var sel = root.querySelector("[data-role=sample]");
    var value = sel ? sel.value : "all";
    var list = scopes();
    for (var i = 0; i < list.length; i += 1) {
      if (list[i].value === value) return list[i];
    }
    return list[0];
  }

  function polar(cx, cy, radius, angle) {
    return [cx + radius * Math.cos(angle), cy + radius * Math.sin(angle)];
  }

  // One TILE_W x TILE_H pie drawn at the origin, without the <svg> wrapper, so a
  // tile can be either a standalone chart or one cell of the condition grid.
  function pieBody(scope) {
    var totals = domainTotals(scope.cols);
    var domains = Object.keys(totals).sort(function (a, b) { return totals[b] - totals[a]; });
    var grand = domains.reduce(function (s, d) { return s + totals[d]; }, 0);

    var cx = 160, cy = 178, R = 130;
    var parts = ['<rect width="' + TILE_W + '" height="' + TILE_H + '" fill="#ffffff"/>'];
    var title = (matrix.preset_label || "eToL") + " — non-human reads by domain";
    parts.push('<text x="16" y="24" font-size="14" font-weight="600">' + esc(title) + "</text>");
    parts.push('<text x="16" y="42" font-size="11" fill="#555">' + esc(scope.subtitle) + "</text>");

    if (grand <= 0) {
      parts.push('<text x="16" y="80" font-size="12" fill="#666">No non-human probe-matched reads to plot.</text>');
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
    parts.push('<text x="16" y="' + (TILE_H - 12) + '" font-size="10" fill="#777">Total non-human reads: ' + grand + "</text>");
    return parts.join("");
  }

  function svgWrap(width, height, body) {
    return '<svg xmlns="http://www.w3.org/2000/svg" width="' + width + '" height="' + height +
      '" viewBox="0 0 ' + width + " " + height + '" font-family="system-ui, sans-serif">' + body + "</svg>";
  }

  function buildSvg(scope) {
    return svgWrap(TILE_W, TILE_H, pieBody(scope));
  }

  // Contact sheet: one labelled pie per condition, two per row. Export PNG/SVG
  // serialise whatever is on the canvas, so this is also the summary image.
  function buildGridSvg() {
    var conditions = conditionNames();
    var ncols = conditions.length > 1 ? 2 : 1;
    var nrows = Math.ceil(conditions.length / ncols);
    var body = [];
    conditions.forEach(function (name, i) {
      var x = (i % ncols) * TILE_W;
      var y = Math.floor(i / ncols) * TILE_H;
      body.push('<g transform="translate(' + x + "," + y + ')">');
      body.push(pieBody(conditionScope(name)));
      body.push('<rect width="' + TILE_W + '" height="' + TILE_H + '" fill="none" stroke="#e5e5e5"/>');
      body.push("</g>");
    });
    return svgWrap(ncols * TILE_W, nrows * TILE_H, body.join(""));
  }

  function render() {
    var scope = currentScope();
    var svg = scope.value === GRID ? buildGridSvg() : buildSvg(scope);
    root.querySelector("[data-role=canvas]").innerHTML = svg;
  }

  function populateScopes() {
    var sel = root.querySelector("[data-role=sample]");
    if (!sel) return;
    var conditions = conditionNames();
    var html = '<option value="all">All samples (whole job)</option>';
    if (conditions.length) {
      html += '<option value="' + GRID + '">All conditions (grid)</option>';
      html += '<optgroup label="Conditions">';
      conditions.forEach(function (name) {
        html += '<option value="' + esc("cond:" + name) + '">' + esc(name) + "</option>";
      });
      html += "</optgroup>";
    }
    html += '<optgroup label="Samples">';
    matrix.cols.forEach(function (col, i) {
      html += '<option value="s:' + i + '">' + esc(col.sample) + "</option>";
    });
    html += "</optgroup>";
    sel.innerHTML = html;
  }

  function currentSvg() {
    var svg = root.querySelector("[data-role=canvas] svg");
    return svg ? svg.outerHTML : null;
  }

  // The grid is a different deliverable from a single pie -- name it so the two
  // don't collide in the download folder.
  function exportName(extension) {
    var stem = currentScope().value === GRID ? "etol_domain_pie_by_condition" : "etol_domain_pie";
    return stem + "." + extension;
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
    var name = exportName("png");
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
        a.download = name;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(a.href);
      });
    };
    img.src = url;
  }

  // One CSV of every pie the drop-down can draw: the whole job, each condition,
  // and each sample, as (scope, domain, reads, percent) rows. Reuses domainTotals
  // so the numbers match each chart's legend exactly (zero-total domains omitted,
  // largest first).
  function exportCsv() {
    if (!matrix) return;
    var lines = ["Scope,Domain,Reads,Percent"];
    scopes().forEach(function (scope) {
      if (scope.value === GRID) return; // a view of the condition scopes, not a scope
      var totals = domainTotals(scope.cols);
      var domains = Object.keys(totals).sort(function (a, b) { return totals[b] - totals[a]; });
      var grand = domains.reduce(function (s, d) { return s + totals[d]; }, 0);
      if (grand <= 0) {
        lines.push([scope.csvLabel, "", 0, "0.0"].map(csvCell).join(","));
        return;
      }
      domains.forEach(function (d) {
        var pct = (totals[d] / grand * 100).toFixed(1);
        lines.push([scope.csvLabel, d, totals[d], pct].map(csvCell).join(","));
      });
    });
    download("etol_domain_pie.csv", "text/csv", lines.join("\r\n"));
  }

  root.querySelector("[data-role=png]").addEventListener("click", exportPng);
  root.querySelector("[data-role=svg]").addEventListener("click", function () {
    var svgText = currentSvg();
    if (svgText) download(exportName("svg"), "image/svg+xml", svgText);
  });
  root.querySelector("[data-role=csv]").addEventListener("click", exportCsv);
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
      populateScopes();
      setStatus("");
      render();
    })
    .catch(function (err) {
      setStatus("Failed to load pie chart: " + err.message);
    });
})();
