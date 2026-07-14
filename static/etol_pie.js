/*
 * eToL / eToL-V domain composition chart (pie + stacked bars).
 *
 * A coarser companion to the heatmap: what share of a job's non-human
 * probe-matched reads falls into each Tree-of-Life domain (Archaea, Bacteria,
 * Fungi, ...). Reuses the same /batch-results/<id>/etol-matrix.json the heatmap
 * loads -- no extra endpoint -- and draws inline SVG so the standalone build
 * stays lean and works offline.
 *
 * Two marks over one aggregation:
 *   - Pie:          proportions within a single scope (one sample/condition/job).
 *   - Stacked bars: one bar per sample (or per pooled condition), so the
 *                   magnitude a pie throws away is comparable ACROSS samples --
 *                   bar height is total abundance, segments are the pie.
 *
 * Bars default to reads per host cell, not raw reads: library depth and host
 * content differ per sample, so raw counts are not comparable across columns --
 * that is what the PGK1/hNSE control probes exist to correct (Hu, Haas & Lathe
 * 2022), and it keeps this chart on the same axis as the heatmap below it. Raw
 * counts and 100%-stacked (every pie side by side) stay available as options.
 *
 * Scopes: the whole job, one design-matrix condition (every sample in that group
 * pooled), one sample, or all conditions at once -- a grid of pies, or one pooled
 * bar per condition. The condition labels are matrix.cols[i].condition, already
 * served for the heatmap's condition strip, so nothing new is fetched.
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
  var GRID = "grid"; // scope value for "all conditions" (pie grid / pooled bars)
  var TILE_W = 540, TILE_H = 340;
  var BAR = { w: 46, gap: 16, left: 78, right: 200, top: 58, bottom: 104, h: 300 };
  var VALUE_LABELS = {
    per_cell: "reads per host cell",
    hits: "matched reads (raw)",
    pct: "% of non-human reads",
  };

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

  // --- Stacked bars ---------------------------------------------------------

  function chartMode() {
    var sel = root.querySelector("[data-role=chart]");
    return sel ? sel.value : "pie";
  }

  function barValue() {
    var sel = root.querySelector("[data-role=value]");
    return sel ? sel.value : "per_cell";
  }

  // One bar per sample; the "all conditions" scope instead pools each condition
  // into a single bar.
  function barSeries(scope) {
    if (scope.value === GRID) {
      return conditionNames().map(function (name) {
        return { label: name, condition: name, cols: conditionScope(name).cols };
      });
    }
    return scope.cols.map(function (c) {
      return {
        label: matrix.cols[c].sample,
        condition: String(matrix.cols[c].condition || "").trim(),
        cols: [c],
      };
    });
  }

  // Stacked segment values for one bar, in the current y-axis units. Pooling a
  // condition sums reads AND host cells (depth-weighted) rather than averaging
  // per-sample ratios, so a deeper library carries proportional weight -- the
  // same aggregation the pooled pie already does.
  // ponytail: depth-weighted pooling; if per-sample spread matters more than the
  // group total, plot samples (scope = All samples) rather than adding error bars.
  function barTotals(series) {
    var raw = domainTotals(series.cols);
    var domains = Object.keys(raw);
    var grand = domains.reduce(function (s, d) { return s + raw[d]; }, 0);
    var mode = barValue();

    if (mode === "hits") {
      return { values: raw, total: grand, na: false };
    }
    if (mode === "pct") {
      var pct = {};
      domains.forEach(function (d) { pct[d] = grand > 0 ? (raw[d] / grand) * 100 : 0; });
      return { values: pct, total: grand > 0 ? 100 : 0, na: false };
    }
    // reads per host cell: undefined when the control probes recovered no reads.
    var hostCells = series.cols.reduce(function (s, c) {
      return s + (matrix.cols[c].host_cells || 0);
    }, 0);
    if (hostCells <= 0) {
      return { values: {}, total: 0, na: true };
    }
    var perCell = {};
    domains.forEach(function (d) { perCell[d] = raw[d] / hostCells; });
    return { values: perCell, total: grand / hostCells, na: false };
  }

  // ~5 axis ticks on the 1/2/5 x 10^n ladder.
  function ticks(max) {
    if (!(max > 0)) return [0, 1];
    var step = Math.pow(10, Math.floor(Math.log10(max / 5)));
    var err = max / 5 / step;
    if (err >= 5) step *= 10;
    else if (err >= 2) step *= 5;
    else if (err >= 1) step *= 2;
    var out = [];
    for (var v = 0; v <= max + step / 2; v += step) out.push(v);
    return out;
  }

  function num(value) {
    return String(parseFloat(value.toPrecision(3)));
  }

  function buildBarsSvg(scope) {
    var series = barSeries(scope);
    var data = series.map(barTotals);

    // One stack order and one legend for every bar: domains ranked by their total
    // over the whole chart, largest at the bottom.
    var everyCol = series.reduce(function (acc, s) { return acc.concat(s.cols); }, []);
    var overall = domainTotals(everyCol);
    var domains = Object.keys(overall).sort(function (a, b) { return overall[b] - overall[a]; });

    var maxTotal = data.reduce(function (m, d) { return Math.max(m, d.total); }, 0);
    var tickVals = ticks(maxTotal);
    var top = tickVals[tickVals.length - 1] || 1;

    var plotW = series.length * (BAR.w + BAR.gap) + BAR.gap;
    var width = BAR.left + plotW + BAR.right;
    var height = BAR.top + BAR.h + BAR.bottom;
    var baseY = BAR.top + BAR.h;

    var parts = ['<rect width="' + width + '" height="' + height + '" fill="#ffffff"/>'];
    var title = (matrix.preset_label || "eToL") + " — non-human reads by domain";
    parts.push('<text x="16" y="24" font-size="14" font-weight="600">' + esc(title) + "</text>");
    parts.push('<text x="16" y="42" font-size="11" fill="#555">' +
      esc(scope.subtitle || scope.label) + " — " + esc(VALUE_LABELS[barValue()]) + "</text>");

    if (!domains.length) {
      parts.push('<text x="16" y="80" font-size="12" fill="#666">No non-human probe-matched reads to plot.</text>');
      return svgWrap(width, height, parts.join(""));
    }

    // Gridlines + y axis.
    tickVals.forEach(function (v) {
      var y = baseY - (v / top) * BAR.h;
      parts.push('<line x1="' + BAR.left + '" y1="' + y.toFixed(1) + '" x2="' + (BAR.left + plotW) +
        '" y2="' + y.toFixed(1) + '" stroke="#e9e9e9" stroke-width="1"/>');
      parts.push('<text x="' + (BAR.left - 8) + '" y="' + (y + 3.5).toFixed(1) +
        '" font-size="10" fill="#555" text-anchor="end">' + esc(num(v)) + "</text>");
    });
    parts.push('<line x1="' + BAR.left + '" y1="' + BAR.top + '" x2="' + BAR.left + '" y2="' + baseY +
      '" stroke="#999" stroke-width="1"/>');
    parts.push('<line x1="' + BAR.left + '" y1="' + baseY + '" x2="' + (BAR.left + plotW) + '" y2="' + baseY +
      '" stroke="#999" stroke-width="1"/>');
    parts.push('<text transform="rotate(-90 18 ' + (BAR.top + BAR.h / 2) + ')" x="18" y="' +
      (BAR.top + BAR.h / 2) + '" font-size="11" fill="#333" text-anchor="middle">' +
      esc(VALUE_LABELS[barValue()]) + "</text>");

    // Bars: stacked bottom-up in the shared domain order.
    series.forEach(function (s, i) {
      var x = BAR.left + BAR.gap + i * (BAR.w + BAR.gap);
      var cursor = baseY;
      var cell = data[i];

      if (cell.na) {
        parts.push('<text x="' + (x + BAR.w / 2) + '" y="' + (baseY - 8) +
          '" font-size="10" fill="#999" text-anchor="middle">n/a</text>');
      }
      domains.forEach(function (d) {
        var value = cell.values[d] || 0;
        if (value <= 0) return;
        var h = (value / top) * BAR.h;
        cursor -= h;
        var tip = esc(s.label + " — " + d + ": " + num(value) + " " + VALUE_LABELS[barValue()]);
        parts.push('<rect x="' + x + '" y="' + cursor.toFixed(2) + '" width="' + BAR.w + '" height="' +
          h.toFixed(2) + '" fill="' + colorFor(d) + '" stroke="#ffffff" stroke-width="0.5"><title>' +
          tip + "</title></rect>");
      });

      // Sample label (rotated) + its condition underneath.
      var lx = x + BAR.w / 2;
      var ly = baseY + 10;
      parts.push('<text x="' + lx + '" y="' + ly + '" font-size="10" fill="#333" text-anchor="end" ' +
        'transform="rotate(-45 ' + lx + " " + ly + ')">' + esc(s.label) + "</text>");
      if (s.condition && scope.value !== GRID) {
        parts.push('<text x="' + lx + '" y="' + (height - 14) + '" font-size="9" fill="#777" ' +
          'text-anchor="middle">' + esc(s.condition) + "</text>");
      }
    });

    // Legend on the right, same order as the stack.
    var lgx = BAR.left + plotW + 20;
    domains.forEach(function (d, i) {
      var y = BAR.top + 10 + i * 22;
      parts.push('<rect x="' + lgx + '" y="' + (y - 11) + '" width="13" height="13" fill="' + colorFor(d) + '"/>');
      parts.push('<text x="' + (lgx + 19) + '" y="' + y + '" font-size="11" fill="#333">' + esc(d) + "</text>");
    });

    var hint = barValue() === "hits"
      ? "Raw counts are not depth-normalized — compare across samples with care."
      : "";
    if (hint) {
      parts.push('<text x="16" y="' + (height - 4) + '" font-size="10" fill="#a15c00">' + esc(hint) + "</text>");
    }
    return svgWrap(width, height, parts.join(""));
  }

  function render() {
    var scope = currentScope();
    var svg;
    if (chartMode() === "bars") {
      svg = buildBarsSvg(scope);
    } else if (scope.value === GRID) {
      svg = buildGridSvg();
    } else {
      svg = buildSvg(scope);
    }
    root.querySelector("[data-role=canvas]").innerHTML = svg;

    // The value axis is a bars-only choice; the pie is always a proportion.
    var valueSel = root.querySelector("[data-role=value]");
    if (valueSel) valueSel.disabled = chartMode() !== "bars";
  }

  function populateScopes() {
    var sel = root.querySelector("[data-role=sample]");
    if (!sel) return;
    var previous = sel.value;
    var conditions = conditionNames();
    var gridLabel = chartMode() === "bars"
      ? "All conditions (one pooled bar each)"
      : "All conditions (grid)";
    var html = '<option value="all">All samples (whole job)</option>';
    if (conditions.length) {
      html += '<option value="' + GRID + '">' + esc(gridLabel) + "</option>";
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
    if (previous) sel.value = previous; // keep the scope when only the mark changed
  }

  function currentSvg() {
    var svg = root.querySelector("[data-role=canvas] svg");
    return svg ? svg.outerHTML : null;
  }

  // Each view is a different deliverable -- name them so they don't collide in
  // the download folder.
  function exportName(extension) {
    var stem = "etol_domain_pie";
    if (chartMode() === "bars") {
      stem = "etol_domain_bars";
    } else if (currentScope().value === GRID) {
      stem = "etol_domain_pie_by_condition";
    }
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

  // One CSV of every scope the drop-down can draw: the whole job, each condition,
  // and each sample, as (scope, domain, reads, percent, reads per host cell) rows
  // -- i.e. the numbers behind both marks, on all three bar axes. Reuses
  // domainTotals so they match each chart exactly (zero-total domains omitted,
  // largest first). Reads per host cell is blank when the control probes recovered
  // no reads, so an unnormalizable sample reads as such instead of as zero.
  function exportCsv() {
    if (!matrix) return;
    var lines = ["Scope,Domain,Reads,Percent,Reads per host cell"];
    scopes().forEach(function (scope) {
      if (scope.value === GRID) return; // a view of the condition scopes, not a scope
      var totals = domainTotals(scope.cols);
      var domains = Object.keys(totals).sort(function (a, b) { return totals[b] - totals[a]; });
      var grand = domains.reduce(function (s, d) { return s + totals[d]; }, 0);
      if (grand <= 0) {
        lines.push([scope.csvLabel, "", 0, "0.0", ""].map(csvCell).join(","));
        return;
      }
      var hostCells = scope.cols.reduce(function (s, c) {
        return s + (matrix.cols[c].host_cells || 0);
      }, 0);
      domains.forEach(function (d) {
        var pct = (totals[d] / grand * 100).toFixed(1);
        var perCell = hostCells > 0 ? num(totals[d] / hostCells) : "";
        lines.push([scope.csvLabel, d, totals[d], pct, perCell].map(csvCell).join(","));
      });
    });
    download("etol_domain_composition.csv", "text/csv", lines.join("\r\n"));
  }

  root.querySelector("[data-role=png]").addEventListener("click", exportPng);
  root.querySelector("[data-role=svg]").addEventListener("click", function () {
    var svgText = currentSvg();
    if (svgText) download(exportName("svg"), "image/svg+xml", svgText);
  });
  root.querySelector("[data-role=csv]").addEventListener("click", exportCsv);
  var sampleSel = root.querySelector("[data-role=sample]");
  if (sampleSel) sampleSel.addEventListener("change", render);
  var chartSel = root.querySelector("[data-role=chart]");
  if (chartSel) {
    chartSel.addEventListener("change", function () {
      populateScopes(); // the "all conditions" option means a grid or pooled bars
      render();
    });
  }
  var valueSel = root.querySelector("[data-role=value]");
  if (valueSel) valueSel.addEventListener("change", render);

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
