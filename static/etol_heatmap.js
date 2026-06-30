/*
 * eToL / eToL-V result heatmap.
 *
 * Renders the rows-by-samples hit matrix served by
 * /batch-results/<id>/etol-matrix.json as an inline SVG, with no external
 * charting library so the standalone build stays lean and works offline.
 *
 * The same component serves both panels; it just picks paper-faithful defaults:
 *   - eToL-V (viral):   raw matched-read counts, blue->white->red (Veso, Fig 8/10)
 *   - eToL (cellular):  log2 reads per host cell with a 3 reads/cell cutoff
 *                       (Hu/Lathe pheatmap/Morpheus)
 * The user can switch the plotted value, the raw vs validated stage (when contig
 * identification ran), the per-row level, and the cutoff, then export PNG/SVG.
 */
(function () {
  "use strict";

  var root = document.getElementById("etol-heatmap");
  if (!root) {
    return;
  }

  var MATRIX_URL = root.dataset.matrixUrl;
  var IS_VIRAL = root.dataset.isViral === "1";

  var CELL_W = 16;
  var CELL_H = 12;
  var MARGIN = { top: 56, right: 84, bottom: 140, left: 178, group: 7 };

  // Diverging blue -> white -> red, matching the eToL-V dissertation figures.
  var COOL = [59, 76, 192];
  var MID = [221, 221, 221];
  var WARM = [180, 4, 38];

  // Stable categorical palette for the row-group (class / viral family) strip.
  var GROUP_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f",
  ];

  var state = {
    level: "species",
    value: IS_VIRAL ? "hits" : "log2",
    stage: "raw",
    cutoff: IS_VIRAL ? 0 : 3,
  };
  var cache = {}; // level -> matrix payload
  var groupColors = {};

  function lerp(a, b, t) {
    return Math.round(a + (b - a) * t);
  }

  function coolwarm(t) {
    t = Math.max(0, Math.min(1, t));
    var from = COOL;
    var to = MID;
    var u = t / 0.5;
    if (t >= 0.5) {
      from = MID;
      to = WARM;
      u = (t - 0.5) / 0.5;
    }
    return (
      "rgb(" +
      lerp(from[0], to[0], u) +
      "," +
      lerp(from[1], to[1], u) +
      "," +
      lerp(from[2], to[2], u) +
      ")"
    );
  }

  // Fixed colours for the disease tags both eToL papers use, so the familiar
  // AD/control palette is preserved. Any other label (e.g. from an uploaded
  // design matrix) gets a stable colour assigned from CONDITION_PALETTE, the
  // same first-seen scheme groupColor() uses for the row-group strip.
  var CONDITION_DEFAULTS = {
    AD: "#c0392b", CTRL: "#2e7d32", CONTROL: "#2e7d32",
    LBD: "#6a1b9a", VaD: "#e65100",
  };
  var CONDITION_PALETTE = [
    "#00897b", "#3949ab", "#c2185b", "#f9a825", "#5d4037",
    "#546e7a", "#ad1457", "#00838f", "#6d4c41", "#283593",
  ];
  var conditionColors = {};

  function conditionColor(cond) {
    if (!cond) return "#cfd8dc";
    if (cond in CONDITION_DEFAULTS) return CONDITION_DEFAULTS[cond];
    // Prefix match so combined tags like "AD/LBD" still read as AD red, matching
    // the original sample-name behaviour.
    var keys = Object.keys(CONDITION_DEFAULTS);
    for (var i = 0; i < keys.length; i += 1) {
      if (cond.indexOf(keys[i]) === 0) return CONDITION_DEFAULTS[keys[i]];
    }
    if (!(cond in conditionColors)) {
      conditionColors[cond] =
        CONDITION_PALETTE[Object.keys(conditionColors).length % CONDITION_PALETTE.length];
    }
    return conditionColors[cond];
  }

  function groupColor(group) {
    if (!(group in groupColors)) {
      groupColors[group] = GROUP_PALETTE[Object.keys(groupColors).length % GROUP_PALETTE.length];
    }
    return groupColors[group];
  }

  function esc(text) {
    return String(text)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  }

  function fmt(value) {
    if (value === null) return "n/a";
    if (state.value === "hits") return String(value);
    if (Math.abs(value) >= 100) return value.toFixed(0);
    return parseFloat(value.toPrecision(3)).toString();
  }

  // Per-cell transformed value used for both colour and the hover tooltip.
  function cellValue(matrix, r, c) {
    var counts = state.stage === "validated" && matrix.confirmed ? matrix.confirmed : matrix.hits;
    var raw = counts[r][c];
    if (raw === null) {
      return { display: null, na: true, passes: false };
    }
    var hostCells = matrix.cols[c].host_cells;
    if (state.value === "hits") {
      return { display: raw, na: false, passes: raw >= state.cutoff && raw > 0 };
    }
    if (hostCells <= 0) {
      return { display: null, na: true, passes: false };
    }
    var rpc = raw / hostCells;
    var passes = rpc >= state.cutoff && rpc > 0;
    if (state.value === "log2") {
      return { display: passes ? Math.log2(rpc) : null, rpc: rpc, na: false, passes: passes };
    }
    return { display: rpc, na: false, passes: passes };
  }

  function buildSvg(matrix) {
    var rows = matrix.rows;
    var cols = matrix.cols;
    var nRows = rows.length;
    var nCols = cols.length;
    var showCellCount = cols.some(function (col) { return col.host_cells > 0; });

    // First pass: transformed values + colour domain over the passing cells.
    var values = [];
    var domainMin = state.value === "log2" ? Infinity : 0;
    var domainMax = -Infinity;
    for (var r = 0; r < nRows; r += 1) {
      values.push([]);
      for (var c = 0; c < nCols; c += 1) {
        var cell = cellValue(matrix, r, c);
        values[r].push(cell);
        if (cell.passes && cell.display !== null) {
          if (cell.display > domainMax) domainMax = cell.display;
          if (state.value === "log2" && cell.display < domainMin) domainMin = cell.display;
        }
      }
    }
    if (domainMax === -Infinity) domainMax = 1;
    if (domainMin === Infinity) domainMin = 0;
    if (domainMin > 0) domainMin = 0;
    var span = domainMax - domainMin || 1;

    function fillFor(cell) {
      if (cell.na) return "#eceff1";
      if (!cell.passes) return coolwarm(0); // below-cutoff / zero -> background
      return coolwarm((cell.display - domainMin) / span);
    }

    var gridW = nCols * CELL_W;
    var gridH = nRows * CELL_H;
    var extra = showCellCount ? CELL_H + 26 : 0;
    var width = MARGIN.left + gridW + MARGIN.right;
    var height = MARGIN.top + gridH + extra + MARGIN.bottom;

    var parts = [];
    parts.push(
      '<svg xmlns="http://www.w3.org/2000/svg" width="' + width + '" height="' + height +
        '" viewBox="0 0 ' + width + " " + height + '" font-family="system-ui, sans-serif">'
    );
    parts.push('<rect width="' + width + '" height="' + height + '" fill="#ffffff"/>');

    var title = (matrix.preset_label || (IS_VIRAL ? "eToL-V" : "eToL")) +
      " probe hits  |  " + valueLabel() + (state.stage === "validated" ? "  (validated)" : "");
    parts.push('<text x="' + MARGIN.left + '" y="22" font-size="14" font-weight="600">' + esc(title) + "</text>");

    // Condition annotation strip above the columns.
    var stripY = MARGIN.top - 14;
    for (var cc = 0; cc < nCols; cc += 1) {
      var x = MARGIN.left + cc * CELL_W;
      parts.push(
        '<rect x="' + x + '" y="' + stripY + '" width="' + CELL_W + '" height="10" fill="' +
          conditionColor(cols[cc].condition) + '"><title>' + esc(cols[cc].sample) +
          (cols[cc].condition ? " (" + esc(cols[cc].condition) + ")" : "") + "</title></rect>"
      );
    }

    // Cells.
    for (var ri = 0; ri < nRows; ri += 1) {
      var y = MARGIN.top + ri * CELL_H;
      parts.push(
        '<rect x="' + (MARGIN.left - MARGIN.group - 1) + '" y="' + y + '" width="' + MARGIN.group +
          '" height="' + CELL_H + '" fill="' + groupColor(rows[ri].group) + '"/>'
      );
      parts.push(
        '<text x="' + (MARGIN.left - MARGIN.group - 5) + '" y="' + (y + CELL_H - 2) +
          '" font-size="9" text-anchor="end" fill="#333">' + esc(rows[ri].label) + "</text>"
      );
      for (var ci = 0; ci < nCols; ci += 1) {
        var cellX = MARGIN.left + ci * CELL_W;
        var cd = values[ri][ci];
        parts.push(
          '<rect x="' + cellX + '" y="' + y + '" width="' + CELL_W + '" height="' + CELL_H +
            '" fill="' + fillFor(cd) + '" stroke="#ffffff" stroke-width="0.5"><title>' +
            esc(rows[ri].label + " × " + cols[ci].sample + " = " + fmt(cd.display)) +
            "</title></rect>"
        );
      }
    }

    // Optional "Cell count" row (host-cell estimate), as in the eToL-V figures.
    var afterGrid = MARGIN.top + gridH;
    if (showCellCount) {
      var ccY = afterGrid + 6;
      var maxCells = Math.max.apply(null, cols.map(function (col) { return col.host_cells; })) || 1;
      parts.push(
        '<text x="' + (MARGIN.left - MARGIN.group - 5) + '" y="' + (ccY + CELL_H - 2) +
          '" font-size="9" text-anchor="end" fill="#333">Cell count</text>'
      );
      for (var k = 0; k < nCols; k += 1) {
        var hc = cols[k].host_cells;
        var shade = 235 - Math.round((hc / maxCells) * 150);
        parts.push(
          '<rect x="' + (MARGIN.left + k * CELL_W) + '" y="' + ccY + '" width="' + CELL_W +
            '" height="' + CELL_H + '" fill="rgb(' + shade + "," + shade + "," + shade +
            ')" stroke="#fff" stroke-width="0.5"><title>' + esc(cols[k].sample + ": " + hc + " est. host cells") +
            "</title></rect>"
        );
      }
      afterGrid = ccY + CELL_H;
    }

    // Column labels, rotated under the grid.
    var labelY = afterGrid + 6;
    for (var cl = 0; cl < nCols; cl += 1) {
      var lx = MARGIN.left + cl * CELL_W + CELL_W / 2;
      parts.push(
        '<text x="' + lx + '" y="' + labelY + '" font-size="9" fill="#333" text-anchor="end" transform="rotate(-60 ' +
          lx + " " + labelY + ')">' + esc(cols[cl].sample) + "</text>"
      );
    }

    parts.push(legendSvg(width, domainMin, domainMax));
    parts.push("</svg>");
    return parts.join("");
  }

  function legendSvg(width, domainMin, domainMax) {
    var x = width - MARGIN.right + 18;
    var top = MARGIN.top;
    var h = 140;
    var out = [
      '<defs><linearGradient id="etolGrad" x1="0" y1="1" x2="0" y2="0">' +
        '<stop offset="0%" stop-color="' + coolwarm(0) + '"/>' +
        '<stop offset="50%" stop-color="' + coolwarm(0.5) + '"/>' +
        '<stop offset="100%" stop-color="' + coolwarm(1) + '"/></linearGradient></defs>',
      '<rect x="' + x + '" y="' + top + '" width="12" height="' + h + '" fill="url(#etolGrad)" stroke="#999" stroke-width="0.5"/>',
      '<text x="' + (x + 16) + '" y="' + (top + 8) + '" font-size="9">' + esc(fmt(domainMax)) + "</text>",
      '<text x="' + (x + 16) + '" y="' + (top + h) + '" font-size="9">' + esc(fmt(domainMin)) + "</text>",
      '<text x="' + (x + 16) + '" y="' + (top + h + 16) + '" font-size="9" fill="#555">' + esc(valueLabel()) + "</text>",
    ];
    return out.join("");
  }

  function valueLabel() {
    if (state.value === "hits") return "matched reads";
    if (state.value === "reads_per_cell") return "reads / host cell";
    return "log2(reads / host cell)";
  }

  function setStatus(message) {
    root.querySelector("[data-role=status]").textContent = message || "";
  }

  // Rebuild the condition legend from the labels actually present, so an uploaded
  // design matrix with arbitrary labels gets matching swatches (the static
  // AD/Control/LBD/VaD list is just the default before any matrix is applied).
  function renderConditionLegend(matrix) {
    var el = root.querySelector("[data-role=condition-legend]");
    if (!el) return;
    var seen = {};
    var labels = [];
    matrix.cols.forEach(function (col) {
      var c = col.condition;
      if (c && !(c in seen)) {
        seen[c] = true;
        labels.push(c);
      }
    });
    if (!labels.length) {
      el.style.display = "none";
      return;
    }
    el.style.display = "";
    var parts = ["Condition: "];
    labels.forEach(function (label) {
      parts.push(
        '<span class="swatch" style="background:' + conditionColor(label) + '"></span>' +
          esc(label) + " "
      );
    });
    el.innerHTML = parts.join("");
  }

  // Note which design matrix (if any) drove the labels, and flag samples it did
  // not cover so the gap is visible rather than silently mislabeled.
  function renderDesignNote(matrix) {
    var el = root.querySelector("[data-role=design-note]");
    if (!el) return;
    var design = matrix.design_matrix;
    if (!design) {
      el.textContent = "";
      el.style.display = "none";
      return;
    }
    el.style.display = "";
    var msg =
      "Condition labels from design matrix “" + design.source + "” (" +
      design.row_count + " sample" + (design.row_count === 1 ? "" : "s") + ").";
    var unmatched = design.unmatched_samples || [];
    if (unmatched.length) {
      msg +=
        " " + unmatched.length + " selected sample" +
        (unmatched.length === 1 ? "" : "s") + " had no matrix row and render unlabeled: " +
        unmatched.join(", ") + ".";
    }
    el.textContent = msg;
  }

  function render() {
    var matrix = cache[state.level];
    if (!matrix) return;
    var hasValidated = !!matrix.confirmed;
    var stageSelect = root.querySelector("[data-role=stage]");
    stageSelect.disabled = !hasValidated;
    if (!hasValidated && state.stage === "validated") {
      state.stage = "raw";
      stageSelect.value = "raw";
    }
    if (!matrix.rows.length || !matrix.cols.length) {
      root.querySelector("[data-role=canvas]").innerHTML =
        '<p class="muted">No matrix to plot for this batch.</p>';
      setStatus("");
      return;
    }
    root.querySelector("[data-role=canvas]").innerHTML = buildSvg(matrix);
    renderConditionLegend(matrix);
    renderDesignNote(matrix);
    setStatus(
      matrix.rows.length + " rows × " + matrix.cols.length + " samples" +
        (hasValidated ? "" : "  —  no validated layer (run contig identification for the raw→validated toggle)")
    );
  }

  function load() {
    if (cache[state.level]) {
      render();
      return;
    }
    setStatus("Loading…");
    fetch(MATRIX_URL + "?level=" + state.level)
      .then(function (resp) {
        if (!resp.ok) throw new Error("HTTP " + resp.status);
        return resp.json();
      })
      .then(function (matrix) {
        cache[state.level] = matrix;
        render();
      })
      .catch(function (err) {
        setStatus("Failed to load heatmap: " + err.message);
      });
  }

  function download(name, mime, data) {
    var blob = new Blob([data], { type: mime });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a");
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  function currentSvg() {
    var svg = root.querySelector("[data-role=canvas] svg");
    return svg ? svg.outerHTML : null;
  }

  function exportPng() {
    var svgText = currentSvg();
    if (!svgText) return;
    var svgEl = root.querySelector("[data-role=canvas] svg");
    var w = parseInt(svgEl.getAttribute("width"), 10);
    var h = parseInt(svgEl.getAttribute("height"), 10);
    var img = new Image();
    var url = URL.createObjectURL(new Blob([svgText], { type: "image/svg+xml" }));
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
        a.download = "etol_heatmap.png";
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(a.href);
      });
    };
    img.src = url;
  }

  function wire(rolename, handler) {
    var el = root.querySelector("[data-role=" + rolename + "]");
    if (el) el.addEventListener("change", handler);
  }

  wire("level", function (e) {
    state.level = e.target.value;
    load();
  });
  wire("value", function (e) {
    state.value = e.target.value;
    state.cutoff = state.value === "hits" ? 0 : IS_VIRAL ? 0.03 : 3;
    root.querySelector("[data-role=cutoff]").value = state.cutoff;
    render();
  });
  wire("stage", function (e) {
    state.stage = e.target.value;
    render();
  });
  wire("cutoff", function (e) {
    var parsed = parseFloat(e.target.value);
    state.cutoff = isNaN(parsed) ? 0 : parsed;
    render();
  });
  root.querySelector("[data-role=png]").addEventListener("click", exportPng);
  root.querySelector("[data-role=svg]").addEventListener("click", function () {
    var svgText = currentSvg();
    if (svgText) download("etol_heatmap.svg", "image/svg+xml", svgText);
  });

  // Reflect the per-preset defaults into the controls before the first render.
  root.querySelector("[data-role=value]").value = state.value;
  root.querySelector("[data-role=cutoff]").value = state.cutoff;
  load();
})();
