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

  // --- Nonparametric two-group test (pure; see tests/test_etol_pie_stats.js) --

  // Midranks of pooled values; tied values share their mean rank.
  function midranks(values) {
    var order = values.map(function (_, i) { return i; })
      .sort(function (i, j) { return values[i] - values[j]; });
    var ranks = new Array(values.length);
    for (var i = 0; i < order.length; ) {
      var j = i;
      while (j + 1 < order.length && values[order[j + 1]] === values[order[i]]) j += 1;
      var mean = (i + j) / 2 + 1;
      for (var k = i; k <= j; k += 1) ranks[order[k]] = mean;
      i = j + 1;
    }
    return ranks;
  }

  function choose(n, k) {
    var c = 1;
    for (var i = 0; i < k; i += 1) c = c * (n - i) / (i + 1);
    return Math.round(c);
  }

  // Fixed-seed LCG: a p-value must not change between two renders of the same
  // figure, so the sampled branch has to be deterministic.
  function lcg(seed) {
    var s = seed >>> 0;
    return function () { s = (Math.imul(s, 1103515245) + 12345) >>> 0; return s / 4294967296; };
  }

  var PERMUTATIONS = 20000;

  // Two-group label spaces this small are enumerated exactly instead of sampled.
  function enumerable(sizes) {
    return sizes.length === 2 && choose(sizes[0] + sizes[1], sizes[0]) <= 100000;
  }

  // Rank statistic S = sum over groups of R_g^2 / n_g.
  //
  // One statistic covers both designs. With two groups S is a parabola in the
  // first group's rank sum whose vertex sits exactly at the null mean, so
  // P(S >= S_obs) IS the two-sided Wilcoxon rank-sum p -- no separate two-sided
  // handling. With more groups S is Kruskal-Wallis H up to constants that are
  // fixed under permutation (the tie correction among them), so the permutation
  // p is Kruskal-Wallis's. Hence one test function, not two.
  function rankStat(ranks, sizes) {
    var s = 0, at = 0;
    for (var g = 0; g < sizes.length; g += 1) {
      var sum = 0;
      for (var i = 0; i < sizes[g]; i += 1) { sum += ranks[at]; at += 1; }
      s += sum * sum / sizes[g];
    }
    return s;
  }

  // p for "are these groups drawn from the same distribution?", by permuting the
  // group labels. Exact by enumeration for a small two-group split; otherwise
  // sampled -- a 4x6 design has ~2.3e12 label assignments, which neither needs
  // nor admits enumeration. Permuting is exact under ties, and eToL group sizes
  // are far too small for the normal/chi-square approximations the textbook
  // formulae rely on.
  function rankTestP(groups) {
    var sizes = groups.map(function (g) { return g.length; });
    if (!sizes.length || sizes.some(function (n) { return n < 1; })) return null;
    var pooled = [];
    groups.forEach(function (g) { pooled = pooled.concat(g); });
    var ranks = midranks(pooled);
    var n = ranks.length;
    var obs = rankStat(ranks, sizes) - 1e-9;
    var hits = 0, trials = 0, j, m, sum;

    if (enumerable(sizes)) {
      var total = n * (n + 1) / 2;
      var idx = [];
      for (j = 0; j < sizes[0]; j += 1) idx.push(j);
      for (;;) {
        sum = 0;
        for (j = 0; j < sizes[0]; j += 1) sum += ranks[idx[j]];
        if (sum * sum / sizes[0] + (total - sum) * (total - sum) / sizes[1] >= obs) hits += 1;
        trials += 1;
        j = sizes[0] - 1; // next combination, lexicographic
        while (j >= 0 && idx[j] === n - sizes[0] + j) j -= 1;
        if (j < 0) break;
        idx[j] += 1;
        for (m = j + 1; m < sizes[0]; m += 1) idx[m] = idx[m - 1] + 1;
      }
      return hits / trials;
    }

    var rand = lcg(20220317); // fixed seed: the same figure must give the same p
    var pool = ranks.slice();
    for (trials = 0; trials < PERMUTATIONS; trials += 1) {
      for (j = n - 1; j > 0; j -= 1) { // Fisher-Yates
        var t = Math.floor(rand() * (j + 1));
        var tmp = pool[j]; pool[j] = pool[t]; pool[t] = tmp;
      }
      if (rankStat(pool, sizes) >= obs) hits += 1;
    }
    return (hits + 1) / (trials + 1);
  }

  // Which test will run and the smallest p it can return, for the caption: the
  // exact floor when enumerated, otherwise the larger of the sampling resolution
  // and the true combinatorial floor (k! of the label assignments reproduce any
  // given grouping).
  function rankTestPlan(sizes) {
    var n = sizes.reduce(function (s, v) { return s + v; }, 0);
    if (enumerable(sizes)) {
      return { name: "Exact Wilcoxon rank-sum test", floor: 2 / choose(n, sizes[0]) };
    }
    var space = 1, left = n, fact = 1, i;
    sizes.forEach(function (size) { space *= choose(left, size); left -= size; });
    for (i = 2; i <= sizes.length; i += 1) fact *= i;
    return {
      name: sizes.length === 2
        ? "Permutation Wilcoxon rank-sum test"
        : "Kruskal-Wallis permutation test",
      floor: Math.max(1 / (PERMUTATIONS + 1), fact / space),
    };
  }

  // PERMANOVA (Anderson 2001) on Bray-Curtis distances between per-sample
  // composition vectors: the single global test of "does composition differ
  // between the groups at all", which the per-domain tests cannot answer jointly
  // -- they ask one question per domain, this asks one question of the whole
  // community. Pseudo-F is built from the distance matrix directly (no
  // ordination), p by permuting group labels, R^2 = share of the total
  // dissimilarity explained by the grouping.
  // ponytail: no betadisper companion -- PERMANOVA cannot tell a shift in group
  // centroid from a difference in within-group spread, so a small p with visibly
  // uneven scatter needs that check before it is called a composition shift.
  function brayCurtis(x, y) {
    var num = 0, den = 0;
    for (var i = 0; i < x.length; i += 1) {
      num += Math.abs(x[i] - y[i]);
      den += x[i] + y[i];
    }
    return den > 0 ? num / den : 0;
  }

  function permanova(vectors, labels) {
    var n = vectors.length;
    var names = labels.filter(function (v, i) { return labels.indexOf(v) === i; });
    var a = names.length;
    if (n < 3 || a < 2 || n - a < 1) return null;

    var d2 = [], i, j;
    for (i = 0; i < n; i += 1) {
      d2.push([]);
      for (j = 0; j < n; j += 1) {
        d2[i][j] = i === j ? 0 : Math.pow(brayCurtis(vectors[i], vectors[j]), 2);
      }
    }

    function stat(lab) {
      var total = 0, within = 0, counts = {}, sums = {}, x, y;
      for (x = 0; x < n; x += 1) {
        counts[lab[x]] = (counts[lab[x]] || 0) + 1;
        for (y = x + 1; y < n; y += 1) {
          total += d2[x][y];
          if (lab[x] === lab[y]) sums[lab[x]] = (sums[lab[x]] || 0) + d2[x][y];
        }
      }
      Object.keys(sums).forEach(function (g) { within += sums[g] / counts[g]; });
      total /= n;
      var between = total - within;
      return { f: (between / (a - 1)) / (within / (n - a)), r2: total > 0 ? between / total : 0 };
    }

    var obs = stat(labels);
    if (!isFinite(obs.f)) return null;
    var rand = lcg(19710407);
    var perm = labels.slice(), hits = 0, trials = 9999, t, k, s, tmp;
    for (t = 0; t < trials; t += 1) {
      for (k = n - 1; k > 0; k -= 1) {
        s = Math.floor(rand() * (k + 1));
        tmp = perm[k]; perm[k] = perm[s]; perm[s] = tmp;
      }
      if (stat(perm).f >= obs.f - 1e-12) hits += 1;
    }
    return { f: obs.f, r2: obs.r2, p: (hits + 1) / (trials + 1), permutations: trials };
  }

  // Benjamini-Hochberg adjusted p-values, returned in the input order.
  function bh(ps) {
    var m = ps.length;
    var order = ps.map(function (_, i) { return i; })
      .sort(function (i, j) { return ps[j] - ps[i]; });
    var out = new Array(m);
    var running = 1;
    order.forEach(function (i, k) {
      running = Math.min(running, ps[i] * m / (m - k));
      out[i] = running;
    });
    return out;
  }

  function stars(q) {
    if (q < 0.001) return "***";
    if (q < 0.01) return "**";
    if (q < 0.05) return "*";
    return "";
  }

  function median(values) {
    var s = values.slice().sort(function (x, y) { return x - y; });
    var mid = s.length >> 1;
    return s.length % 2 ? s[mid] : (s[mid - 1] + s[mid]) / 2;
  }

  if (typeof module !== "undefined" && module.exports) {
    module.exports = {
      midranks: midranks, rankTestP: rankTestP, rankTestPlan: rankTestPlan,
      permanova: permanova, brayCurtis: brayCurtis, bh: bh, stars: stars, choose: choose,
    };
  }

  var root = typeof document === "undefined" ? null : document.getElementById("etol-pie");
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

  // Domains present in a set of columns, largest total first (stack + legend order).
  function orderedDomains(cols) {
    var overall = domainTotals(cols);
    return Object.keys(overall).sort(function (a, b) { return overall[b] - overall[a]; });
  }

  // --- Significance ----------------------------------------------------------
  //
  // A bar segment is one pooled number; a p-value needs the spread that pooling
  // threw away. So the test ignores the drawn bar and goes back to the per-sample
  // values behind it: group them by design-matrix condition, compare the two
  // groups per domain with an exact Wilcoxon rank-sum test (read counts are
  // skewed, zero-inflated and n is small, so no distributional assumption is
  // safe), Benjamini-Hochberg corrected across the domains in the legend.
  //
  // Refused rather than approximated:
  //   - raw matched reads: not depth-normalized, so a "difference" can just be a
  //     difference in library size. No test on that axis.
  //   - the chi-square/Fisher-on-the-stack that a stacked bar invites: reads
  //     within a sample are not independent draws, so it returns p ~ 0 for any
  //     two libraries and means nothing.
  // ponytail: two groups only; add Kruskal-Wallis if a 3-condition design shows up.
  function domainStats(cols, domains) {
    var mode = barValue();
    if (mode === "hits") {
      return { note: "Raw counts are not depth-normalized — comparing them across samples is " +
        "unreliable, so no significance test is offered on this axis." };
    }
    var groups = {};
    cols.forEach(function (c) {
      var cond = String(matrix.cols[c].condition || "").trim();
      if (cond) (groups[cond] = groups[cond] || []).push(c);
    });
    var names = Object.keys(groups);
    if (names.length < 2) {
      return { note: "No significance test: needs at least two design-matrix conditions." };
    }

    // Per-sample value on the displayed axis; null when the sample cannot be
    // normalized (no control-probe reads / no reads at all) so it drops out of
    // the group rather than counting as a zero it never measured.
    var perCol = matrix.cols.map(function (_, c) { return domainTotals([c]); });
    function valueOf(c, domain) {
      var totals = perCol[c];
      var raw = totals[domain] || 0;
      if (mode === "pct") {
        var grand = Object.keys(totals).reduce(function (s, d) { return s + totals[d]; }, 0);
        return grand > 0 ? (raw / grand) * 100 : null;
      }
      var cells = matrix.cols[c].host_cells || 0;
      return cells > 0 ? raw / cells : null;
    }
    function valuesFor(name, domain) {
      return groups[name].map(function (c) { return valueOf(c, domain); })
        .filter(function (v) { return v !== null; });
    }

    var rows = domains.map(function (d) {
      return { domain: d, values: names.map(function (name) { return valuesFor(name, d); }) };
    }).filter(function (r) {
      return r.values.every(function (v) { return v.length >= 3; });
    });
    if (!rows.length) {
      return { note: "No significance test: needs at least 3 normalizable samples per condition." };
    }

    var ps = rows.map(function (r) { return rankTestP(r.values); });
    var qs = bh(ps);
    var by = {};
    rows.forEach(function (r, i) {
      by[r.domain] = {
        p: ps[i], q: qs[i], stars: stars(qs[i]),
        n: r.values.map(function (v) { return v.length; }),
        medians: r.values.map(median),
      };
    });

    // Smallest q a single differing domain could reach: the test's own p floor
    // after correcting over the tested domains. Above 0.05 the design is
    // underpowered and an unmarked domain means "undetectable at any effect
    // size", not "no difference" -- the caption has to say which, or a null
    // result gets read as a negative finding.
    var plans = rows.map(function (r) {
      return rankTestPlan(r.values.map(function (v) { return v.length; }));
    });
    var floorQ = rows.length * Math.min.apply(null, plans.map(function (p) { return p.floor; }));

    // Global test: one p for the whole community, on the same per-sample vectors
    // over the same domains, so it cannot disagree with the per-domain results.
    var vectors = [], labels = [];
    names.forEach(function (name) {
      groups[name].forEach(function (c) {
        var vector = domains.map(function (d) { return valueOf(c, d); });
        if (vector.every(function (v) { return v !== null; })) {
          vectors.push(vector);
          labels.push(name);
        }
      });
    });

    return {
      groups: names, by: by, axis: VALUE_LABELS[mode],
      counts: rows[0].values.map(function (v) { return v.length; }),
      tested: rows.length, testName: plans[0].name, floorQ: floorQ,
      global: permanova(vectors, labels),
    };
  }

  function buildBarsSvg(scope) {
    var series = barSeries(scope);
    var data = series.map(barTotals);

    // One stack order and one legend for every bar: domains ranked by their total
    // over the whole chart, largest at the bottom.
    var everyCol = series.reduce(function (acc, s) { return acc.concat(s.cols); }, []);
    var domains = orderedDomains(everyCol);
    var stats = domainStats(everyCol, domains);

    var maxTotal = data.reduce(function (m, d) { return Math.max(m, d.total); }, 0);
    var tickVals = ticks(maxTotal);
    var top = tickVals[tickVals.length - 1] || 1;

    var plotW = series.length * (BAR.w + BAR.gap) + BAR.gap;
    var width = BAR.left + plotW + BAR.right;
    var baseY = BAR.top + BAR.h;

    // Caption under the axis: either why there is no test, or which test the
    // legend's stars come from. Drawn into the SVG, so the PNG/SVG exports carry
    // it without the reader having to find this page again.
    var foot;
    if (stats.note) {
      foot = [stats.note];
    } else {
      foot = [
        stats.testName + " per domain on per-sample " + stats.axis + ": " +
          stats.groups.map(function (name, i) { return name + " (n=" + stats.counts[i] + ")"; })
            .join(", ") + "; Benjamini-Hochberg FDR over " + stats.tested + " domains.",
        "*** q<0.001    ** q<0.01    * q<0.05    unmarked = not significant. " +
          "Replicates are samples, not reads.",
      ];
      if (stats.global) {
        foot.splice(1, 0, "Whole-community test — PERMANOVA on Bray-Curtis composition: " +
          "pseudo-F = " + num(stats.global.f) + ", R² = " + num(stats.global.r2) +
          ", p = " + num(stats.global.p) + " (" + stats.global.permutations + " permutations).");
      }
      if (stats.floorQ > 0.05) {
        foot.push("Underpowered: at these group sizes a perfectly separated domain reaches only " +
          "q = " + num(stats.floorQ) + ", so unmarked means undetectable, not absent. " +
          "Needs about 6 samples per condition.");
      }
    }
    var height = baseY + BAR.bottom + foot.length * 14;

    var parts = ['<rect width="' + width + '" height="' + height + '" fill="#ffffff"/>'];
    var title = (matrix.preset_label || "eToL") + " — non-human reads by domain";
    parts.push('<text x="16" y="24" font-size="14" font-weight="600">' + esc(title) + "</text>");
    parts.push('<text x="16" y="42" font-size="11" fill="#555">' +
      esc(scope.subtitle || scope.label) + " — " + esc(VALUE_LABELS[barValue()]) + "</text>");

    if (!domains.length) {
      parts.push('<text x="16" y="80" font-size="12" fill="#666">No non-human probe-matched reads to plot.</text>');
      return svgWrap(width, baseY + BAR.bottom, parts.join(""));
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
        parts.push('<text x="' + lx + '" y="' + (baseY + 90) + '" font-size="9" fill="#777" ' +
          'text-anchor="middle">' + esc(s.condition) + "</text>");
      }
    });

    // Legend on the right, same order as the stack; a tested domain carries its
    // significance stars here rather than as brackets over the bars -- with one
    // test per segment, brackets over a stack would be ambiguous about which
    // segment they mark.
    var lgx = BAR.left + plotW + 20;
    domains.forEach(function (d, i) {
      var y = BAR.top + 10 + i * 22;
      var mark = (stats.by && stats.by[d] && stats.by[d].stars) || "";
      parts.push('<rect x="' + lgx + '" y="' + (y - 11) + '" width="13" height="13" fill="' + colorFor(d) + '"/>');
      parts.push('<text x="' + (lgx + 19) + '" y="' + y + '" font-size="11" fill="#333">' + esc(d) +
        (mark ? '<tspan font-weight="700">&#160;&#160;' + mark + "</tspan>" : "") + "</text>");
    });

    foot.forEach(function (line, i) {
      var warn = stats.note || i === foot.length - 1 && stats.floorQ > 0.05;
      parts.push('<text x="16" y="' + (baseY + BAR.bottom + 2 + i * 14) + '" font-size="10" fill="' +
        (warn ? "#a15c00" : "#555") + '">' + esc(line) + "</text>");
    });
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

    // The numbers behind the figure's stars, on the axis currently selected --
    // stars alone are not citable. Appended as a second block so it stays one
    // download; blank line separates it for a spreadsheet import.
    var everyCol = matrix.cols.map(function (_, i) { return i; });
    var stats = domainStats(everyCol, orderedDomains(everyCol));
    lines.push("");
    if (stats.note) {
      lines.push(csvCell(stats.note));
    } else {
      lines.push(csvCell(stats.testName + " per domain on per-sample " + stats.axis +
        ", Benjamini-Hochberg FDR over " + stats.tested + " domains"));
      if (stats.global) {
        lines.push(csvCell("PERMANOVA (Bray-Curtis, " + stats.global.permutations +
          " permutations): pseudo-F=" + num(stats.global.f) + ", R2=" + num(stats.global.r2) +
          ", p=" + num(stats.global.p)));
      }
      // Long format: one row per domain x group, p/q repeated. Keeps a single
      // table shape whether the design has two conditions or ten.
      lines.push(["Domain", "Group", "n", "Median", "p", "q (BH)", "Significance"]
        .map(csvCell).join(","));
      Object.keys(stats.by).forEach(function (d) {
        var s = stats.by[d];
        stats.groups.forEach(function (name, i) {
          lines.push([d, name, s.n[i], num(s.medians[i]),
            s.p.toExponential(3), s.q.toExponential(3), s.stars || "ns"].map(csvCell).join(","));
        });
      });
    }
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
