/*
 * Check the stacked-bar significance maths in static/etol_pie.js.
 *
 *     node tests/test_etol_pie_stats.js
 *
 * The expected p-values are the exact Wilcoxon rank-sum two-sided values, which
 * for these sizes are just counted by hand: with complete separation only the
 * two extreme label assignments are as extreme as the observed one, so
 * p = 2 / C(n1+n2, n1).
 */
const assert = require("assert");
const { midranks, wilcoxonP, bh, stars, choose } = require("../static/etol_pie.js");

const near = (a, b, msg) => assert.ok(Math.abs(a - b) < 1e-9, `${msg}: ${a} != ${b}`);

// Midranks: ties share the mean of the ranks they span.
assert.deepStrictEqual(midranks([10, 20, 30]), [1, 2, 3]);
assert.deepStrictEqual(midranks([5, 5, 9]), [1.5, 1.5, 3]);
assert.deepStrictEqual(midranks([0, 0, 0, 0]), [2.5, 2.5, 2.5, 2.5]);

// Complete separation, 4 vs 4: 2 / C(8,4) = 2/70.
near(wilcoxonP([1, 2, 3, 4], [5, 6, 7, 8]), 2 / 70, "4v4 separated");
// Direction must not matter (two-sided).
near(wilcoxonP([5, 6, 7, 8], [1, 2, 3, 4]), 2 / 70, "4v4 reversed");
// 3 vs 3 cannot reach 0.05 even when perfectly separated: 2 / C(6,3) = 0.1.
near(wilcoxonP([1, 2, 3], [4, 5, 6]), 0.1, "3v3 floor");
// No difference at all, and the all-tied case (every domain absent), are p = 1.
near(wilcoxonP([1, 2, 3, 4], [1, 2, 3, 4]), 1, "identical groups");
near(wilcoxonP([0, 0, 0], [0, 0, 0]), 1, "all ties");
assert.strictEqual(wilcoxonP([], [1, 2]), null, "empty group");
// Sampled branch (C(30,15) > 100k) still lands near the analytic answer.
const big = (n, off) => Array.from({ length: n }, (_, i) => i + off);
assert.ok(wilcoxonP(big(15, 0), big(15, 100)) < 0.001, "large separated groups");
assert.ok(wilcoxonP(big(15, 0), big(15, 0)) > 0.9, "large identical groups");

// Benjamini-Hochberg, in input order, monotone after the step-up.
assert.deepStrictEqual(bh([0.01, 0.02, 0.03]).map((q) => +q.toFixed(6)), [0.03, 0.03, 0.03]);
assert.deepStrictEqual(bh([0.001, 0.5]).map((q) => +q.toFixed(6)), [0.002, 0.5]);
assert.deepStrictEqual(bh([0.6, 0.7]).map((q) => +q.toFixed(6)), [0.7, 0.7]); // capped at 1

// Power floor behind the chart's "underpowered" caption: the best q a single
// separated domain can reach is 2 / C(n1+n2, n1) times the number tested.
assert.strictEqual(choose(6, 3), 20);
assert.strictEqual(choose(10, 5), 252);
assert.ok((2 / choose(6, 3)) * 7 > 0.05, "3v3 over 7 domains is unreachable");
assert.ok((2 / choose(10, 5)) * 7 > 0.05, "5v5 over 7 domains is unreachable");
assert.ok((2 / choose(12, 6)) * 7 < 0.05, "6v6 over 7 domains is reachable");

assert.strictEqual(stars(0.0009), "***");
assert.strictEqual(stars(0.009), "**");
assert.strictEqual(stars(0.049), "*");
assert.strictEqual(stars(0.051), "");

console.log("ok - etol_pie stats");
