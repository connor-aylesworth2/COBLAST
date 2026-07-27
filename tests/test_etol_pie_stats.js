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
const {
  midranks, rankTestP, rankTestPlan, permanova, brayCurtis, bh, stars, choose,
} = require("../static/etol_pie.js");

const near = (a, b, msg) => assert.ok(Math.abs(a - b) < 1e-9, `${msg}: ${a} != ${b}`);

// Midranks: ties share the mean of the ranks they span.
assert.deepStrictEqual(midranks([10, 20, 30]), [1, 2, 3]);
assert.deepStrictEqual(midranks([5, 5, 9]), [1.5, 1.5, 3]);
assert.deepStrictEqual(midranks([0, 0, 0, 0]), [2.5, 2.5, 2.5, 2.5]);

// --- two groups: the exact enumerated branch, equal to two-sided Wilcoxon ----
// Complete separation, 4 vs 4: 2 / C(8,4) = 2/70.
near(rankTestP([[1, 2, 3, 4], [5, 6, 7, 8]]), 2 / 70, "4v4 separated");
// Direction must not matter (the rank statistic is two-sided by construction).
near(rankTestP([[5, 6, 7, 8], [1, 2, 3, 4]]), 2 / 70, "4v4 reversed");
// 3 vs 3 cannot reach 0.05 even when perfectly separated: 2 / C(6,3) = 0.1.
near(rankTestP([[1, 2, 3], [4, 5, 6]]), 0.1, "3v3 floor");
// Unequal group sizes: 2 vs 4 separated is 2 / C(6,2) = 2/15.
near(rankTestP([[1, 2], [3, 4, 5, 6]]), 2 / 15, "2v4 separated");
// No difference at all, and the all-tied case (every domain absent), are p = 1.
near(rankTestP([[1, 2, 3, 4], [1, 2, 3, 4]]), 1, "identical groups");
near(rankTestP([[0, 0, 0], [0, 0, 0]]), 1, "all ties");
assert.strictEqual(rankTestP([[], [1, 2]]), null, "empty group");
// Sampled branch (C(30,15) > 100k) still lands near the analytic answer.
const big = (n, off) => Array.from({ length: n }, (_, i) => i + off);
assert.ok(rankTestP([big(15, 0), big(15, 100)]) < 0.001, "large separated groups");
assert.ok(rankTestP([big(15, 0), big(15, 0)]) > 0.9, "large identical groups");

// --- more than two groups: the same call is Kruskal-Wallis -------------------
// 4 groups of 6, cleanly ordered: far below any threshold. 4 groups of 6 with
// identical values must not manufacture a signal.
const g4 = (off) => big(6, off);
assert.ok(rankTestP([g4(0), g4(100), g4(200), g4(300)]) < 0.001, "4x6 separated");
assert.ok(rankTestP([g4(0), g4(0), g4(0), g4(0)]) > 0.9, "4x6 identical");
// One group shifted, three overlapping, is still detectable at 4x6.
assert.ok(rankTestP([g4(0), g4(1), g4(2), g4(500)]) < 0.01, "4x6 one group apart");

// --- what the test can resolve, behind the "underpowered" caption -----------
assert.strictEqual(choose(6, 3), 20);
assert.strictEqual(choose(10, 5), 252);
near(rankTestPlan([3, 3]).floor, 0.1, "3v3 floor");
near(rankTestPlan([6, 6]).floor, 2 / 924, "6v6 floor");
assert.ok(rankTestPlan([3, 3]).floor * 7 > 0.05, "3v3 over 7 domains is unreachable");
assert.ok(rankTestPlan([5, 5]).floor * 7 > 0.05, "5v5 over 7 domains is unreachable");
assert.ok(rankTestPlan([6, 6]).floor * 7 < 0.05, "6v6 over 7 domains is reachable");
// 4x6 is sampled, so the floor is the sampling resolution, not combinatorics.
assert.ok(rankTestPlan([6, 6, 6, 6]).floor * 7 < 0.001, "4x6 over 7 domains is ample");
assert.strictEqual(rankTestPlan([6, 6]).name, "Exact Wilcoxon rank-sum test");
assert.strictEqual(rankTestPlan([6, 6, 6, 6]).name, "Kruskal-Wallis permutation test");

// --- PERMANOVA: the whole-community question --------------------------------
near(brayCurtis([1, 0], [1, 0]), 0, "identical composition");
near(brayCurtis([1, 0], [0, 1]), 1, "disjoint composition");
near(brayCurtis([0, 0], [0, 0]), 0, "empty composition");
// Two compositions, 6 samples each, cleanly distinguishable -> small p, high R2.
const compA = Array.from({ length: 6 }, (_, i) => [90 + i, 10 - i, 5]);
const compB = Array.from({ length: 6 }, (_, i) => [10 + i, 90 - i, 5]);
const labels2 = [...Array(6).fill("AD"), ...Array(6).fill("CTRL")];
const sep = permanova([...compA, ...compB], labels2);
assert.ok(sep.p < 0.01, `separated composition p=${sep.p}`);
assert.ok(sep.r2 > 0.8, `separated composition R2=${sep.r2}`);
// Same composition in both groups must not come out significant.
const same = permanova([...compA, ...compA], labels2);
assert.ok(same.p > 0.2, `identical composition p=${same.p}`);
assert.strictEqual(permanova([[1, 2]], ["AD"]), null, "too few samples");

// Benjamini-Hochberg, in input order, monotone after the step-up.
assert.deepStrictEqual(bh([0.01, 0.02, 0.03]).map((q) => +q.toFixed(6)), [0.03, 0.03, 0.03]);
assert.deepStrictEqual(bh([0.001, 0.5]).map((q) => +q.toFixed(6)), [0.002, 0.5]);
assert.deepStrictEqual(bh([0.6, 0.7]).map((q) => +q.toFixed(6)), [0.7, 0.7]); // capped at 1

assert.strictEqual(stars(0.0009), "***");
assert.strictEqual(stars(0.009), "**");
assert.strictEqual(stars(0.049), "*");
assert.strictEqual(stars(0.051), "");

console.log("ok - etol_pie stats");
