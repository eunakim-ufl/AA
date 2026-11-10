# Assignment 2: Dynamic Programming

Eun A Kim (GatorID: eunakim, UFID: 42686893)

All code is C++17 (`g++ -O2 -std=c++17`). Every algorithm is iterative (no recursion) and stores its tables in contiguous arrays (`std::vector`, a dynamic array).

**Test machine.** CPU: Intel Xeon E5-2697 v4 @ 2.30 GHz, OS: Debian Linux (kernel 5.10.0-37-amd64), compiler: g++ 10.2.1 (Debian 10.2.1-6).


| File | Content |
|---|---|
| `problem1.cpp` | Problem 1: DP solvers, brute-force checker, data generator, experiments |
| `problem2.cpp` | Problem 2: DP solver, prefix-sum checker, data generator, experiments |
| `plot_p1.py`, `plot_p2.py` | Plotting scripts |
| `results_p1.csv`, `results_p2.csv` | Running-time experiments |
| `mem_p1.csv`, `mem_p2.csv` | Memory experiments (one process per configuration) |


```bash
# compile
g++ -O2 -std=c++17 problem1.cpp -o p1
g++ -O2 -std=c++17 problem2.cpp -o p2

# running-time experiments
./p1 > results_p1.csv
./p2 > results_p2.csv

# plots
python plot_p1.py
python plot_p2.py
```

---

## Problem 1: Weighted Approximate Common Substring

### 1.1 Problem

X = x<sub>1</sub>…x<sub>n</sub> and Y = y<sub>1</sub>…y<sub>m</sub> are strings over [A–Z]. We choose two substrings of the same length L, X[i−L+1..i] and Y[j−L+1..j], and compare them position by position (no gaps). Each position scores

```
s(a, b) = w_a      if a = b
        = -delta   if a != b
```

and the score of the pair is the sum over its L positions. We want the pair with the largest score. The empty pair (L = 0) scores 0, so the optimum is never negative.

Because there are no gaps, the k-th character of one substring is always compared with the k-th character of the other. In the n × m grid of all comparisons (x<sub>i</sub>, y<sub>j</sub>), a candidate pair is therefore a contiguous segment of one diagonal (j − i constant). The problem becomes a maximum-sum contiguous segment problem on each of the n + m − 1 diagonals, and the DP below solves all of them at once.

### 1.2 Optimization function

For 0 ≤ i ≤ n and 0 ≤ j ≤ m, `D[i][j]` is the largest score of a pair (possibly empty) that **ends** at x<sub>i</sub> and y<sub>j</sub>.

### 1.3 Recurrence (Bellman equations)

```
D[i][0] = D[0][j] = 0
D[i][j] = max( 0,  D[i-1][j-1] + s(x_i, y_j) )      for 1 <= i <= n, 1 <= j <= m
OPT     = max over all i, j of D[i][j]
```

### 1.4 Correctness

We prove by induction on i + j that `D[i][j]` is the best score of a pair ending at (i, j).

If i = 0 or j = 0, only the empty pair ends there, so the value is 0.

Otherwise, a pair ending at (i, j) is either empty (score 0) or non-empty. A non-empty pair is a pair ending at (i−1, j−1), possibly empty, followed by the comparison of x<sub>i</sub> with y<sub>j</sub>. Its score is (score of that shorter pair) + s(x<sub>i</sub>, y<sub>j</sub>). The term s(x<sub>i</sub>, y<sub>j</sub>) is the same for every such pair, so the score is maximized by taking the best pair ending at (i−1, j−1). By the induction hypothesis that best score is `D[i-1][j-1]`. Taking the larger of the empty pair and the best non-empty pair gives the recurrence.

This is the optimal-substructure property: the prefix of an optimal pair ending at (i, j) is itself optimal among pairs ending at (i−1, j−1). If it were not, swapping in a better prefix would raise the total score, which contradicts optimality.

Every pair ends at some cell (or is empty), so the largest `D[i][j]` is the optimal score. The cells are filled row by row, so `D[i-1][j-1]` is always ready when needed.

Handout example: "CAABC" (from position 3 of X) against "CAACC" (from position 4 of Y) gives w<sub>C</sub> + w<sub>A</sub> + w<sub>A</sub> − δ + w<sub>C</sub> = 2w<sub>C</sub> + 2w<sub>A</sub> − δ. With the Scenario 2 weights and δ = 5814 the program prints 2·2782 + 2·8167 − 5814 = 16084.

### 1.5 Extracting the optimal pair

While filling the table we keep the cell `(i*, j*)` with the largest value. With ties, the first one reached in row order is kept.

If `D[i][j] > 0`, the maximum in the recurrence was reached by the second term. The best pair ending at (i, j) is then the best pair ending at (i−1, j−1) plus one position. So starting from `(i*, j*)` we step to `(i*-1, j*-1)`, `(i*-2, j*-2)`, … while D > 0, and count the steps as L. The first cell with D = 0 is where the optimal pair starts over. The answer is `X[i*-L+1 .. i*]` and `Y[j*-L+1 .. j*]`, with score `D[i*][j*]`.

Two-row version: with each value we also store the length of its pair, `R[i][j] = R[i-1][j-1] + 1` if `D[i][j] > 0` and 0 otherwise. Only two rows of D and R are kept, and R at the best cell gives L directly, so extraction does not need the full table. On every input, the two versions returned the same pair.

### 1.6 Complexity

**Time Θ(nm).** There are (n+1)(m+1) cells, each computed with a constant number of operations. The walk back takes at most min(n, m) steps.

| Version | Arrays | Space |
|---|---|---|
| Full table | (n+1)(m+1) four-byte integers | Θ(nm) |
| Two-row | 4 arrays of m+1 four-byte integers | Θ(m); Θ(min(n, m)) if Y is the shorter string |
| Input | n + m bytes | Θ(n + m) |

No score can exceed max<sub>l</sub> w<sub>l</sub> · min(n, m) ≤ 12702 · 8000 = 101,616,000 < 2<sup>31</sup> − 1, so `int` cannot overflow. The program checks this bound before each run.

### 1.7 Validation

`solveFull` (full table + walk back) and `solveLinear` (two rows) were compared with `bruteForce`, which tries every start pair and every length in O(nm · min(n, m)) time. The test used 500 random cases with n, m ≤ 60, covering both scenarios and random δ. Half of the cases used a 3-letter alphabet so that long pairs containing mismatches occur often. All 500 passed.

In every experiment the program also recomputed the score of the reported pair from scratch (`score_ok`) and compared the two solvers (`full_eq_lin`). All checks passed.

Handout strings: in Scenario 1 the optimum is "BCAA" (position 2 in X, position 3 in Y, score 4). "BCAA" also appears at position 6 of X, so the optimum is not unique, and the program reports the first one. In Scenario 2 with δ = 5814 the optimum is "AABCAA" / "ABBCAA" with score 3·8167 + 1492 + 2782 − 5814 = 22961. It contains one mismatch, because the high weights on the surrounding matches outweigh the penalty.

### 1.8 Experimental setup

**Scenario 1:** w<sub>l</sub> = 1 for all letters, δ = 10. One mismatch cancels ten matches, so including a mismatch only helps if it brings in more than ten extra matches.

**Scenario 2:** w<sub>l</sub> is proportional to the frequency of letter l in English. We used the standard table by Lewand (*Cryptological Mathematics*), in units of 0.001 %:

| Letter | A | B | C | D | E | F | G | H | I | J | K | L | M |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| w | 8167 | 1492 | 2782 | 4253 | 12702 | 2228 | 2015 | 6094 | 6966 | 153 | 772 | 4025 | 2406 |

| Letter | N | O | P | Q | R | S | T | U | V | W | X | Y | Z |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| w | 6749 | 7507 | 1929 | 95 | 5987 | 6327 | 9056 | 2758 | 978 | 2360 | 150 | 1974 | 74 |

Scaling all weights and δ by the same factor c > 0 scales every score by c, so the optimal pair depends only on the ratios. Integer units keep all scores exact. The 10 intermediate penalties split [w<sub>min</sub>, w<sub>max</sub>] = [74, 12702] into 11 equal steps:

```
delta_t = 74 + t * (12702 - 74) / 11 = 74 + 1148 t,     t = 1..10
        = 1222, 2370, 3518, 4666, 5814, 6962, 8110, 9258, 10406, 11554
```

The endpoints δ = 74 (t = 0) and δ = 12702 (t = 11) were also run.

**Data:** random strings with letters drawn either uniformly from the 26 letters or from the English distribution above. Sizes were n = m ∈ {100, 500, 1000, 2000, 4000, 8000}, plus one unequal case, 1000 × 4000. The seed is fixed (20261103). The full table was run up to 4000 × 4000 (64 MB); for 8000 × 8000 only the two-row version was run.

**Measurement:** time with `std::chrono::steady_clock`, using one warm-up run and then repetitions until at least 200 ms and 3 runs had passed. The reported number is the average per run. Memory was measured in a separate process for each configuration, in two ways:

- **heap peak:** bytes allocated by the solver, counted exactly by overriding `operator new`;
- **ΔRSS:** change in the resident memory reported by the OS (`/proc/self/statm`, `getrusage`).

### 1.9 Results

**Table 1.1** Running time, Scenario 1, English letters.

| n × m | cells | full table (µs) | two-row (µs) | two-row, ns per cell |
|---|---|---|---|---|
| 100 × 100 | 10<sup>4</sup> | 18.80 | 22.72 | 2.27 |
| 500 × 500 | 2.5·10<sup>5</sup> | 527.22 | 519.56 | 2.08 |
| 1000 × 1000 | 10<sup>6</sup> | 2353.20 | 2415.01 | 2.42 |
| 2000 × 2000 | 4·10<sup>6</sup> | 9132.71 | 8558.43 | 2.14 |
| 1000 × 4000 | 4·10<sup>6</sup> | 9019.38 | 9455.98 | 2.36 |
| 4000 × 4000 | 1.6·10<sup>7</sup> | 43435.57 | 34973.76 | 2.19 |
| 8000 × 8000 | 6.4·10<sup>7</sup> | — | 150102.80 | 2.35 |

**Table 1.2** Memory. Theory = DP arrays + input strings (n + m bytes); heap peak = DP arrays only.

| n × m | full: theory (B) | full: heap peak (B) | two-row: theory (B) | two-row: heap peak (B) | full: ΔRSS (KB) |
|---|---|---|---|---|---|
| 100 × 100 | 41,004 | 40,804 | 1,816 | 1,616 | n/r |
| 500 × 500 | 1,005,004 | 1,004,004 | 9,016 | 8,016 | n/r |
| 1000 × 1000 | 4,010,004 | 4,008,004 | 18,016 | 16,016 | 5,256 |
| 2000 × 2000 | 16,020,004 | 16,016,004 | 36,016 | 32,016 | 17,064 |
| 4000 × 4000 | 64,040,004 | 64,032,004 | 72,016 | 64,016 | 63,800 |
| 8000 × 8000 | — | — | 144,016 | 128,016 | — |

n/r = not resolvable: the array is too small to raise the process's start-up peak (see 1.10).

**Table 1.3** Scenario 2, δ sweep, n = m = 2000, English letters.

| t | δ | length | mismatches | score | optimal pair |
|---|---|---|---|---|---|
| 0 | 74 | 1852 | 1702 | 1,155,080 | most of a diagonal |
| 1 | 1222 | 74 | 57 | 105,713 | (74 letters) |
| 2 | 2370 | 18 | 10 | 63,360 | ESEOGAEBEATHAERRET / ESECLRESEAANFEGHIT |
| 3 | 3518 | 14 | 7 | 53,378 | ESEOGAEBEATHAE / ESECLRESEAANFE |
| 4–6 | 4666–6962 | 7 | 1 | 51,651 / 50,503 / 49,355 | EEYTMEN / EETTMEN |
| 7–11 | 8110–12702 | 5 | 0 | 48,294 | ESOTE / ESOTE |

Uniform letters, same sizes: δ = 74 gave length 1635 with 1543 mismatches; δ = 1222 and 2370 gave the 8-letter pair EEBQFJJE / EEBKUGRE with 4 mismatches; every δ ≥ 3518 gave the exact match "EEB" (2·12702 + 1492 = 26,896).

**Table 1.4** Scenario 1, length of the optimal pair (always 0 mismatches).

| n = m | 100 | 500 | 1000 | 2000 | 4000 | 8000 |
|---|---|---|---|---|---|---|
| uniform | 2 | 3 | 3 | 4 | 4 | 5 |
| English | 3 | 4 | 5 | 5 | 6 | 6 |

![Figure 1](p1_plots.png)

**Figure 1.** Left: running time vs. n·m (log–log). Middle: memory vs. n·m; ΔRSS is shown only for the full table, because the two-row arrays are too small to change RSS. Right: length and number of mismatches of the optimal pair vs. δ (n = m = 2000); a mismatch count of 0 is plotted at 0.5 because 0 cannot be shown on a log axis.

### 1.10 Analysis

**Time.** Over a 6400-fold increase in n·m, the time per cell stays between 2.1 and 2.4 ns (Table 1.1), and the log–log plot is a line of slope 1. This agrees with Θ(nm). The 2000 × 2000 and 1000 × 4000 inputs have the same number of cells and similar times, so the cost depends on n·m, not on the shape. The full-table and two-row versions run at similar speed. Their difference stays within about 20 % and its direction is not consistent across experiments, so we draw no conclusion from it.

**Memory.** The heap peak equals the size of the DP arrays exactly. For example, (4000+1)<sup>2</sup> · 4 = 64,032,004 bytes for the full table versus 4 · 4001 · 4 = 64,016 bytes for the two-row version, a factor of 1000. For the full table, ΔRSS minus the heap peak is 1342, 1423 and 1269 KB for the three largest sizes. This is a roughly constant overhead whose source we did not identify. Between 2000 × 2000 and 4000 × 4000, ΔRSS grows by 46,736 KB against 46,891 KB for the array, a 99.7 % match. For small arrays ΔRSS shows nothing: `ru_maxrss` is the highest RSS the process has ever had, about 4.2 MB right after start-up, and arrays that do not exceed that peak leave it unchanged. For the two-row version, the heap peak is therefore the meaningful measurement.

**Scenario 1.** No optimum contained a mismatch, so in practice the answer is the longest exact common substring. Its length grows very slowly (Table 1.4). Let q be the probability that two random letters are equal: q = 1/26 ≈ 0.038 for uniform letters and q = Σ p<sub>l</sub><sup>2</sup> ≈ 0.0655 for English. A rough estimate of the longest exact match is log<sub>1/q</sub>(nm). It gives 4.2 (uniform) and 5.1 (English) at 1000 × 1000, and 5.5 and 6.6 at 8000 × 8000. The observed lengths are 0 to 1.3 below these values, so the estimate is slightly high but follows the logarithmic growth. Because English has the larger q, its matches are longer.

**Scenario 2.** The change with δ can be understood from the expected score of one step along a diagonal of two random strings:

```
E[s] = sum_l p_l^2 * w_l  -  delta * (1 - sum_l p_l^2)
```

If E[s] > 0, the running sum tends to grow, and the best pair can cover most of a diagonal. If E[s] < 0, it tends to fall, and only short segments survive. Setting E[s] = 0 gives the threshold δ<sub>c</sub> ≈ 154 for uniform letters and δ<sub>c</sub> ≈ 573 for English letters (calculated from the weight table). δ = 74 is below both thresholds, and indeed the optimum covered 1852 (English) and 1635 (uniform) of the 2000 positions. All ten intermediate values are above both thresholds and gave short pairs. This E[s] argument is a heuristic explanation, not a proof.

The number of mismatches in the optimum can never increase when δ increases. Let A be optimal for δ<sub>1</sub> and B optimal for δ<sub>2</sub> > δ<sub>1</sub>, with total match weights W<sub>A</sub>, W<sub>B</sub> and mismatch counts k<sub>A</sub>, k<sub>B</sub>. Optimality gives

```
W_A - k_A * delta1  >=  W_B - k_B * delta1
W_B - k_B * delta2  >=  W_A - k_A * delta2
```

Adding the two inequalities and cancelling W<sub>A</sub> + W<sub>B</sub> leaves (k<sub>A</sub> − k<sub>B</sub>)(δ<sub>2</sub> − δ<sub>1</sub>) ≥ 0, so k<sub>A</sub> ≥ k<sub>B</sub>. The data agree: 1702 → 57 → 10 → 7 → 1 → 0 (English) and 1543 → 4 → 0 (uniform). The length of the optimum is not covered by this argument, although it also decreased here.

When the optimal pair stays the same between two δ values, its score changes by −(mismatches) × Δδ. For t = 4 → 5 → 6 (1 mismatch) the score drops by exactly 1148 each step. For the uniform strings at t = 1 → 2 (4 mismatches) it drops by 4592 = 4 · 1148. Once the optimum has no mismatches, its score no longer depends on δ, and it stays fixed for all larger δ (t ≥ 7). The surviving pairs are made of high-weight letters (E, T, S, O), which is what frequency-based weights should produce.

δ also changed the running time by a constant factor. The two-row solver took 22.6 ms at δ = 1222 and 8.7 ms at δ = 12702, while the full-table solver stayed between 8.0 and 9.2 ms. We suspect branch misprediction in the two-row loop, but did not verify it. The asymptotic cost is Θ(nm) in every case.

---




## Problem 2: Largest All-Zero Square Sub-matrix

### 2.1 Problem

Given an m × n Boolean matrix B (non-zero entries stored as 1), find the largest k such that B contains a k × k block of consecutive rows and columns that is all 0, and return its position.

### 2.2 Optimization function

`D[i][j]` is the side length of the largest all-zero square whose **bottom-right corner** is (i, j).

### 2.3 Recurrence (Bellman equations)

```
D[i][j] = 0                                               if B[i][j] = 1
        = 1                                               if B[i][j] = 0 and (i = 0 or j = 0)
        = 1 + min( D[i-1][j], D[i][j-1], D[i-1][j-1] )    otherwise
k*      = max over all i, j of D[i][j]
```

The code adds a row and a column of zeros in front, so the boundary case is handled by the general formula.

### 2.4 Correctness

If B[i][j] = 1, no all-zero square can end at (i, j), so D = 0. In the boundary case with B[i][j] = 0, the only square ending at (i, j) is the cell itself, so D = 1.

Now let B[i][j] = 0 with i, j ≥ 1, and let t = min(D[i−1][j], D[i][j−1], D[i−1][j−1]). We show D[i][j] = t + 1.

*D[i][j] ≤ t + 1.* Take an all-zero k × k square ending at (i, j) with k ≥ 2. Removing its last row and its last column in different ways gives three all-zero (k−1) × (k−1) squares, ending at (i−1, j), (i, j−1) and (i−1, j−1). So each of the three neighbor values is at least k − 1. This gives t ≥ k − 1, i.e. k ≤ t + 1. For k = 1 the bound holds trivially.

*D[i][j] ≥ t + 1.* Each neighbor has an all-zero square of side at least t ending at it, so it also has one of side exactly t (shrink toward the corner). Consider the (t+1) × (t+1) square with rows i−t..i and columns j−t..j, and a cell (r, c) in it:
- (r, c) = (i, j): zero by assumption;
- r = i, c < j: covered by the t × t square ending at (i, j−1) (rows i−t+1..i, columns j−t..j−1);
- r < i, c = j: covered by the square ending at (i−1, j) (rows i−t..i−1, columns j−t+1..j);
- r < i, c < j: covered by the square ending at (i−1, j−1) (rows i−t..i−1, columns j−t..j−1).

So every cell is 0, and a (t+1) × (t+1) all-zero square ends at (i, j).

Every square has a bottom-right corner, so `k*` = max D is the answer. Row-major order computes the three neighbors before (i, j).

### 2.5 Extracting the optimal square

We keep `k*` and the first cell `(i*, j*)` where it occurs. The square is rows `i*-k*+1 .. i*`, columns `j*-k*+1 .. j*`. No walk back is needed, because the corner and the side length determine the square.

### 2.6 Complexity

**Time Θ(mn):** one test and one minimum of three values per cell.

**Space.** Following the hint, B is stored with one byte per entry (`uint8_t`). The DP values do not fit in a byte, since k can be as large as min(m, n) = 1000 > 255. Each row of D depends only on the row above, so we keep two rows of four-byte integers.

| Part | Size | Space |
|---|---|---|
| Matrix B | mn bytes | Θ(mn) |
| Two DP rows | 2(n+1) · 4 bytes | Θ(n) |
| Total | mn + 8(n+1) bytes | Θ(mn) |

The extra space is Θ(n), or Θ(min(m, n)) if the matrix is processed so that rows are the longer side.

### 2.7 Validation

`verify` checks every result independently with 2D prefix sums. It confirms that (1) the reported square lies inside B and contains only zeros, and (2) no all-zero (k*+1) × (k*+1) square exists anywhere in B. Together these prove the answer is optimal. All experiments passed. A hand-made 5 × 6 matrix with known answer k = 3 at top-left (0, 1) (0-indexed) also passed.

### 2.8 Experimental setup

The sizes were 10 × 10, 10 × 100, 10 × 1000, 100 × 1000 and 1000 × 1000. Each entry is set to 1 independently with probability p ∈ {0.5, 0.1, 0.01}. Several p values are used because at p = 0.5 the answer is only 2–4. Timing is the same as in Problem 1, except at least 5 runs are averaged. Memory was measured in a separate process for each configuration.

### 2.9 Results

**Table 2.1** Running time in µs, with k* in parentheses.

| m × n | cells | p = 0.5 | p = 0.1 | p = 0.01 |
|---|---|---|---|---|
| 10 × 10 | 10<sup>2</sup> | 0.196 (2) | 0.272 (3) | 0.309 (9) |
| 10 × 100 | 10<sup>3</sup> | 1.209 (2) | 2.524 (8) | 3.413 (10) |
| 10 × 1000 | 10<sup>4</sup> | 43.712 (3) | 26.194 (8) | 34.236 (10) |
| 100 × 1000 | 10<sup>5</sup> | 511.615 (4) | 292.002 (9) | 342.590 (33) |
| 1000 × 1000 | 10<sup>6</sup> | 5118.774 (4) | 2909.975 (10) | 3418.036 (30) |

**Table 2.2** Memory (the same for all p).

| m × n | theory (B) | heap peak (B) | ΔRSS (KB) |
|---|---|---|---|
| 10 × 10 | 188 | 188 | 0 |
| 10 × 100 | 1,808 | 1,808 | 0 |
| 10 × 1000 | 18,008 | 18,008 | 0 |
| 100 × 1000 | 108,008 | 108,008 | 0 |
| 1000 × 1000 | 1,008,008 | 1,008,008 | 2,292 – 2,416 (three runs) |

![Figure 2](p2_plots.png)

**Figure 2.** Left: running time vs. m·n for three values of p (log–log). Right: memory vs. m·n; theory and heap peak coincide. ΔRSS values of 0 are plotted at 1 KB because 0 cannot be shown on a log axis.

### 2.10 Analysis

**Time.** All three lines have slope about 1 on the log–log plot. For p = 0.1 the time per cell stays between 2.5 and 2.9 ns from 10<sup>2</sup> to 10<sup>6</sup> cells, which agrees with Θ(mn).

The content of the matrix changes only the constant factor. On the large matrices one cell costs about 5.1 ns at p = 0.5, 2.9 ns at p = 0.1 and 3.4 ns at p = 0.01. At p = 0.01 more cells take the zero branch, which does more work (a minimum of three values and a comparison with the best). At p = 0.5 the two branches are equally likely, so the CPU cannot predict them well. This explanation is our guess and was not measured. On the two smallest matrices p = 0.5 is the fastest, so the lines cross. These runs take under 4 µs, and the same tiny matrix is processed many times in a row, so they are not representative.

**Memory.** The heap peak equals mn + 8(n+1) bytes exactly for every size. ΔRSS is 0 up to 10<sup>5</sup> cells. Our guess is that these allocations were served from memory the process already had resident, or that Linux updates its RSS counters in batches; we did not confirm the cause. At 10<sup>6</sup> cells, ΔRSS exceeds the array size (984 KB) by 1308–1432 KB, close to the overhead seen in Problem 1.

**Size of k\*.** Setting the expected number of all-zero k × k squares, (m−k+1)(n−k+1)(1−p)<sup>k²</sup>, equal to 1 and approximating (m−k+1)(n−k+1) by mn gives a first-moment estimate:

```
k ≈ sqrt( ln(mn) / -ln(1-p) )
```

For 1000 × 1000 this gives 4.5, 11.4 and 37.1 for p = 0.5, 0.1 and 0.01. The observed values are 4, 10 and 30. The estimate is somewhat high, but it shows that k* grows only like the square root of log(mn).

---

## Conclusion

Both problems are solved by a table DP in which each cell takes constant time, so the running time is Θ(nm). In both cases the recurrence only reads the previous row, so the extra memory can be reduced to one dimension. In Problem 1, storing the pair length in each cell makes it possible to extract the answer without the full table. The experiments agree with the analysis: time per cell stays constant as the input grows, the measured heap memory equals the computed array sizes, and every answer passed an independent check (brute force for Problem 1, prefix sums for Problem 2).
