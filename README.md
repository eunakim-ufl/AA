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

X = x<sub>1</sub>…x<sub>n</sub> and Y = y<sub>1</sub>…y<sub>m</sub> are strings over [A–Z]. The two input strings may have different lengths (n = m and n ≠ m are both allowed). We choose two substrings of the same length L, X[i−L+1..i] and Y[j−L+1..j], and compare them position by position (since gaps are not allowed). 

Each position scores

```
s(a, b) = w_a      if a = b
        = -delta   if a != b
```

and the score of the pair is the sum over its L positions. We want the pair with the largest score. The empty pair (L = 0) scores 0, so the optimum is never negative.

Because there are no gaps, the k-th character of one substring is always compared with the k-th character of the other. In the n × m grid of all comparisons (x<sub>i</sub>, y<sub>j</sub>), a candidate pair is therefore a contiguous segment of one diagonal (j − i constant). The problem becomes a maximum-sum contiguous segment problem on each of the n + m − 1 diagonals, and the DP below solves all of them at once.

### 1.2 Optimization function

For 0 ≤ i ≤ n and 0 ≤ j ≤ m, `D[i][j]` is the largest score of a pair that ends at
position i of X and position j of Y, that is, X[i−L+1 .. i] against Y[j−L+1 .. j]
for some 0 ≤ L ≤ min(i, j):

```
D[i][j] = max over 0 <= L <= min(i,j) of  sum_{k=0}^{L-1} s(x_{i-k}, y_{j-k})
```

L = 0 is the empty pair with score 0, so `D[i][j] ≥ 0`.

### 1.3 Recurrence (Bellman equations)

```
D[i][0] = D[0][j] = 0
D[i][j] = max( 0,  D[i-1][j-1] + s(x_i, y_j) )      for 1 <= i <= n, 1 <= j <= m
OPT = max over all i, j of D[i][j]
```
OPT = Optimal Solution(score)


### 1.4 Correctness

We prove by induction on i + j that `D[i][j]` is the best score of a pair ending at (i, j).

If i = 0 or j = 0, only the empty pair ends there, so the value is 0.

Otherwise, a pair ending at (i, j) is either empty (score 0) or non-empty. A non-empty pair of length L is a pair of length L−1 ending at (i−1, j−1), possibly empty, followed
by the comparison of x<sub>i</sub> with y<sub>j</sub>. Its score is (score of the shorter pair) +
s(x<sub>i</sub>, y<sub>j</sub>). The last term is the same for all such pairs, so the best non-empty
pair uses the best pair ending at (i−1, j−1), whose score is `D[i-1][j-1]` by the
induction hypothesis. Taking the larger of 0 and `D[i-1][j-1] + s(x_i, y_j)` gives the recurrence.

This is the optimal-substructure property: if an optimal pair ending at (i, j) is
non-empty, removing its last position leaves an optimal pair ending at (i−1, j−1).
Otherwise a better prefix could be swapped in and the total score would increase.

Each `D[i][j]` is the score of an actual pair, so max D ≤ OPT. Conversely, every
non-empty pair ends at some cell (i, j) with i, j ≥ 1, so its score is at most `D[i][j]`,
and the empty pair scores 0 = `D[0][0]`. Hence OPT = max D. The cells are filled row
by row, so `D[i-1][j-1]` is always computed before `D[i][j]`.

Handout example: "CAABC" (from position 3 of X) against "CAACC" (from position 4 of Y)
gives w<sub>C</sub> + w<sub>A</sub> + w<sub>A</sub> − δ + w<sub>C</sub> = 2w<sub>C</sub> + 2w<sub>A</sub> − δ,
since only the fourth position (B against C) is a mismatch. This matches the score
defined in 1.1. A numerical check with the actual weights is given in Section 1.8 (Experimental setup).


### 1.5 Extracting the optimal pair

While filling the table we keep the cell `(i*, j*)` with the largest value; with ties, the first one reached in row order is kept. If the largest value is 0, the answer is the empty pair (L = 0).

If `D[i][j] > 0`, the maximum in the recurrence came from the second term, so `D[i][j] = D[i-1][j-1] + s(x_i, y_j)`. Starting from `(i*, j*)` we step to `(i*-1, j*-1)`, `(i*-2, j*-2)`, … while D > 0, and count the visited cells as L.
The walk stops at the cell `(i*-L, j*-L)`, which has D = 0; this happens at row or
column 0 at the latest. The pair starts right after that cell:
`X[i*-L+1 .. i*]` and `Y[j*-L+1 .. j*]`. Adding the L equations
`D[i][j] = D[i-1][j-1] + s(x_i, y_j)` along the walk, and using D = 0 at the stopping
cell, shows that the score of this pair is exactly `D[i*][j*]`.

**Optional space-saving variant (not required by the assignment).** The full table above is enough for both the DP and the extraction. 
To reduce memory, I also implemented a two-row version: with each value we also store the length of its pair, `R[i][j] = R[i-1][j-1] + 1` if `D[i][j] > 0` and 0 otherwise. 

R counts the consecutive positive cells on the diagonal, so `R[i*][j*]` equals the L found by the walk above. Only two rows of D and R are kept, plus the best value, `(i*, j*)` and its R. Since the tie rule is also the same, both versions return the same pair. This was confirmed on every test where both were run (`full_eq_lin = yes`).


### 1.6 Complexity

**Time: Θ(nm).** The running time is O(nm); in fact it is Θ(nm). The DP table has (n+1)(m+1) cells. Row 0 and column 0 are set to 0. Each of
the other nm cells needs one addition, one comparison with 0, and one comparison with the best
value so far (to record `(i*, j*)`), so O(1) per cell. All nm cells are computed for every
input, so filling the table takes Θ(nm). The traceback takes at most min(n, m) steps, and
reading the input takes Θ(n + m). The total is Θ(nm) for n, m ≥ 1.

**Space.**

| Version | Arrays | Auxiliary space |
|---|---|---|
| Full table | (n+1)(m+1) four-byte integers | Θ(nm) |
| Two-row (optional) | 2 rows of D + 2 rows of R, each m+1 four-byte integers | Θ(m) |

Scores are stored as `int`, which is 4 bytes on the test platform (x86-64, g++); the measured heap peak (4000+1)² · 4 = 64,032,004 bytes confirms this.


Y (length m) is the string along the columns. The two-row auxiliary space can be reduced to
Θ(min(n, m)) by placing the shorter string along the columns. Both versions also store the
input strings, n + m bytes, i.e. Θ(n + m).

**Overflow.** Every `D[i][j]` is nonnegative and at most w<sub>max</sub> · min(n, m), where
w<sub>max</sub> = max<sub>l</sub> w<sub>l</sub>. For the longest input used in the experiments (8000),
this is 12702 · 8000 = 101,616,000 < 2<sup>31</sup> − 1, so 32-bit signed integers are enough for all
experiments. The program checks this bound before each run.

### 1.7 Validation
`bruteForce` tries every start pair and every length in O(nm · min(n, m))
time and serves as an independent reference. `selfTest` compares it with both DP versions
on 500 random cases (fixed seed 7), with n and m drawn independently from 1..60. The cases
are split evenly (125 each) between Scenario 1 (δ = 10) and Scenario 2 (δ chosen at random
from δ<sub>0</sub>, …, δ<sub>11</sub>), and between English-frequency letters and a 3-letter alphabet.
The 3-letter alphabet makes matches frequent, which produces long pairs containing
mismatches and many tied optima, the cases where the recurrence or the tie rule is most
likely to go wrong. A case passes only if
(1) the scores of `solveFull` and `solveLinear` both equal the brute-force optimum,
(2) the two DP versions return the same start positions and length, and
(3) the score of that pair, recomputed directly from the strings, equals the reported score.
(1) and (3) together show that the returned pair is optimal. Positions are not compared with
the brute force, because it breaks ties by the earliest start while the DP breaks ties by the
earliest end cell, so the two may return different pairs that are equally optimal.
All 500 cases passed. In the main experiments, `score_ok` (recomputed score) and
`full_eq_lin` (same pair from both DP versions) were "yes" on every row where they apply;
for 8000 × 8000 only the two-row version was run, so `full_eq_lin` is n/a there.


### 1.8 Experimental setup

Problem 1: positions are 1-based, as in the handout.

**Scenario 1:** w<sub>l</sub> = 1 for all letters, δ = 10. One mismatch cancels ten matches, so including a mismatch only helps if it brings in more than ten extra matches.

**Scenario 2:** w<sub>l</sub> is proportional to the frequency of letter l in English. I used the standard table by Lewand (*Cryptological Mathematics*), in units of 0.001 %:

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

**Check on the handout strings.** For X = "ABCAABCAA" and Y = "ABBCAACCBBBBBB":

- Scenario 2, δ = 5814 (t = 5): the handout pair "CAABC" / "CAACC" scores
  2·2782 + 2·8167 − 5814 = 16084, as printed by the program. This pair is not optimal.
  The optimum is "AABCAA" (position 4 of X) / "ABBCAA" (position 1 of Y), with one mismatch
  and score 3·8167 + 1492 + 2782 − 5814 = 22961.
- Scenario 1 (w = 1, δ = 10): one mismatch costs more than any common substring of length
  up to 9 can gain, so the optimum has no mismatches. It is "BCAA" (position 2 of X, position 3
  of Y) with score 4. This optimum is not unique, since "BCAA" also starts at position 6 of X.

**Data:** random strings with letters drawn either uniformly from the 26 letters or from the English distribution above. Sizes were n = m ∈ {100, 500, 1000, 2000, 4000, 8000}, plus one unequal case, 1000 × 4000. The seed is fixed (20261103). Each configuration uses one random instance; reported times are averages over repeated runs on that instance. The δ sweep reuses the same pair of 2000 × 2000 strings for every δ. The full table was run up to 4000 × 4000 (64 MB); for 8000 × 8000 only the two-row version was run.

**Measurement:** time with `std::chrono::steady_clock`, using one warm-up run and then repetitions until at least 200 ms and 3 runs had passed. The reported number is the average per run. Memory was measured in a separate process for each configuration, in two ways:

- **heap peak:** the peak number of bytes requested through the overridden operator new; it excludes allocator metadata and rounding. In Problem 1 the inputs are created before the measurement starts, so theory − heap peak = n + m (the input bytes). In Problem 2 the matrix B is included in both columns, so they coincide;

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

**Time.** From 100 × 100 to 8000 × 8000 (a 6400-fold increase in n·m), the time per cell
stays between 2.1 and 2.4 ns (Table 1.1), and the log–log plot has slope close to 1,
consistent with Θ(nm). The 2000 × 2000 and 1000 × 4000 inputs have the same number of cells
and similar times, so the cost depends on n·m, not on the two lengths separately. In these size experiments the two versions ran at similar speed. The largest gap was at
4000 × 4000: 43.4 ms for the full table versus 35.0 ms for the two-row version (24 % more
relative to the two-row time, 19.5 % less relative to the full-table time).


**Memory.** The heap peak equals the size of the DP arrays exactly: (4000+1)<sup>2</sup> · 4 =
64,032,004 bytes for the full table versus 4 · 4001 · 4 = 64,016 bytes for the two-row
version (two rows each of D and R), a factor of 1000. For the full table, ΔRSS grows by
46,736 KB from 2000 × 2000 to 4000 × 4000, while the array grows by 46,891 KB (99.7 %).
On top of the array there is a roughly constant extra 1269–1423 KB whose source I did not
identify. The two-row arrays do not show up in ΔRSS: `ru_maxrss` is the highest RSS since
start-up, already 4244 KB, and these arrays never push it higher. For them the heap peak
is the meaningful measurement.

**Scenario 1.** ("Uniform" strings draw each letter with probability 1/26; "English" strings
use the English letter frequencies.) Every optimum had 0 mismatches, and this is forced here. Consider the prefix of an optimal pair that ends at, and includes, its first mismatch. If r matches precede that mismatch, this prefix scores r − 10. It cannot be negative, since deleting it would strictly raise the score, so r ≥ 10. But then those r matches form an exact common substring of score r ≥ 10, while every optimum scored at most 6 (Table 1.4). So no optimum contains a mismatch, and the answer is the longest exact common substring.

Let q be the probability that two random letters are equal: q = 1/26 ≈ 0.038 (uniform) and
q = Σ p<sub>l</sub><sup>2</sup> ≈ 0.0655 (English). The rough estimate log<sub>1/q</sub>(nm) gives 4.2 and
5.1 at 1000 × 1000, and 5.5 and 6.6 at 8000 × 8000. The observed lengths are at most 1.3 below
these values and grow logarithmically in the same way. English has the larger q, so its
matches are longer.

**Scenario 2.** The weights are the English-frequency weights for both string distributions.
The effect of δ can be explained by the expected score of one step along a diagonal of two
random strings:

```
E[s] = sum_l p_l^2 * w_l  -  delta * (1 - sum_l p_l^2)
```

If E[s] > 0, the running sum tends to grow, so long segments are favored. If E[s] < 0, long
segments tend to accumulate negative score, so shorter segments are generally favored.
E[s] = 0 gives δ<sub>c</sub> ≈ 154 (uniform) and δ<sub>c</sub> ≈ 573 (English). At δ = 74, below both
thresholds, the optimum covered 1852 (English) and 1635 (uniform) of the 2000 positions.
All ten intermediate values are above both thresholds and gave much shorter pairs. This is
a heuristic explanation, not a guarantee about the optimal length.

The number of mismatches in the optimum cannot increase with δ. Let A be optimal for
δ<sub>1</sub> and B optimal for δ<sub>2</sub> > δ<sub>1</sub>, with match weights W<sub>A</sub>, W<sub>B</sub> and
mismatch counts k<sub>A</sub>, k<sub>B</sub>. Then

```
W_A - k_A * delta1  >=  W_B - k_B * delta1
W_B - k_B * delta2  >=  W_A - k_A * delta2
```

Adding the two gives (k<sub>A</sub> − k<sub>B</sub>)(δ<sub>2</sub> − δ<sub>1</sub>) ≥ 0, so k<sub>A</sub> ≥ k<sub>B</sub>.
The data agree: 1702 → 57 → 10 → 7 → 1 → 0 (English) and 1543 → 4 → 0 (uniform).
The length of the optimum also decreased, but this argument does not prove that.

While the optimal pair stays the same, its score changes by −(mismatches) × Δδ: it drops by
exactly 1148 per step for English at t = 4 → 5 → 6 (1 mismatch), and by 4592 = 4 · 1148 for
uniform at t = 1 → 2 (4 mismatches). Once the optimum has no mismatches, its score no longer
depends on δ while every other pair's score can only fall, so it stays optimal for all larger
δ: "EEB" with score 2 · 12702 + 1492 = 26896 for uniform (t ≥ 3), and "ESOTE" with score
48294 for English (t ≥ 7), a pair made of frequent, high-weight letters.

δ also changed the running time of the two-row solver: 22.6 ms at δ = 1222 versus 8.7 ms at
δ = 12702, while the full-table solver stayed between 8.0 and 9.2 ms. The cause may be
branch misprediction, but I did not verify this. The cost is Θ(nm) in every case.




## Problem 2: Largest All-Zero Square Sub-matrix

### 2.1 Problem

Given an m × n Boolean matrix B (non-zero entries stored as 1), find the largest k such that
B contains a k × k block of consecutive rows and columns that is all 0, and return its
position (top-left corner, 0-indexed). If B has no zero entry, the answer is k = 0.

### 2.2 Optimization function

For 0 ≤ i < m and 0 ≤ j < n, `D[i][j]` is the side length of the largest all-zero square whose
**bottom-right corner** is (i, j), and 0 if B[i][j] = 1.

### 2.3 Recurrence (Bellman equations)

```
D[i][j] = 0                                               if B[i][j] = 1
        = 1                                               if B[i][j] = 0 and (i = 0 or j = 0)
        = 1 + min( D[i-1][j], D[i][j-1], D[i-1][j-1] )    otherwise
k*      = max over all i, j of D[i][j]
```
k* = optimal solution k

The code adds an extra row and column of D values equal to 0 in front of the table, so the boundary case is handled by the general formula (1 + min(0, …) = 1).

### 2.4 Correctness

We prove the recurrence by induction in row-major order; the three neighbors of (i, j) come
earlier in this order, so their D values are correct by the induction hypothesis.

If B[i][j] = 1, no all-zero square can end at (i, j), so D = 0. In the boundary case with
B[i][j] = 0, the only square ending at (i, j) is the cell itself, so D = 1.

Now let B[i][j] = 0 with i, j ≥ 1, and let t = min(D[i−1][j], D[i][j−1], D[i−1][j−1]).
We show D[i][j] = t + 1.

*D[i][j] ≤ t + 1.* Take an all-zero k × k square ending at (i, j) with k ≥ 2. Removing one row
(top or bottom) and one column (left or right) gives three all-zero (k−1) × (k−1) squares,
ending at (i−1, j), (i, j−1) and (i−1, j−1). So each neighbor value is at least k − 1, which
gives t ≥ k − 1, i.e. k ≤ t + 1. For k = 1 the bound holds because t ≥ 0.

*D[i][j] ≥ t + 1.*  If t = 0, the zero cell B[i][j] itself is a 1 × 1 all-zero square, so the claim holds. Now suppose t ≥ 1. Each neighbor has an all-zero square of side at least t ending at it, so it also has one of side exactly t by taking the t × t sub-square with the same bottom-right corner. Since D[i−1][j−1] ≥ t, we have i − t ≥ 0 and j − t ≥ 0, so the (t+1) × (t+1) square with rows i−t..i and columns j−t..j lies inside B.

A cell (r, c) of this square is

- (r, c) = (i, j): zero by assumption;
- r = i, c < j: covered by the t × t square ending at (i, j−1) (rows i−t+1..i, columns j−t..j−1);
- r < i, c = j: covered by the square ending at (i−1, j) (rows i−t..i−1, columns j−t+1..j);
- r < i, c < j: covered by the square ending at (i−1, j−1) (rows i−t..i−1, columns j−t..j−1).

So every cell is 0, and a (t+1) × (t+1) all-zero square ends at (i, j).

Each `D[i][j]` is the side of an actual all-zero square, so max D is at most the answer.
Conversely, every all-zero square ends at some cell (i, j), so its side is at most `D[i][j]`.
Hence k* = max D.

### 2.5 Extracting the optimal square

While filling the table row by row (only the current and previous rows of D are kept), we keep
`k*` and the first cell `(i*, j*)` in row-major order where it occurs (strict comparison, so
ties keep the earlier cell). The square is rows `i*-k*+1 .. i*`, columns `j*-k*+1 .. j*`. It is all zero by the
second half of 2.4, which uses exactly these rows and columns. No walk back is needed, because
the corner and the side length determine the square. If k* = 0, there is no all-zero square.

### 2.6 Complexity

**Time: Θ(mn).** The running time is O(nm); in fact it is Θ(nm). Each cell needs one test of B[i][j], one minimum of three values, and one
comparison with the best value so far, so O(1) per cell. All mn cells are computed for every
input, so the total is Θ(mn). Extraction is O(1), since the corner and side are already known.

**Space.** Following the hint, B is stored with one byte per entry (`uint8_t`). The DP values
need more than a byte, because k can be as large as min(m, n), which is 1000 > 255 in the
largest experiment. Each row of D depends only on the row above, so only two rows of
four-byte integers are kept.

| Part | Size | Space |
|---|---|---|
| Matrix B | mn bytes | Θ(mn) |
| Two DP rows | 2(n+1) · 4 bytes | Θ(n) extra |
| Total | mn + 8(n+1) bytes | Θ(mn) |

The extra space could be reduced to Θ(min(m, n)) by running the DP along the shorter
dimension; the implementation does not do this. The total stays Θ(mn) because of B.

### 2.7 Validation

`verify` checks every result independently with 2D prefix sums. It confirms that (1) the
reported square lies inside B and contains only zeros, and (2) no all-zero (k*+1) × (k*+1)
square exists in B. Any larger all-zero square would contain one of size (k*+1) × (k*+1), so
(1) and (2) together prove that k* is optimal. All experiments passed. A hand-made 5 × 6 matrix
with known answer k = 3 at top-left (0, 1) (0-indexed) also passed. `verify` runs after the
measurements, so its memory is not included in them.

### 2.8 Experimental setup

Row and column indices in Problem 2 are 0-based.

The sizes were 10 × 10, 10 × 100, 10 × 1000, 100 × 1000 and 1000 × 1000. Each entry is set to 1
independently with probability p ∈ {0.5, 0.1, 0.01} (fixed seed 20261103). Several values of p
are used because at p = 0.5 the answer is only 2–4. Timing is the same as in Problem 1, except
that at least 5 runs are averaged. Memory was measured in a separate process for each configuration.

### 2.9 Results

**Table 2.1** Running time in µs, with k* in parentheses (from `results_p2.csv`).

| m × n | cells | p = 0.5 | p = 0.1 | p = 0.01 |
|---|---|---|---|---|
| 10 × 10 | 10<sup>2</sup> | 0.196 (2) | 0.272 (3) | 0.309 (9) |
| 10 × 100 | 10<sup>3</sup> | 1.209 (2) | 2.524 (8) | 3.413 (10) |
| 10 × 1000 | 10<sup>4</sup> | 43.712 (3) | 26.194 (8) | 34.236 (10) |
| 100 × 1000 | 10<sup>5</sup> | 511.615 (4) | 292.002 (9) | 342.590 (33) |
| 1000 × 1000 | 10<sup>6</sup> | 5118.774 (4) | 2909.975 (10) | 3418.036 (30) |

Table 2.1 comes from a single run of all configurations (`results_p2.csv`). The memory runs
(`mem_p2.csv`) start a separate process for each configuration, which generates different
random matrices, so their k* values can differ; the memory figures do not depend on k*.


**Table 2.2** Memory. Theory and heap peak do not depend on p. (from `mem_p2.csv`)

| m × n | theory (B) | heap peak (B) | ΔRSS (KB) |
|---|---|---|---|
| 10 × 10 | 188 | 188 | 0 |
| 10 × 100 | 1,808 | 1,808 | 0 |
| 10 × 1000 | 18,008 | 18,008 | 0 |
| 100 × 1000 | 108,008 | 108,008 | 0 |
| 1000 × 1000 | 1,008,008 | 1,008,008 | 2,292 – 2,416 (one run per p) |

![Figure 2](p2_plots.png)

**Figure 2.** Left: running time vs. m·n for three values of p (log–log). Right: memory vs. m·n
for p = 0.1; theory and heap peak coincide. ΔRSS values of 0 are plotted at 1 KB because 0
cannot be shown on a log axis.

### 2.10 Analysis

**Time.** For p = 0.1 and p = 0.01 the time per cell stays within 2.5–2.9 ns and 3.1–3.4 ns
from 10<sup>2</sup> to 10<sup>6</sup> cells, and the log–log lines have slope close to 1, consistent
with Θ(mn). For p = 0.5 the time per cell is 1.2–2.0 ns on the two smallest matrices,
4.4 ns at 10<sup>4</sup> cells and 5.1 ns from 10<sup>5</sup> cells on, so its line is steeper in
between and the lines cross. The two smallest runs take under 4 µs and process the same tiny
matrix many times in a row, so I do not treat them as representative.

The matrix content changes only the constant factor; k* has no effect, since every cell is
computed. On the two largest matrices one cell costs 5.1 ns (p = 0.5), 2.9 ns (p = 0.1) and
3.4 ns (p = 0.01). A possible explanation, which I did not measure: at p = 0.5 the branch on
B[i][j] is unpredictable, and at p = 0.01 more cells are zero (99 % instead of 90 %) and need
the minimum of three neighbors.

**Memory.** The heap peak equals mn + 8(n+1) bytes exactly for every size. ΔRSS is 0 up to
10<sup>5</sup> cells (heap peak at most 106 KB). Possible causes are that these allocations used
pages the process already had resident, or that the kernel updates its RSS counter in
batches; I did not confirm either. At 10<sup>6</sup> cells, ΔRSS exceeds the heap peak (984 KB)
by 1308–1432 KB, similar to the extra 1269–1423 KB seen in Problem 1.

**Size of k\*.** Let N<sub>k</sub> be the number of all-zero k × k squares. Then
E[N<sub>k</sub>] = (m−k+1)(n−k+1)(1−p)<sup>k²</sup>, and since k* ≥ k exactly when N<sub>k</sub> ≥ 1,
P(k* ≥ k) ≤ E[N<sub>k</sub>]. Setting E[N<sub>k</sub>] = 1 and approximating (m−k+1)(n−k+1) by mn gives

```
k ≈ sqrt( ln(mn) / -ln(1-p) )
```

Above this value E[N<sub>k</sub>] < 1 and keeps decreasing. Thus the estimate gives a first-moment threshold for the typical scale of k*, rather than a deterministic upper bound. For 1000 × 1000, the threshold estimate gives 4.5, 11.5 and 37.1 for p = 0.5, 0.1 and 0.01; the observed values are 4, 10 and 30, respectively. (When m = 10, k* is capped by min(m, n) = 10; this cap is reached at p = 0.01 for the 10 × 100 and 10 × 1000 instances, while the 10 × 10 instance has k* = 9.) The estimate grows only like sqrt(ln(mn)), so k* increases very
slowly with the matrix size but strongly as p decreases.

---

## Conclusion

Both problems are solved by a table DP in which each cell takes constant time, so the running
time is Θ(nm) for Problem 1 and Θ(mn) for Problem 2. Each recurrence reads only the current
and previous rows, so the extra memory can be reduced to one row's worth; in Problem 1,
storing the pair length in each cell makes the answer extractable without the full table.
The experiments agree with the analysis. The time per cell stays nearly constant as the input
grows (except for the tiny p = 0.5 matrices), and the measured heap memory equals the computed
array sizes. All checks passed: in Problem 1, both DP versions matched a brute force on 500
random small cases, and in every experiment the recomputed score matched and both versions
returned the same pair where both were run; in Problem 2, every result passed the prefix-sum
verification.
