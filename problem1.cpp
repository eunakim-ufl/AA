// problem2.cpp  —  Largest all-zero square sub-matrix (iterative DP)
// Build:  g++ -O2 -std=c++17 problem2.cpp -o p2
// Run all experiments:   ./p2 > results_p2.csv
// Single run (memory):   ./p2 <m> <n> <p>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

// ===== memory measurement =====
#include <new>
#include <fstream>
#include <unistd.h>
static std::size_t g_heapCur = 0, g_heapPeak = 0;
// Count every heap allocation (std::vector uses this)
void* operator new(std::size_t sz) {
    void* p = std::malloc(sz + 16);
    if (!p) throw std::bad_alloc();
    *static_cast<std::size_t*>(p) = sz;
    g_heapCur += sz;
    if (g_heapCur > g_heapPeak) g_heapPeak = g_heapCur;
    return static_cast<char*>(p) + 16;
}
void operator delete(void* p) noexcept {
    if (!p) return;
    char* b = static_cast<char*>(p) - 16;
    g_heapCur -= *reinterpret_cast<std::size_t*>(b);
    std::free(b);
}
void operator delete(void* p, std::size_t) noexcept { operator delete(p); }
void resetHeapPeak() { g_heapPeak = g_heapCur; }
long currentRssKB() {   // Linux: resident memory RIGHT NOW (not the historical peak)
    std::ifstream f("/proc/self/statm");
    long size = 0, res = 0;
    if (!(f >> size >> res)) return -1;
    return res * (sysconf(_SC_PAGESIZE) / 1024);
}

using Matrix = std::vector<uint8_t>;          // row-major, B[i*n + j], 1 byte per cell

struct Result { int k; int row; int col; };    // side length + top-left corner

// ---------- synthetic data: each cell is 1 with probability p ----------
Matrix generate(int m, int n, double p, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::bernoulli_distribution one(p);
    Matrix B((size_t)m * n);
    for (size_t t = 0; t < B.size(); ++t) B[t] = one(rng) ? 1 : 0;
    return B;
}

// ---------- DP: D[i][j] = side of largest zero square with bottom-right (i,j) ----------
// Only two rows of D are kept (O(n) extra space). Column 0 is a sentinel (= 0),
// and the initial 'prev' row of zeros is the sentinel row above row 0.
Result largestZeroSquare(const Matrix& B, int m, int n) {
    std::vector<int> prev(n + 1, 0), cur(n + 1, 0);
    int best = 0, bi = -1, bj = -1;
    for (int i = 0; i < m; ++i) {
        cur[0] = 0;
        const uint8_t* row = &B[(size_t)i * n];
        for (int j = 1; j <= n; ++j) {
            if (row[j - 1]) {
                cur[j] = 0;
            } else {
                cur[j] = 1 + std::min(prev[j], std::min(cur[j - 1], prev[j - 1]));
                if (cur[j] > best) { best = cur[j]; bi = i; bj = j - 1; }
            }
        }
        std::swap(prev, cur);                 // O(1) pointer swap
    }
    if (best == 0) return {0, -1, -1};
    // Solution extraction: bottom-right (bi,bj) and side best -> top-left corner
    return {best, bi - best + 1, bj - best + 1};
}

// ---------- independent verifier using 2D prefix sums ----------
// (1) the reported square is inside B and all zeros, (2) no (k+1)x(k+1) zero square exists.
bool verify(const Matrix& B, int m, int n, const Result& r) {
    std::vector<int> S((size_t)(m + 1) * (n + 1), 0);
    auto at = [&](int i, int j) -> int& { return S[(size_t)i * (n + 1) + j]; };
    for (int i = 1; i <= m; ++i)
        for (int j = 1; j <= n; ++j)
            at(i, j) = B[(size_t)(i - 1) * n + (j - 1)] + at(i - 1, j) + at(i, j - 1) - at(i - 1, j - 1);
    auto sum = [&](int r0, int c0, int k) {
        return at(r0 + k, c0 + k) - at(r0, c0 + k) - at(r0 + k, c0) + at(r0, c0);
    };
    if (r.k > 0) {
        if (r.row < 0 || r.col < 0 || r.row + r.k > m || r.col + r.k > n) return false;
        if (sum(r.row, r.col, r.k) != 0) return false;
    }
    int K = r.k + 1;
    if (K <= std::min(m, n))
        for (int r0 = 0; r0 + K <= m; ++r0)
            for (int c0 = 0; c0 + K <= n; ++c0)
                if (sum(r0, c0, K) == 0) return false;
    return true;
}

// ---------- measurement helpers ----------
long peakRssKB() {
#if defined(__APPLE__)
    struct rusage ru; getrusage(RUSAGE_SELF, &ru); return ru.ru_maxrss / 1024;  // bytes on macOS
#elif defined(__unix__)
    struct rusage ru; getrusage(RUSAGE_SELF, &ru); return ru.ru_maxrss;         // KB on Linux
#else
    return -1;
#endif
}

// Repeat until >= 200 ms total and >= 5 reps, return average microseconds per run
double timeIt(const Matrix& B, int m, int n, Result& out, int& reps) {
    using clk = std::chrono::steady_clock;
    out = largestZeroSquare(B, m, n);          // warm-up (caches)
    reps = 0;
    double total = 0.0;
    auto start = clk::now();
    do {
        out = largestZeroSquare(B, m, n);
        ++reps;
        total = std::chrono::duration<double, std::micro>(clk::now() - start).count();
    } while (total < 200000.0 || reps < 5);
    return total / reps;
}

void runOne(int m, int n, double p, uint64_t seed) {
    long rssBefore = currentRssKB();
    std::size_t heapBase = g_heapCur; resetHeapPeak();
    Matrix B = generate(m, n, p, seed);
    Result r; int reps;
    double us = timeIt(B, m, n, r, reps);
    std::size_t heapPeak = g_heapPeak - heapBase;      // B + two DP rows
    long rssAfter = currentRssKB();
    size_t theory = (size_t)m * n * sizeof(uint8_t) + 2 * (size_t)(n + 1) * sizeof(int);
    bool ok = verify(B, m, n, r);
    std::printf("%d,%d,%.3f,%d,%d,%d,%.3f,%d,%zu,%zu,%ld,%ld,%s\n",
                m, n, p, r.k, r.row, r.col, us, reps, theory, heapPeak,
                rssBefore, rssAfter, ok ? "yes" : "NO");
}

// ---------- hand-made test with known answer (k=3, top-left (0,1)) ----------
void demo() {
    const int m = 5, n = 6;
    Matrix B = { 1,0,0,0,1,0,
                 0,0,0,0,1,0,
                 1,0,0,0,0,1,
                 0,1,0,0,0,0,
                 0,0,1,0,0,1 };
    Result r = largestZeroSquare(B, m, n);
    std::fprintf(stderr, "[demo] k=%d top-left=(%d,%d) verified=%s\n",
                 r.k, r.row, r.col, verify(B, m, n, r) ? "yes" : "NO");
    for (int i = r.row; i < r.row + r.k; ++i) {
        std::fprintf(stderr, "       ");
        for (int j = r.col; j < r.col + r.k; ++j) std::fprintf(stderr, "%d ", B[i * n + j]);
        std::fprintf(stderr, "\n");
    }
}

int main(int argc, char** argv) {
    if (argc == 4) {                           // single run -> clean per-process memory
        runOne(std::atoi(argv[1]), std::atoi(argv[2]), std::atof(argv[3]), 20261103ULL);
        return 0;
    }
    demo();                                    // goes to stderr, CSV stays clean
    std::printf("m,n,p,k,top_row,top_col,avg_time_us,reps,theory_bytes,heap_peak_bytes,"
                "rss_before_kb,rss_after_kb,verified\n");
    const int sizes[5][2] = {{10,10},{10,100},{10,1000},{100,1000},{1000,1000}};
    const double dens[3] = {0.5, 0.1, 0.01};
    uint64_t seed = 20261103ULL;
    for (auto& s : sizes)
        for (double p : dens)
            runOne(s[0], s[1], p, seed++);
    return 0;
}
