// problem1.cpp — Weighted approximate common substring (no gaps), iterative DP
// Build:  g++ -O2 -std=c++17 problem1.cpp -o p1
// Run all experiments:  ./p1 > results_p1.csv     (demo + self-test go to stderr)
// Memory (one process per run):  ./p1 mem <n> <m> <full|lin>
#include <array>
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>
#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

// ===== memory measurement =====
#include <new>
#include <fstream>
#include <unistd.h>
static std::size_t g_heapCur = 0, g_heapPeak = 0;
// Count every heap allocation (std::vector, std::string use this)
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



using Weights = std::array<int, 26>;

// Lewand English letter frequencies in units of 0.001 %  (E = 12.702 % -> 12702).
// Integers keep w_l exactly proportional to frequency and all scores exact.
const int FREQ[26] = {
    8167, 1492, 2782, 4253, 12702, 2228, 2015, 6094, 6966, 153, 772, 4025, 2406,  // A-M
    6749, 7507, 1929,   95,  5987, 6327, 9056, 2758,  978, 2360, 150, 1974,   74  // N-Z
};

struct Result {
    int score;   // optimal score (0 = empty substring)
    int s1, s2;  // 0-indexed start positions in X and Y
    int len;     // common length of both substrings
};

inline int sc(char a, char b, const Weights& w, int delta) {
    return a == b ? w[a - 'A'] : -delta;
}

// ---------------- DP with full table + traceback extraction ----------------
// D[i][j] = best score of a pair of equal-length substrings ENDING at X[i-1], Y[j-1]
// D[i][j] = max(0, D[i-1][j-1] + s(X[i-1], Y[j-1])),  D[0][*] = D[*][0] = 0
Result solveFull(const std::string& X, const std::string& Y, const Weights& w, int delta) {
    const int n = (int)X.size(), m = (int)Y.size();
    std::vector<int> D((size_t)(n + 1) * (m + 1), 0);
    auto at = [&](int i, int j) -> int& { return D[(size_t)i * (m + 1) + j]; };
    int best = 0, bi = 0, bj = 0;
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= m; ++j) {
            int v = at(i - 1, j - 1) + sc(X[i - 1], Y[j - 1], w, delta);
            if (v < 0) v = 0;
            at(i, j) = v;
            if (v > best) { best = v; bi = i; bj = j; }
        }
    // Extraction: walk back along the diagonal while D > 0 (iterative, no recursion)
    int len = 0, i = bi, j = bj;
    while (i > 0 && j > 0 && at(i, j) > 0) { ++len; --i; --j; }
    return {best, bi - len, bj - len, len};
}

// ---------------- Space-optimized DP: two rows + run length ----------------
// L[j] = length of the run ending at (i,j); allows extraction without the full table.
Result solveLinear(const std::string& X, const std::string& Y, const Weights& w, int delta) {
    const int n = (int)X.size(), m = (int)Y.size();
    std::vector<int> pD(m + 1, 0), cD(m + 1, 0), pL(m + 1, 0), cL(m + 1, 0);
    int best = 0, bi = 0, bj = 0, blen = 0;
    for (int i = 1; i <= n; ++i) {
        cD[0] = 0; cL[0] = 0;
        for (int j = 1; j <= m; ++j) {
            int v = pD[j - 1] + sc(X[i - 1], Y[j - 1], w, delta);
            if (v > 0) { cD[j] = v; cL[j] = pL[j - 1] + 1; }
            else       { cD[j] = 0; cL[j] = 0; }
            if (v > best) { best = v; bi = i; bj = j; blen = cL[j]; }
        }
        std::swap(pD, cD); std::swap(pL, cL);
    }
    return {best, bi - blen, bj - blen, blen};
}

// ---------------- Brute force O(n*m*min(n,m)) for verification ----------------
Result bruteForce(const std::string& X, const std::string& Y, const Weights& w, int delta) {
    Result r{0, 0, 0, 0};
    const int n = (int)X.size(), m = (int)Y.size();
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < m; ++b) {
            int s = 0;
            for (int L = 1; a + L <= n && b + L <= m; ++L) {
                s += sc(X[a + L - 1], Y[b + L - 1], w, delta);
                if (s > r.score) r = {s, a, b, L};
            }
        }
    return r;
}

int scoreOf(const std::string& X, const std::string& Y, int s1, int s2, int len,
            const Weights& w, int delta) {
    int s = 0;
    for (int k = 0; k < len; ++k) s += sc(X[s1 + k], Y[s2 + k], w, delta);
    return s;
}
int mismatches(const std::string& X, const std::string& Y, int s1, int s2, int len) {
    int c = 0;
    for (int k = 0; k < len; ++k) c += (X[s1 + k] != Y[s2 + k]);
    return c;
}

// ---------------- Synthetic data ----------------
std::string genAlpha(int n, int k, std::mt19937_64& rng) {      // uniform over first k letters
    std::uniform_int_distribution<int> d(0, k - 1);
    std::string s(n, 'A');
    for (auto& c : s) c = char('A' + d(rng));
    return s;
}
std::string genEnglish(int n, std::mt19937_64& rng) {           // English letter distribution
    std::discrete_distribution<int> d(FREQ, FREQ + 26);
    std::string s(n, 'A');
    for (auto& c : s) c = char('A' + d(rng));
    return s;
}

// ---------------- Scenarios ----------------
Weights scenario1() { Weights w; w.fill(1); return w; }          // delta = 10
Weights scenario2() { Weights w; for (int l = 0; l < 26; ++l) w[l] = FREQ[l]; return w; }
int deltaS2(int t) {  // t = 0..11 ; t = 1..10 are the 10 intermediate values
    const int wmin = *std::min_element(FREQ, FREQ + 26), wmax = *std::max_element(FREQ, FREQ + 26);
    return (int)std::lround(wmin + t * (wmax - wmin) / 11.0);
}

// ---------------- Measurement ----------------
volatile int g_sink;
template <class F>
double timeIt(F f) {                       // avg microseconds; >= 200 ms total and >= 3 reps
    using clk = std::chrono::steady_clock;
    f();                                   // warm-up
    int reps = 0; double total = 0;
    auto start = clk::now();
    do { f(); ++reps; total = std::chrono::duration<double, std::micro>(clk::now() - start).count(); }
    while (total < 200000.0 || reps < 3);
    return total / reps;
}
long peakRssKB() {
#if defined(__APPLE__)
    struct rusage ru; getrusage(RUSAGE_SELF, &ru); return ru.ru_maxrss / 1024;
#elif defined(__unix__)
    struct rusage ru; getrusage(RUSAGE_SELF, &ru); return ru.ru_maxrss;
#else
    return -1;
#endif
}
size_t theoryFull(int n, int m) { return (size_t)(n + 1) * (m + 1) * sizeof(int) + n + m; }
size_t theoryLin(int n, int m)  { return 4 * (size_t)(m + 1) * sizeof(int) + n + m; }
std::string show(const std::string& s, int st, int len) {
    if (len == 0) return "-";
    return len <= 30 ? s.substr(st, len) : s.substr(st, 30) + "...";
}

void runCase(const char* exp, const char* dist, const std::string& X, const std::string& Y,
             const Weights& w, int delta, int t) {
    const int n = (int)X.size(), m = (int)Y.size();
    if ((long long)*std::max_element(w.begin(), w.end()) * std::min(n, m) > INT_MAX) {
        std::fprintf(stderr, "score could overflow int\n"); std::exit(1);
    }
    const bool doFull = (size_t)n * m <= 16000000ULL;            // full table <= 64 MB
    Result r = solveLinear(X, Y, w, delta);
    double tL = timeIt([&] { g_sink = solveLinear(X, Y, w, delta).score; });
    double tF = -1; bool same = true;
    if (doFull) {
        Result f = solveFull(X, Y, w, delta);
        same = f.score == r.score && f.s1 == r.s1 && f.s2 == r.s2 && f.len == r.len;
        tF = timeIt([&] { g_sink = solveFull(X, Y, w, delta).score; });
    }
    bool ok = scoreOf(X, Y, r.s1, r.s2, r.len, w, delta) == r.score;
    std::printf("%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.2f,%.2f,%zu,%zu,%s,%s,%s,%s\n",
                exp, dist, n, m, delta, t, r.score, r.len,
                r.len ? r.s1 + 1 : 0, r.len ? r.s2 + 1 : 0,           // 1-indexed like the handout
                mismatches(X, Y, r.s1, r.s2, r.len), tF, tL,
                doFull ? theoryFull(n, m) : 0, theoryLin(n, m),
                doFull ? (same ? "yes" : "NO") : "n/a", ok ? "yes" : "NO",
                show(X, r.s1, r.len).c_str(), show(Y, r.s2, r.len).c_str());
}

// ---------------- Correctness checks ----------------
void demo() {
    const std::string X = "ABCAABCAA", Y = "ABBCAACCBBBBBB";
    Weights w2 = scenario2(); int d = deltaS2(5);
    std::fprintf(stderr, "[demo] handout alignment %s / %s : 2*w_C + 2*w_A - delta = 2*%d + 2*%d - %d = %d\n",
                 X.substr(2, 5).c_str(), Y.substr(3, 5).c_str(), w2['C' - 'A'], w2['A' - 'A'], d,
                 scoreOf(X, Y, 2, 3, 5, w2, d));
    Result a = solveFull(X, Y, scenario1(), 10), b = solveFull(X, Y, w2, d);
    std::fprintf(stderr, "[demo] Scenario 1 optimum: score=%d  %s (pos %d) / %s (pos %d)\n",
                 a.score, X.substr(a.s1, a.len).c_str(), a.s1 + 1, Y.substr(a.s2, a.len).c_str(), a.s2 + 1);
    std::fprintf(stderr, "[demo] Scenario 2 (delta=%d) optimum: score=%d  %s (pos %d) / %s (pos %d)\n",
                 d, b.score, X.substr(b.s1, b.len).c_str(), b.s1 + 1, Y.substr(b.s2, b.len).c_str(), b.s2 + 1);
}
bool selfTest() {
    std::mt19937_64 rng(7);
    int fails = 0;
    for (int it = 0; it < 500; ++it) {
        int n = 1 + (int)(rng() % 60), m = 1 + (int)(rng() % 60);
        bool smallAlpha = it % 2;                          // 3-letter alphabet -> long runs w/ mismatches
        std::string X = smallAlpha ? genAlpha(n, 3, rng) : genEnglish(n, rng);
        std::string Y = smallAlpha ? genAlpha(m, 3, rng) : genEnglish(m, rng);
        Weights w; int delta;
        if (it % 4 < 2) { w = scenario1(); delta = 10; }
        else            { w = scenario2(); delta = deltaS2((int)(rng() % 12)); }
        Result b = bruteForce(X, Y, w, delta), f = solveFull(X, Y, w, delta), l = solveLinear(X, Y, w, delta);
        bool ok = b.score == f.score && f.score == l.score && f.s1 == l.s1 && f.s2 == l.s2 &&
                  f.len == l.len && scoreOf(X, Y, f.s1, f.s2, f.len, w, delta) == f.score;
        if (!ok) ++fails;
    }
    std::fprintf(stderr, "[self-test] 500 random cases vs brute force: %s (%d failures)\n",
                 fails ? "FAILED" : "all passed", fails);
    return fails == 0;
}

int main(int argc, char** argv) {
    if (argc == 5 && std::strcmp(argv[1], "mem") == 0) {
        int n = std::atoi(argv[2]), m = std::atoi(argv[3]);
        bool full = std::strcmp(argv[4], "full") == 0;
        std::mt19937_64 rng(1);
        std::string X = genEnglish(n, rng), Y = genEnglish(m, rng);
        long rssBefore = currentRssKB();
        std::size_t heapBase = g_heapCur; resetHeapPeak();
        g_sink = full ? solveFull(X, Y, scenario1(), 10).score
                      : solveLinear(X, Y, scenario1(), 10).score;
        std::size_t heapPeak = g_heapPeak - heapBase;  // DP arrays only
        long rssPeak = peakRssKB();
        std::printf("%d,%d,%s,%zu,%zu,%ld,%ld\n", n, m, argv[4],
                    full ? theoryFull(n, m) : theoryLin(n, m), heapPeak, rssBefore, rssPeak);
        return 0;
    }

    demo();
    if (!selfTest()) return 1;
    std::printf("experiment,dist,n,m,delta,t,score,len,start1,start2,mismatches,"
                "time_full_us,time_lin_us,theory_full_bytes,theory_lin_bytes,full_eq_lin,score_ok,sub1,sub2\n");
    std::mt19937_64 rng(20261103ULL);
    const int sizes[][2] = {{100,100},{500,500},{1000,1000},{2000,2000},{4000,4000},{1000,4000},{8000,8000}};
    const char* dn[2] = {"uniform", "english"};
    auto gen = [&](int n, int d) { return d ? genEnglish(n, rng) : genAlpha(n, 26, rng); };

    // Exp 1: Scenario 1 (w=1, delta=10) — running time vs size
    for (auto& s : sizes) for (int d = 0; d < 2; ++d) {
        std::string X = gen(s[0], d), Y = gen(s[1], d);
        runCase("S1", dn[d], X, Y, scenario1(), 10, -1);
    }
    // Exp 2: Scenario 2 — delta sweep on fixed strings (t=1..10 intermediate, 0/11 endpoints)
    for (int d = 0; d < 2; ++d) {
        std::string X = gen(2000, d), Y = gen(2000, d);
        for (int t = 0; t <= 11; ++t) runCase("S2sweep", dn[d], X, Y, scenario2(), deltaS2(t), t);
    }
    // Exp 3: Scenario 2 — running time vs size at middle delta (t=5)
    for (auto& s : sizes) for (int d = 0; d < 2; ++d) {
        std::string X = gen(s[0], d), Y = gen(s[1], d);
        runCase("S2size", dn[d], X, Y, scenario2(), deltaS2(5), 5);
    }
    return 0;
}
