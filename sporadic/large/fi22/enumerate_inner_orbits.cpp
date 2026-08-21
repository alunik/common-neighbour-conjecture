// orbit_bfs.cpp — generic eigenline-orbit BFS for the 3.Fi22 census.
//
// Input: a cls_XX.txt file describing an F4-linear action of centralizer
// generators on a disjoint union of eigenspace blocks (dims d_b <= 21),
// possibly permuting blocks.  Enumerates ALL orbits of the generated group
// on the union of the projective-line sets P(F4^{d_b}), by BFS over
// canonical line representatives (first nonzero coordinate = 1).
//
// Output: one line per orbit: "ORBIT <block> <seedIndexInBlock> <size> <digits>"
// plus a final "TOTAL <orbits> <lines>" (asserted == sum of block line counts).
//
// GF(4) code: 0,1,2,3 = 0, 1, w, w^2  (bit0 = a, bit1 = b for x = a + b w).
// Packed vectors: coordinate i in bits (2i, 2i+1), dim <= 21 -> 42 bits.
//
// Usage: orbit_bfs <clsfile> <outfile> <nthreads>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <chrono>
using namespace std;

static uint8_t MUL[4][4];
static void mulInit() {
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 4; y++) {
      int a = x & 1, b = x >> 1, c = y & 1, d = y >> 1;
      int ra = (a & c) ^ (b & d);
      int rb = (a & d) ^ (b & c) ^ (b & d);
      MUL[x][y] = (uint8_t)(ra | (rb << 1));
    }
}

static void die(const char *m) { fprintf(stderr, "FATAL: %s\n", m); exit(1); }

// ---------- packed-vector helpers ----------
static const uint64_t LOWMASK = 0x5555555555555555ULL;
static inline uint64_t mulW(uint64_t v) {   // multiply every coord by w
  uint64_t A = v & LOWMASK, B = (v >> 1) & LOWMASK;
  return B | ((A ^ B) << 1);
}
static inline uint64_t mulW2(uint64_t v) {  // multiply every coord by w^2
  uint64_t A = v & LOWMASK, B = (v >> 1) & LOWMASK;
  return (A ^ B) | (A << 1);
}
// canonicalize nonzero v (first nonzero coord -> 1); also return pivot k
static inline uint64_t canon(uint64_t v, int *kout) {
  uint64_t nz = (v | (v >> 1)) & LOWMASK;
  int tz = __builtin_ctzll(nz);
  int k = tz >> 1;
  *kout = k;
  int c = (int)((v >> (2 * k)) & 3);
  if (c == 2) v = mulW2(v);       // c = w   -> multiply by w^-1 = w^2
  else if (c == 3) v = mulW(v);   // c = w^2 -> multiply by w
  return v;
}

// ---------- class-file structures ----------
struct Block { int dim; uint64_t nlines, offset; vector<uint64_t> pref; };
struct GenBlock { int tgt; uint64_t tab[3][16384]; int nchunk; };
struct Gen { vector<GenBlock> gb; };

static int NB, NG;
static vector<Block> B;
static vector<Gen> G;
static string clsName;
static long long centOrder = -1;

// canonical line index within block (dim d): pivot k, tail = v >> (2k+2)
// idx = pref[k] + tail ; pref[k] = sum_{j<k} 4^{d-1-j}
static inline uint64_t encodeIdx(const Block &b, uint64_t v, int k) {
  return b.pref[k] + (v >> (2 * k + 2));
}
static inline uint64_t decodeIdx(const Block &b, uint64_t idx, int *kout) {
  int k = 0;
  // linear scan fine (d <= 21); pref increasing
  while (k + 1 < b.dim && idx >= b.pref[k + 1]) k++;
  *kout = k;
  uint64_t tail = idx - b.pref[k];
  return (1ULL << (2 * k)) | (tail << (2 * k + 2));
}

// ---------- parsing ----------
static string tok(FILE *f) {
  string s; int c;
  while ((c = fgetc(f)) != EOF && isspace(c)) {}
  if (c == EOF) return s;
  s.push_back((char)c);
  while ((c = fgetc(f)) != EOF && !isspace(c)) s.push_back((char)c);
  return s;
}

static void loadCls(const char *path) {
  FILE *f = fopen(path, "r");
  if (!f) die("cannot open cls file");
  string t;
  vector<vector<vector<uint8_t>>> mats; // per gen per block rows(strings)
  vector<vector<int>> bperm;
  while (!(t = tok(f)).empty()) {
    if (t == "CLASS") { clsName = tok(f); }
    else if (t == "CENT") { centOrder = atoll(tok(f).c_str()); }
    else if (t == "NBLOCKS") {
      NB = atoi(tok(f).c_str());
      B.resize(NB);
    } else if (t == "DIMS") {
      for (int j = 0; j < NB; j++) {
        B[j].dim = atoi(tok(f).c_str());
        if (B[j].dim < 1 || B[j].dim > 21) die("dim out of range");
      }
    } else if (t == "NGENS") {
      NG = atoi(tok(f).c_str());
      mats.assign(NG, {});
      bperm.assign(NG, {});
    } else if (t == "BASIS") {
      int j = atoi(tok(f).c_str());
      // skip dim rows of 27 digits
      for (int r = 0; r < B[j - 1].dim; r++) tok(f);
    } else if (t == "GEN") {
      int gi = atoi(tok(f).c_str());
      string bp = tok(f);
      if (bp != "BLOCKPERM") die("expected BLOCKPERM");
      bperm[gi - 1].resize(NB);
      for (int j = 0; j < NB; j++) bperm[gi - 1][j] = atoi(tok(f).c_str()) - 1;
      mats[gi - 1].resize(NB);
      for (int j = 0; j < NB; j++) {
        int ds = B[j].dim, dt = B[bperm[gi - 1][j]].dim;
        if (ds != dt) die("block perm dims mismatch");
        mats[gi - 1][j].resize(ds * dt);
        for (int r = 0; r < ds; r++) {
          string row = tok(f);
          if ((int)row.size() != dt) die("matrix row length");
          for (int c = 0; c < dt; c++)
            mats[gi - 1][j][r * dt + c] = (uint8_t)(row[c] - '0');
        }
      }
    } else die("unknown token in cls file");
  }
  fclose(f);
  if (NB <= 0 || NG <= 0) die("missing NBLOCKS/NGENS");
  // append inverse generators (forward-only BFS is correct but can be slow):
  // invert each block matrix by Gauss-Jordan over GF(4); inverse blockperm.
  {
    static const uint8_t INV4[4] = {0, 1, 3, 2};
    int ng0 = NG;
    for (int g = 0; g < ng0; g++) {
      vector<vector<uint8_t>> im(NB);
      vector<int> ibp(NB);
      for (int j = 0; j < NB; j++) {
        int tgt = bperm[g][j];
        ibp[tgt] = j;
        int d = B[j].dim;
        // invert mats[g][j] (d x d): result maps block tgt back to block j
        vector<uint8_t> A(mats[g][j]);
        vector<uint8_t> I(d * d, 0);
        for (int r = 0; r < d; r++) I[r * d + r] = 1;
        for (int col = 0; col < d; col++) {
          int piv = -1;
          for (int r = col; r < d; r++) if (A[r * d + col]) { piv = r; break; }
          if (piv < 0) die("singular generator matrix");
          if (piv != col) {
            for (int c = 0; c < d; c++) {
              swap(A[piv * d + c], A[col * d + c]);
              swap(I[piv * d + c], I[col * d + c]);
            }
          }
          uint8_t pv = INV4[A[col * d + col]];
          for (int c = 0; c < d; c++) {
            A[col * d + c] = MUL[pv][A[col * d + c]];
            I[col * d + c] = MUL[pv][I[col * d + c]];
          }
          for (int r = 0; r < d; r++) {
            if (r == col) continue;
            uint8_t f2 = A[r * d + col];
            if (!f2) continue;
            for (int c = 0; c < d; c++) {
              A[r * d + c] ^= MUL[f2][A[col * d + c]];
              I[r * d + c] ^= MUL[f2][I[col * d + c]];
            }
          }
        }
        im[tgt] = I;   // acts on block tgt with target j
      }
      if (ibp == bperm[g] && im == mats[g]) continue;  // involution: skip
      bperm.push_back(ibp);
      mats.push_back(im);
      NG++;
    }
    fprintf(stderr, "generators: %d (with inverses; %d original)\n", NG, ng0);
  }
  // offsets and prefix tables
  uint64_t off = 0;
  for (int j = 0; j < NB; j++) {
    int d = B[j].dim;
    B[j].pref.resize(d);
    uint64_t s = 0;
    for (int k = 0; k < d; k++) { B[j].pref[k] = s; s += 1ULL << (2 * (d - 1 - k)); }
    B[j].nlines = s / 3 * 3 == s ? s : s; // s = (4^d-1)/3? no: s = sum 4^{d-1-k} = (4^d-1)/3
    B[j].nlines = s;
    B[j].offset = off;
    off += s;
  }
  // build generator tables
  G.resize(NG);
  for (int g = 0; g < NG; g++) {
    G[g].gb.resize(NB);
    for (int j = 0; j < NB; j++) {
      GenBlock &gb = G[g].gb[j];
      gb.tgt = bperm[g][j];
      int ds = B[j].dim, dt = B[gb.tgt].dim;
      gb.nchunk = (ds + 6) / 7;
      // packed rows times scalar
      vector<uint64_t> rowmul(ds * 4);
      for (int r = 0; r < ds; r++) {
        for (int c = 0; c < 4; c++) {
          uint64_t p = 0;
          for (int cc = 0; cc < dt; cc++) {
            uint8_t val = MUL[c][mats[g][j][r * dt + cc]];
            p |= ((uint64_t)val) << (2 * cc);
          }
          rowmul[r * 4 + c] = p;
        }
      }
      for (int ch = 0; ch < gb.nchunk; ch++) {
        uint64_t *T = gb.tab[ch];
        T[0] = 0;
        int base = ch * 7;
        int width = ds - base; if (width > 7) width = 7;
        for (int p = 1; p < (1 << (2 * width)); p++) {
          int j2 = __builtin_ctz(p) >> 1;      // lowest nonzero coord slot
          int digit = (p >> (2 * j2)) & 3;
          int q = p & ~(3 << (2 * j2));
          T[p] = T[q] ^ rowmul[(base + j2) * 4 + digit];
        }
        for (int p = 1 << (2 * width); p < 16384; p++) T[p] = ~0ULL; // poison
      }
    }
  }
}

// apply generator g to (block j, packed v) -> (block, packed) not canonical
static inline uint64_t applyGen(int g, int j, uint64_t v, int *jout) {
  const GenBlock &gb = G[g].gb[j];
  *jout = gb.tgt;
  uint64_t r = gb.tab[0][v & 16383];
  if (gb.nchunk > 1) r ^= gb.tab[1][(v >> 14) & 16383];
  if (gb.nchunk > 2) r ^= gb.tab[2][(v >> 28) & 16383];
  return r;
}

// ---------- BFS ----------
static uint64_t TOTAL;
static uint64_t *visited;   // bitmap over global indices
static inline bool testAndSet(uint64_t gidx) {
  uint64_t bit = 1ULL << (gidx & 63);
  uint64_t old = __atomic_fetch_or(&visited[gidx >> 6], bit, __ATOMIC_RELAXED);
  return (old & bit) != 0;   // true if already visited
}

struct Item { uint32_t blk; uint64_t vec; };

int main(int argc, char **argv) {
  mulInit();
  if (argc != 4) { fprintf(stderr, "usage: %s cls.txt out.txt nthreads\n", argv[0]); return 1; }
  loadCls(argv[1]);
  int NT = atoi(argv[3]);
  if (NT < 1) NT = 1;
  TOTAL = 0;
  for (int j = 0; j < NB; j++) TOTAL += B[j].nlines;
  fprintf(stderr, "class %s: %d blocks, %d gens, TOTAL lines %llu\n",
          clsName.c_str(), NB, NG, (unsigned long long)TOTAL);
  uint64_t words = (TOTAL + 63) / 64;
  visited = (uint64_t *)calloc(words, 8);
  if (!visited) die("visited alloc failed");
  FILE *out = fopen(argv[2], "w");
  if (!out) die("outfile");
  fprintf(out, "# orbit_bfs class %s cent %lld total %llu\n",
          clsName.c_str(), centOrder, (unsigned long long)TOTAL);

  uint64_t nOrbits = 0, nVisited = 0;
  vector<Item> cur, nxt;
  auto t0 = chrono::steady_clock::now();

  // sweep global indices for unvisited seeds
  uint64_t sweepWord = 0;
  int curBlk = 0;
  while (true) {
    // find next unvisited global index
    while (sweepWord < words && visited[sweepWord] == ~0ULL) sweepWord++;
    if (sweepWord >= words) break;
    uint64_t base = sweepWord << 6;
    uint64_t w = visited[sweepWord];
    int bit = __builtin_ctzll(~w);
    uint64_t gseed = base + bit;
    if (gseed >= TOTAL) break;
    // decode seed
    while (curBlk + 1 < NB && gseed >= B[curBlk + 1].offset) curBlk++;
    while (gseed < B[curBlk].offset) curBlk--;
    int blk = curBlk;
    int kk;
    uint64_t v = decodeIdx(B[blk], gseed - B[blk].offset, &kk);
    // BFS this orbit
    testAndSet(gseed);
    uint64_t osize = 1;
    cur.clear(); nxt.clear();
    cur.push_back({(uint32_t)blk, v});
    while (!cur.empty()) {
      size_t n = cur.size();
      nxt.clear();
      if (n < 20000 || NT == 1) {
        for (size_t i = 0; i < n; i++) {
          for (int g = 0; g < NG; g++) {
            int jb;
            uint64_t img = applyGen(g, cur[i].blk, cur[i].vec, &jb);
            int k2;
            img = canon(img, &k2);
            uint64_t gidx = B[jb].offset + encodeIdx(B[jb], img, k2);
            if (!testAndSet(gidx)) nxt.push_back({(uint32_t)jb, img});
          }
        }
      } else {
        vector<vector<Item>> parts(NT);
        vector<thread> th;
        size_t chunk = (n + NT - 1) / NT;
        for (int tI = 0; tI < NT; tI++) {
          th.emplace_back([&, tI]() {
            size_t lo = tI * chunk, hi = lo + chunk;
            if (hi > n) hi = n;
            auto &mine = parts[tI];
            for (size_t i = lo; i < hi; i++) {
              for (int g = 0; g < NG; g++) {
                int jb;
                uint64_t img = applyGen(g, cur[i].blk, cur[i].vec, &jb);
                int k2;
                img = canon(img, &k2);
                uint64_t gidx = B[jb].offset + encodeIdx(B[jb], img, k2);
                if (!testAndSet(gidx)) mine.push_back({(uint32_t)jb, img});
              }
            }
          });
        }
        for (auto &x : th) x.join();
        size_t tot = 0;
        for (auto &p : parts) tot += p.size();
        nxt.reserve(tot);
        for (auto &p : parts) nxt.insert(nxt.end(), p.begin(), p.end());
      }
      osize += nxt.size();
      cur.swap(nxt);
    }
    nOrbits++; nVisited += osize;
    // rep digits
    string digits;
    for (int i = 0; i < B[blk].dim; i++)
      digits.push_back((char)('0' + ((v >> (2 * i)) & 3)));
    fprintf(out, "ORBIT %d %llu %llu %s\n", blk + 1,
            (unsigned long long)(gseed - B[blk].offset),
            (unsigned long long)osize, digits.c_str());
    fflush(out);
    if (nOrbits % 50 == 0) {
      auto dt = chrono::duration_cast<chrono::seconds>(
        chrono::steady_clock::now() - t0).count();
      fprintf(stderr, "progress: orbits=%llu visited=%llu (%.2f%%) t=%llds\n",
              (unsigned long long)nOrbits, (unsigned long long)nVisited,
              100.0 * nVisited / TOTAL, (long long)dt);
    }
  }
  if (nVisited != TOTAL) {
    fprintf(out, "ERROR visited %llu != total %llu\n",
            (unsigned long long)nVisited, (unsigned long long)TOTAL);
    die("coverage mismatch");
  }
  fprintf(out, "TOTAL %llu %llu\n", (unsigned long long)nOrbits,
          (unsigned long long)nVisited);
  fclose(out);
  auto dt = chrono::duration_cast<chrono::seconds>(
    chrono::steady_clock::now() - t0).count();
  fprintf(stderr, "DONE class %s: %llu orbits, %llu lines, %llds\n",
          clsName.c_str(), (unsigned long long)nOrbits,
          (unsigned long long)nVisited, (long long)dt);
  return 0;
}
