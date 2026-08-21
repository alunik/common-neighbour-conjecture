// Deterministic CREATE-side witness search for 3.Fi22 <= GL(27,4).
//
// Targets are the exact C_Fi22(x)-orbit representatives on the eigenlines
// of one representative x of every prime-order class.  For each target w,
// search along deterministic images r of one regular seed checked by GAP.
// Thus r is regular by group invariance.  A candidate u = w + r is retained
// when it has the reference Hermitian norm and the reference zero count in the
// aligned degree-61776 dual-line orbit.  The retained u is *not* trusted:
// check_second_summands.g subsequently computes its stabilizer exactly.
//
// Usage:
//   find_inner_sums <root> <outfile> <threads> <max-attempts> <survivors>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

using namespace std;

static constexpr int DIM = 27;
static constexpr int NPTS = 61776;
static constexpr int NWORDS = (NPTS + 63) / 64;
static uint8_t MUL[4][4];
static uint8_t CONJ4[4] = {0, 1, 3, 2};
static uint8_t G[2][DIM][DIM], H[DIM][DIM];
static vector<uint64_t> bitA[DIM], bitB[DIM], bitC[DIM];
static uint8_t r0[DIM];
static int referenceNorm = -1;
static int referenceZeros = -1;

[[noreturn]] static void die(const string &message) {
  fprintf(stderr, "FATAL: %s\n", message.c_str());
  exit(1);
}

static void mulInit() {
  for (int x = 0; x < 4; ++x)
    for (int y = 0; y < 4; ++y) {
      int a = x & 1, b = x >> 1, c = y & 1, d = y >> 1;
      int ra = (a & c) ^ (b & d);
      int rb = (a & d) ^ (b & c) ^ (b & d);
      MUL[x][y] = static_cast<uint8_t>(ra | (rb << 1));
    }
}

static string readLine(FILE *file) {
  string line;
  int ch;
  while ((ch = fgetc(file)) != EOF && ch != '\n')
    if (ch != '\r') line.push_back(static_cast<char>(ch));
  return line;
}

static void decodeDigits(const string &digits, uint8_t *vector, int dim) {
  if (static_cast<int>(digits.size()) != dim) die("digit-vector length");
  for (int i = 0; i < dim; ++i) {
    if (digits[i] < '0' || digits[i] > '3') die("invalid GF(4) digit");
    vector[i] = static_cast<uint8_t>(digits[i] - '0');
  }
}

static string digitsOf(const uint8_t *vector, int dim) {
  string result;
  result.reserve(dim);
  for (int i = 0; i < dim; ++i)
    result.push_back(static_cast<char>('0' + vector[i]));
  return result;
}

static void loadCoreData(const string &root) {
  string path = root + "/generated/matrix_generators.txt";
  FILE *file = fopen(path.c_str(), "r");
  if (!file) die("cannot open matrix_generators.txt");
  for (int generator = 0; generator < 2; ++generator)
    for (int row = 0; row < DIM; ++row) {
      string line = readLine(file);
      decodeDigits(line, G[generator][row], DIM);
    }
  fclose(file);

  path = root + "/generated/hermitian_form.txt";
  file = fopen(path.c_str(), "r");
  if (!file) die("cannot open hermitian_form.txt");
  for (int row = 0; row < DIM; ++row) {
    string line = readLine(file);
    decodeDigits(line, H[row], DIM);
  }
  fclose(file);

  for (int coordinate = 0; coordinate < DIM; ++coordinate) {
    bitA[coordinate].assign(NWORDS, 0);
    bitB[coordinate].assign(NWORDS, 0);
    bitC[coordinate].assign(NWORDS, 0);
  }
  path = root + "/generated/dual_orbit.tsv";
  file = fopen(path.c_str(), "r");
  if (!file) die("cannot open dual_orbit.tsv");
  for (int point = 0; point < NPTS; ++point) {
    string line = readLine(file);
    size_t tab = line.find('\t');
    if (tab == string::npos || atoi(line.c_str()) != point + 1)
      die("dual-orbit index mismatch");
    string digits = line.substr(tab + 1);
    uint8_t vector[DIM];
    decodeDigits(digits, vector, DIM);
    for (int coordinate = 0; coordinate < DIM; ++coordinate) {
      if (vector[coordinate] & 1)
        bitA[coordinate][point >> 6] |= 1ULL << (point & 63);
      if (vector[coordinate] & 2)
        bitB[coordinate][point >> 6] |= 1ULL << (point & 63);
    }
  }
  if (!readLine(file).empty()) die("extra dual-orbit row");
  fclose(file);
  for (int coordinate = 0; coordinate < DIM; ++coordinate)
    for (int word = 0; word < NWORDS; ++word)
      bitC[coordinate][word] =
          bitA[coordinate][word] ^ bitB[coordinate][word];

  path = root + "/generated/regular_seed.txt";
  file = fopen(path.c_str(), "r");
  if (!file) die("cannot open regular_seed.txt");
  string seedDigits = readLine(file);
  if (!readLine(file).empty()) die("extra row in regular_seed.txt");
  fclose(file);
  decodeDigits(seedDigits, r0, DIM);
}

static int hermitianNorm(const uint8_t *vector) {
  uint8_t result = 0;
  for (int i = 0; i < DIM; ++i) {
    if (!vector[i]) continue;
    uint8_t inner = 0;
    for (int j = 0; j < DIM; ++j)
      inner ^= MUL[H[i][j]][CONJ4[vector[j]]];
    result ^= MUL[vector[i]][inner];
  }
  return result;
}

struct ZeroWorkspace {
  vector<uint64_t> a, b;
  ZeroWorkspace() : a(NWORDS), b(NWORDS) {}
};

static int zeroCount(const uint8_t *candidate, ZeroWorkspace &workspace) {
  fill(workspace.a.begin(), workspace.a.end(), 0);
  fill(workspace.b.begin(), workspace.b.end(), 0);
  for (int coordinate = 0; coordinate < DIM; ++coordinate) {
    const vector<uint64_t> *sourceA = nullptr;
    const vector<uint64_t> *sourceB = nullptr;
    switch (candidate[coordinate]) {
      case 0:
        continue;
      case 1:
        sourceA = &bitA[coordinate];
        sourceB = &bitB[coordinate];
        break;
      case 2:
        sourceA = &bitB[coordinate];
        sourceB = &bitC[coordinate];
        break;
      case 3:
        sourceA = &bitC[coordinate];
        sourceB = &bitA[coordinate];
        break;
    }
    for (int word = 0; word < NWORDS; ++word) {
      workspace.a[word] ^= (*sourceA)[word];
      workspace.b[word] ^= (*sourceB)[word];
    }
  }
  int count = 0;
  for (int word = 0; word < NWORDS - 1; ++word)
    count += __builtin_popcountll(~(workspace.a[word] | workspace.b[word]));
  constexpr int tail = NPTS & 63;
  uint64_t mask = tail ? ((1ULL << tail) - 1) : ~0ULL;
  count += __builtin_popcountll(
      ~(workspace.a[NWORDS - 1] | workspace.b[NWORDS - 1]) & mask);
  return count;
}

static void stepVector(uint8_t *vector, int generator) {
  uint8_t image[DIM] = {};
  for (int row = 0; row < DIM; ++row) {
    uint8_t scalar = vector[row];
    if (!scalar) continue;
    for (int column = 0; column < DIM; ++column)
      image[column] ^= MUL[scalar][G[generator][row][column]];
  }
  memcpy(vector, image, DIM);
}

struct Target {
  string cls;
  int block;
  string coordinates;
  string digits;
  uint8_t vector[DIM];
};

static string nextToken(FILE *file) {
  string result;
  int ch;
  do ch = fgetc(file); while (ch != EOF && isspace(ch));
  while (ch != EOF && !isspace(ch)) {
    result.push_back(static_cast<char>(ch));
    ch = fgetc(file);
  }
  return result;
}

static vector<vector<string>> loadBases(
    const string &path, vector<int> &dimensions) {
  FILE *file = fopen(path.c_str(), "r");
  if (!file) die("cannot open class file " + path);
  vector<vector<string>> bases;
  int blocks = 0;
  for (;;) {
    string token = nextToken(file);
    if (token.empty()) break;
    if (token == "NBLOCKS") {
      blocks = atoi(nextToken(file).c_str());
      bases.resize(blocks);
    } else if (token == "DIMS") {
      dimensions.resize(blocks);
      for (int block = 0; block < blocks; ++block)
        dimensions[block] = atoi(nextToken(file).c_str());
    } else if (token == "BASIS") {
      int block = atoi(nextToken(file).c_str()) - 1;
      if (block < 0 || block >= blocks || dimensions.empty())
        die("malformed BASIS section");
      for (int row = 0; row < dimensions[block]; ++row) {
        string digits = nextToken(file);
        if (static_cast<int>(digits.size()) != DIM)
          die("basis-row length in " + path);
        bases[block].push_back(digits);
      }
    } else {
      // All other records occupy their current line.  BASIS is the only
      // multiline record needed by this reader.
      readLine(file);
    }
  }
  fclose(file);
  return bases;
}

static void coordinatesToVector(
    const vector<string> &basis, const string &coordinates,
    uint8_t *vector) {
  memset(vector, 0, DIM);
  if (basis.size() != coordinates.size())
    die("coordinate/basis dimension mismatch");
  for (size_t row = 0; row < basis.size(); ++row) {
    uint8_t scalar = static_cast<uint8_t>(coordinates[row] - '0');
    if (!scalar) continue;
    for (int column = 0; column < DIM; ++column)
      vector[column] ^=
          MUL[scalar][static_cast<uint8_t>(basis[row][column] - '0')];
  }
}

static void addTarget(
    vector<Target> &targets, unordered_set<string> &seen,
    const string &cls, int block, const string &coordinates,
    const vector<string> &basis) {
  Target target;
  target.cls = cls;
  target.block = block;
  target.coordinates = coordinates;
  coordinatesToVector(basis, coordinates, target.vector);
  target.digits = digitsOf(target.vector, DIM);
  if (target.digits.find_first_not_of('0') == string::npos)
    die("zero target");
  if (seen.insert(target.digits).second) targets.push_back(target);
}

static vector<Target> collectTargets(const string &root) {
  vector<Target> targets;
  unordered_set<string> seen;
  const vector<string> classes =
      {"2A", "2B", "2C", "3A", "3B", "3C", "3D", "5A"};
  for (const string &cls : classes) {
    string classPath = root + "/generated/cls_" + cls + ".txt";
    string orbitPath = root + "/generated/orbits_" + cls + ".txt";
    vector<int> dimensions;
    vector<vector<string>> bases = loadBases(classPath, dimensions);
    FILE *file = fopen(orbitPath.c_str(), "r");
    if (!file) die("cannot open orbit file " + orbitPath);
    for (;;) {
      string line = readLine(file);
      if (line.empty() && feof(file)) break;
      if (line.rfind("ORBIT", 0) != 0) continue;
      int block;
      unsigned long long seedIndex, orbitSize;
      char coordinates[128];
      if (sscanf(line.c_str(), "ORBIT %d %llu %llu %127s",
                 &block, &seedIndex, &orbitSize, coordinates) != 4)
        die("malformed orbit row in " + orbitPath);
      if (block < 1 || block > static_cast<int>(bases.size()))
        die("orbit block out of range");
      addTarget(targets, seen, cls, block, coordinates, bases[block - 1]);
    }
    fclose(file);
  }

  const vector<string> tiny = {"7A", "11A", "11B", "13A", "13B"};
  for (const string &cls : tiny) {
    vector<int> dimensions;
    vector<vector<string>> bases = loadBases(
        root + "/generated/cls_" + cls + "_tiny.txt", dimensions);
    for (size_t block = 0; block < bases.size(); ++block) {
      int dimension = dimensions[block];
      for (int pivot = 0; pivot < dimension; ++pivot) {
        uint64_t tails = 1ULL << (2 * (dimension - 1 - pivot));
        for (uint64_t tail = 0; tail < tails; ++tail) {
          string coordinates(dimension, '0');
          coordinates[pivot] = '1';
          uint64_t value = tail;
          for (int position = pivot + 1; position < dimension; ++position) {
            coordinates[position] =
                static_cast<char>('0' + (value & 3));
            value >>= 2;
          }
          addTarget(targets, seen, cls, static_cast<int>(block) + 1,
                    coordinates, bases[block]);
        }
      }
    }
  }
  if (targets.size() != 5363)
    die("target count is " + to_string(targets.size()) +
        ", expected 5363");
  return targets;
}

static uint64_t splitmix64(uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

struct Survivor {
  int epoch;
  int attempt;
  int globalAttempt;
  string word;
  string r;
  string u;
  int norm;
  int zeros;
};

struct SearchResult {
  vector<Survivor> survivors;
  int attempts = 0;
};

static SearchResult searchTarget(
    const Target &target, size_t targetIndex, int maxAttempts, int want,
    ZeroWorkspace &workspace) {
  constexpr int EPOCH = 250;
  constexpr int BURN = 48;
  SearchResult result;
  uint8_t r[DIM], u[DIM];
  for (int epoch = 0;
       result.attempts < maxAttempts &&
           static_cast<int>(result.survivors.size()) < want;
       ++epoch) {
    memcpy(r, r0, DIM);
    uint64_t state = splitmix64(
        0x464932325749544eULL ^
        (static_cast<uint64_t>(targetIndex) << 24) ^
        static_cast<uint64_t>(epoch));
    string word;
    word.reserve(BURN + 2 * EPOCH);
    for (int step = 0; step < BURN; ++step) {
      state = splitmix64(state);
      int generator = static_cast<int>(state & 1);
      word.push_back(static_cast<char>('1' + generator));
      stepVector(r, generator);
    }
    for (int attempt = 1;
         attempt <= EPOCH && result.attempts < maxAttempts &&
             static_cast<int>(result.survivors.size()) < want;
         ++attempt) {
      for (int step = 0; step < 2; ++step) {
        state = splitmix64(state);
        int generator = static_cast<int>(state & 1);
        word.push_back(static_cast<char>('1' + generator));
        stepVector(r, generator);
      }
      ++result.attempts;
      bool nonzero = false;
      for (int coordinate = 0; coordinate < DIM; ++coordinate) {
        u[coordinate] = r[coordinate] ^ target.vector[coordinate];
        nonzero |= u[coordinate] != 0;
      }
      if (!nonzero) continue;
      int norm = hermitianNorm(u);
      if (norm != referenceNorm) continue;
      int zeros = zeroCount(u, workspace);
      if (zeros != referenceZeros) continue;
      result.survivors.push_back(
          {epoch, attempt, result.attempts, word,
           digitsOf(r, DIM), digitsOf(u, DIM), norm, zeros});
    }
  }
  return result;
}

int main(int argc, char **argv) {
  if (argc != 6) {
    fprintf(stderr,
        "usage: %s <root> <outfile> <threads> <max-attempts> <survivors>\n",
        argv[0]);
    return 1;
  }
  const string root = argv[1];
  const string outfile = argv[2];
  int threadCount = atoi(argv[3]);
  int maxAttempts = atoi(argv[4]);
  int want = atoi(argv[5]);
  if (threadCount < 1 || maxAttempts < 1 || want < 1)
    die("invalid numeric argument");

  mulInit();
  loadCoreData(root);
  ZeroWorkspace seedWorkspace;
  referenceNorm = hermitianNorm(r0);
  referenceZeros = zeroCount(r0, seedWorkspace);
  if (referenceNorm != 1 || referenceZeros != 15424)
    die("reference regular-seed invariants mismatch");
  vector<Target> targets = collectTargets(root);
  fprintf(stderr,
      "loaded %zu targets; reference seed norm=%d zeros=%d\n",
      targets.size(), referenceNorm, referenceZeros);

  vector<SearchResult> results(targets.size());
  atomic<size_t> next{0};
  atomic<size_t> completed{0};
  mutex logMutex;
  vector<thread> threads;
  for (int threadIndex = 0; threadIndex < threadCount; ++threadIndex) {
    threads.emplace_back([&]() {
      ZeroWorkspace workspace;
      for (;;) {
        size_t targetIndex = next.fetch_add(1);
        if (targetIndex >= targets.size()) break;
        results[targetIndex] = searchTarget(
            targets[targetIndex], targetIndex, maxAttempts, want, workspace);
        size_t done = completed.fetch_add(1) + 1;
        if (done % 100 == 0 || done == targets.size()) {
          lock_guard<mutex> lock(logMutex);
          fprintf(stderr, "progress %zu/%zu\n", done, targets.size());
        }
      }
    });
  }
  for (thread &worker : threads) worker.join();

  FILE *output = fopen(outfile.c_str(), "w");
  if (!output) die("cannot open output file");
  fprintf(output, "# FI22_WITNESS_SCREEN_V1\n");
  fprintf(output,
      "# tid cls block coordinates target survivor epoch attempt "
      "global_attempt norm zeros word r u\n");
  size_t failures = 0, survivorCount = 0;
  long long totalAttempts = 0;
  for (size_t targetIndex = 0; targetIndex < targets.size(); ++targetIndex) {
    const Target &target = targets[targetIndex];
    const SearchResult &result = results[targetIndex];
    totalAttempts += result.attempts;
    if (result.survivors.empty()) ++failures;
    for (size_t survivorIndex = 0;
         survivorIndex < result.survivors.size(); ++survivorIndex) {
      const Survivor &survivor = result.survivors[survivorIndex];
      fprintf(output,
          "%zu %s %d %s %s %zu %d %d %d %d %d %s %s %s\n",
          targetIndex, target.cls.c_str(), target.block,
          target.coordinates.c_str(), target.digits.c_str(), survivorIndex,
          survivor.epoch, survivor.attempt, survivor.globalAttempt,
          survivor.norm, survivor.zeros, survivor.word.c_str(),
          survivor.r.c_str(), survivor.u.c_str());
      ++survivorCount;
    }
  }
  fprintf(output,
      "# done targets=%zu survivors=%zu failures=%zu attempts=%lld\n",
      targets.size(), survivorCount, failures, totalAttempts);
  fclose(output);
  fprintf(stderr,
      "DONE targets=%zu survivors=%zu failures=%zu attempts=%lld\n",
      targets.size(), survivorCount, failures, totalAttempts);
  return failures == 0 ? 0 : 2;
}
