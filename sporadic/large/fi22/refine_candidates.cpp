// Exact equivariant color-refinement result builder for
// H = 3.Fi22.2 <= GL(54,2), using the degree-3510 Fischer action.
//
// The GAP exporter proves that the complete initial 0/2016/2079 value-orbit
// labels are
// equivariant for the full preimage of every point stabilizer, including
// the normal C3 kernel, and checks both standard-generator edges at all
// 3510 points.  This program:
//   * transports each F4^27 survivor to the ATLAS F2^54 module using the
//     exact exported restriction-of-scalars contribution table;
//   * forms the 3-coloring of the 3510 Fischer points;
//   * repeatedly refines by the equivariant key
//       (own color, sum of neighbor colors, sum of their squares);
//   * exports singleton refined cells.
//
// Every vector stabilizer fixes each singleton point.  GAP downstream
// computes that pointwise container in Fi22.2, lifts/enumerates it in
// 3.Fi22.2, and checks the original 54-bit vector element by element.

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using namespace std;

static constexpr int F4DIM = 27;
static constexpr int H54DIM = 54;
static constexpr int NPTS = 3510;
static constexpr int LOCALDIM = 12;
static constexpr int NROWS = NPTS * LOCALDIM;
static constexpr int NWORDS = (NROWS + 63) / 64 + 1;
static constexpr int VALENCY = 693;
static constexpr int MIN_SINGLETONS = 20;

static string DATA = "sporadic/large/fi22/generated/fischer/";
static uint64_t columnBits[H54DIM][NWORDS];
static uint16_t valueLabel[1 << LOCALDIM];
static int32_t neighbors[NPTS][VALENCY];
static uint64_t bridgeContribution[F4DIM][4];
static vector<int> quotientBase;

static void die(const char *message) {
  fprintf(stderr, "FATAL: %s\n", message);
  exit(1);
}

static string readLine(FILE *file) {
  string line;
  int character;
  while ((character = fgetc(file)) != EOF && character != '\n')
    line.push_back((char)character);
  return line;
}

static void initializePath() {
  const char *value = getenv("FI22_FISCHER_DATA");
  if (value && *value) DATA = value;
  if (DATA.back() != '/') DATA.push_back('/');
}

static uint64_t bits54(const string &digits) {
  if ((int)digits.size() != H54DIM) die("54-bit string has wrong length");
  uint64_t value = 0;
  for (int position = 0; position < H54DIM; position++) {
    if (digits[position] != '0' && digits[position] != '1')
      die("non-binary digit in 54-bit string");
    if (digits[position] == '1') value |= 1ULL << position;
  }
  return value;
}

static void loadData() {
  char path[1024];
  snprintf(path, sizeof path, "%sbases.txt", DATA.c_str());
  FILE *file = fopen(path, "r");
  if (!file) die("cannot open bases.txt");
  memset(columnBits, 0, sizeof columnBits);
  for (int point = 0; point < NPTS; point++) {
    string line = readLine(file);
    if ((int)line.size() != LOCALDIM * H54DIM)
      die("ordered-basis row has wrong length");
    for (int local = 0; local < LOCALDIM; local++) {
      int packedRow = point * LOCALDIM + local;
      for (int coordinate = 0; coordinate < H54DIM; coordinate++) {
        char digit = line[local * H54DIM + coordinate];
        if (digit != '0' && digit != '1') die("non-binary basis digit");
        if (digit == '1')
          columnBits[coordinate][packedRow >> 6] |=
              1ULL << (packedRow & 63);
      }
    }
  }
  fclose(file);

  snprintf(path, sizeof path, "%sfunctionals.txt", DATA.c_str());
  file = fopen(path, "r");
  if (!file) die("cannot open functionals.txt");
  fill(begin(valueLabel), end(valueLabel), UINT16_MAX);
  valueLabel[0] = 0;
  int loadedValues = 0;
  for (;;) {
    string line = readLine(file);
    if (line.empty() && feof(file)) break;
    if (line.empty()) continue;
    char digits[32];
    int label;
    if (sscanf(line.c_str(), "%31s %d", digits, &label) != 2)
      die("functional-label row parse");
    if ((int)strlen(digits) != LOCALDIM || label < 1 || label > 2)
      die("functional-label row range");
    int code = 0;
    for (int position = 0; position < LOCALDIM; position++) {
      if (digits[position] != '0' && digits[position] != '1')
        die("non-binary functional digit");
      if (digits[position] == '1') code |= 1 << position;
    }
    if (code == 0 || valueLabel[code] != UINT16_MAX)
      die("duplicate or zero functional label");
    valueLabel[code] = (uint16_t)label;
    loadedValues++;
  }
  fclose(file);
  if (loadedValues != 4095) die("functional-label table is incomplete");
  for (int code = 0; code < (1 << LOCALDIM); code++)
    if (valueLabel[code] == UINT16_MAX)
      die("functional-label lookup has a hole");

  snprintf(path, sizeof path, "%sneighbors.txt", DATA.c_str());
  file = fopen(path, "r");
  if (!file) die("cannot open neighbors.txt");
  for (int point = 0; point < NPTS; point++) {
    for (int position = 0; position < VALENCY; position++) {
      int neighbor;
      if (fscanf(file, "%d", &neighbor) != 1)
        die("Fischer-neighbor row parse");
      if (neighbor < 0 || neighbor >= NPTS)
        die("Fischer-neighbor index out of range");
      neighbors[point][position] = neighbor;
    }
  }
  fclose(file);

  memset(bridgeContribution, 0, sizeof bridgeContribution);
  snprintf(path, sizeof path, "%sf4_to_h54.tsv", DATA.c_str());
  ifstream bridge(path);
  if (!bridge) die("cannot open f4_to_h54.tsv");
  string line;
  if (!getline(bridge, line)) die("bridge table has no header");
  int loadedContributions = 0;
  while (getline(bridge, line)) {
    if (line.empty()) continue;
    istringstream input(line);
    int coordinate, digit;
    string hbits;
    if (!(input >> coordinate >> digit >> hbits))
      die("bridge contribution row parse");
    if (coordinate < 1 || coordinate > F4DIM ||
        digit < 1 || digit > 3 ||
        bridgeContribution[coordinate - 1][digit] != 0)
      die("bridge contribution row range or duplicate");
    bridgeContribution[coordinate - 1][digit] = bits54(hbits);
    if (bridgeContribution[coordinate - 1][digit] == 0)
      die("bridge contribution unexpectedly zero");
    loadedContributions++;
  }
  if (loadedContributions != F4DIM * 3)
    die("bridge contribution table is incomplete");

  snprintf(path, sizeof path, "%sbaseB.txt", DATA.c_str());
  file = fopen(path, "r");
  if (!file) die("cannot open baseB.txt");
  int point;
  while (fscanf(file, "%d", &point) == 1) {
    if (point < 1 || point > NPTS) die("quotient-base point out of range");
    quotientBase.push_back(point);
  }
  fclose(file);
  if (quotientBase.empty()) die("quotient base is empty");
}

static uint64_t h54FromF4(const string &digits) {
  if ((int)digits.size() != F4DIM) die("F4 vector has wrong length");
  uint64_t vector = 0;
  for (int coordinate = 0; coordinate < F4DIM; coordinate++) {
    if (digits[coordinate] < '0' || digits[coordinate] > '3')
      die("F4 vector digit out of range");
    int digit = digits[coordinate] - '0';
    if (digit) vector ^= bridgeContribution[coordinate][digit];
  }
  if (vector >> H54DIM) die("bridge produced a bit beyond dimension 54");
  return vector;
}

struct WorkBuffers {
  uint64_t accumulated[NWORDS];
  uint32_t initial[NPTS], current[NPTS], next[NPTS];
};

static void initialColors(uint64_t vector, WorkBuffers &work,
                          int counts[3]) {
  memset(work.accumulated, 0, sizeof work.accumulated);
  for (int coordinate = 0; coordinate < H54DIM; coordinate++) {
    if (!((vector >> coordinate) & 1ULL)) continue;
    for (int word = 0; word < NWORDS; word++)
      work.accumulated[word] ^= columnBits[coordinate][word];
  }
  counts[0] = counts[1] = counts[2] = 0;
  for (int point = 0; point < NPTS; point++) {
    int offset = point * LOCALDIM;
    int word = offset >> 6;
    int shift = offset & 63;
    uint64_t code = work.accumulated[word] >> shift;
    if (shift > 64 - LOCALDIM)
      code |= work.accumulated[word + 1] << (64 - shift);
    code &= (1ULL << LOCALDIM) - 1;
    uint16_t label = valueLabel[code];
    counts[label]++;
    work.initial[point] = label + 1;
  }
}

struct RefinementKey {
  uint32_t own;
  uint64_t sum, squareSum;
  uint32_t point;
};

static int refine(WorkBuffers &work) {
  bool initialSeen[4] = {false, false, false, false};
  for (int point = 0; point < NPTS; point++) {
    work.current[point] = work.initial[point];
    initialSeen[work.initial[point]] = true;
  }
  int previousCells =
      initialSeen[1] + initialSeen[2] + initialSeen[3];
  static thread_local vector<RefinementKey> keys(NPTS);
  for (int round = 0; round < 8; round++) {
    for (int point = 0; point < NPTS; point++) {
      uint64_t sum = 0, squareSum = 0;
      for (int position = 0; position < VALENCY; position++) {
        uint64_t color = work.current[neighbors[point][position]];
        sum += color;
        squareSum += color * color;
      }
      keys[point] = {work.current[point], sum, squareSum,
                     (uint32_t)point};
    }
    sort(keys.begin(), keys.end(),
         [](const RefinementKey &left, const RefinementKey &right) {
           if (left.own != right.own) return left.own < right.own;
           if (left.sum != right.sum) return left.sum < right.sum;
           return left.squareSum < right.squareSum;
         });
    int cells = 0;
    for (int position = 0; position < NPTS; position++) {
      if (position == 0 ||
          keys[position].own != keys[position - 1].own ||
          keys[position].sum != keys[position - 1].sum ||
          keys[position].squareSum != keys[position - 1].squareSum)
        cells++;
      work.next[keys[position].point] = cells;
    }
    if (cells == previousCells) return cells;
    for (int point = 0; point < NPTS; point++)
      work.current[point] = work.next[point];
    previousCells = cells;
  }
  return previousCells;
}

static void singletonCells(WorkBuffers &work, int cells,
                           vector<int> &singletonPoints) {
  static thread_local vector<uint32_t> sizes;
  static thread_local vector<int> cellPoint;
  sizes.assign(cells + 1, 0);
  cellPoint.assign(cells + 1, -1);
  for (int point = 0; point < NPTS; point++) {
    sizes[work.current[point]]++;
    cellPoint[work.current[point]] = point + 1;
  }
  singletonPoints.clear();
  for (int cell = 1; cell <= cells; cell++)
    if (sizes[cell] == 1) singletonPoints.push_back(cellPoint[cell]);
}

struct ScreenCandidate {
  string survivor, epoch, attempt, globalAttempt;
  string word, r, u;
};

struct ScreenTarget {
  string tid, cls, block, coordinates, target;
  vector<ScreenCandidate> candidates;
};

struct RefinementResult {
  bool found = false;
  int candidatePosition = -1;
  int cells = 0;
  int colorCounts[3] = {0, 0, 0};
  bool strictBase = false;
  vector<int> singletonPoints;
};

static vector<string> splitWords(const string &line) {
  istringstream input(line);
  vector<string> words;
  string word;
  while (input >> word) words.push_back(word);
  return words;
}

static vector<ScreenTarget> loadScreen(const char *path) {
  ifstream input(path);
  if (!input) die("cannot open witness screen");
  vector<ScreenTarget> targets;
  string line;
  long long rows = 0;
  while (getline(input, line)) {
    if (line.empty() || line[0] == '#') continue;
    vector<string> words = splitWords(line);
    if (words.size() != 14) die("witness-screen row must have 14 fields");
    if (targets.empty() || targets.back().tid != words[0]) {
      ScreenTarget target;
      target.tid = words[0];
      target.cls = words[1];
      target.block = words[2];
      target.coordinates = words[3];
      target.target = words[4];
      targets.push_back(std::move(target));
    } else if (targets.back().cls != words[1] ||
               targets.back().block != words[2] ||
               targets.back().coordinates != words[3] ||
               targets.back().target != words[4]) {
      die("inconsistent witness-screen target rows");
    }
    ScreenCandidate candidate;
    candidate.survivor = words[5];
    candidate.epoch = words[6];
    candidate.attempt = words[7];
    candidate.globalAttempt = words[8];
    candidate.word = words[11];
    candidate.r = words[12];
    candidate.u = words[13];
    if ((int)candidate.r.size() != F4DIM ||
        (int)candidate.u.size() != F4DIM)
      die("witness-screen vector length");
    targets.back().candidates.push_back(std::move(candidate));
    rows++;
  }
  if (targets.empty()) die("witness screen is empty");
  fprintf(stderr, "witness screen: %zu targets, %lld candidates\n",
          targets.size(), rows);
  return targets;
}

static bool containsBase(const vector<int> &singletonPoints) {
  static thread_local vector<uint8_t> isSingleton(NPTS + 1);
  fill(isSingleton.begin(), isSingleton.end(), 0);
  for (int point : singletonPoints) isSingleton[point] = 1;
  for (int point : quotientBase)
    if (!isSingleton[point]) return false;
  return true;
}

static int threadCount() {
  const char *value = getenv("SLURM_CPUS_PER_TASK");
  int count = value ? atoi(value) : 0;
  if (count < 1) count = (int)thread::hardware_concurrency();
  return count > 0 ? count : 1;
}

static void refineAllCandidates(const char *screenPath,
                                const char *outputPath,
                                size_t limit, size_t firstTarget) {
  vector<ScreenTarget> allTargets = loadScreen(screenPath);
  if (firstTarget > allTargets.size()) die("first target is out of range");
  size_t end = allTargets.size();
  if (limit > 0 && firstTarget + limit < end) end = firstTarget + limit;
  vector<ScreenTarget> targets(
      make_move_iterator(allTargets.begin() + firstTarget),
      make_move_iterator(allTargets.begin() + end));
  vector<vector<RefinementResult>> results(targets.size());
  for (size_t targetPosition = 0; targetPosition < targets.size();
       targetPosition++)
    results[targetPosition].resize(targets[targetPosition].candidates.size());

  atomic<size_t> next{0};
  atomic<long long> failures{0}, strict{0}, candidatesDone{0};
  vector<thread> workers;
  int threads = threadCount();
  time_t started = time(nullptr);
  for (int threadPosition = 0; threadPosition < threads; threadPosition++) {
    workers.emplace_back([&]() {
      WorkBuffers *work = new WorkBuffers;
      for (;;) {
        size_t targetPosition = next.fetch_add(1);
        if (targetPosition >= targets.size()) break;
        const ScreenTarget &target = targets[targetPosition];
        for (size_t candidatePosition = 0;
             candidatePosition < target.candidates.size();
             candidatePosition++) {
          RefinementResult &result =
              results[targetPosition][candidatePosition];
          uint64_t hVector =
              h54FromF4(target.candidates[candidatePosition].u);
          initialColors(hVector, *work, result.colorCounts);
          result.cells = refine(*work);
          singletonCells(*work, result.cells, result.singletonPoints);
          result.found = !result.singletonPoints.empty();
          result.candidatePosition = (int)candidatePosition;
          result.strictBase = containsBase(result.singletonPoints);
          if (!result.found) failures++;
          if (result.strictBase) strict++;
          candidatesDone++;
        }
        size_t done = targetPosition + 1;
        if (done % 250 == 0)
          fprintf(stderr,
              "H-WL-all progress %zu/%zu candidates=%lld failures=%lld "
              "strict=%lld t=%llds\n",
              done, targets.size(), (long long)candidatesDone.load(),
              (long long)failures.load(), (long long)strict.load(),
              (long long)difftime(time(nullptr), started));
      }
      delete work;
    });
  }
  for (auto &worker : workers) worker.join();

  FILE *output = fopen(outputPath, "w");
  if (!output) die("cannot open all-candidate result output");
  fprintf(output, "# FI22_ALL_REFINED_CANDIDATES_V1\n");
  fprintf(output, "# first_target=%zu targets=%zu min_singletons=%d\n",
          firstTarget, targets.size(), MIN_SINGLETONS);
  fprintf(output, "# tid cls block coordinates target survivor epoch attempt "
                  "global_attempt word r u value0 value_orbit1 value_orbit2 "
                  "cells nsingle strict_base singleton_points\n");
  for (size_t targetPosition = 0; targetPosition < targets.size();
       targetPosition++) {
    const ScreenTarget &target = targets[targetPosition];
    for (size_t candidatePosition = 0;
         candidatePosition < target.candidates.size();
         candidatePosition++) {
      const ScreenCandidate &candidate =
          target.candidates[candidatePosition];
      const RefinementResult &result =
          results[targetPosition][candidatePosition];
      if (!result.found) {
        fprintf(output,
                "%s %s %s %s %s FAIL %s %s %s %s %s %s "
                "0 0 0 0 0 0 -\n",
                target.tid.c_str(), target.cls.c_str(),
                target.block.c_str(), target.coordinates.c_str(),
                target.target.c_str(), candidate.epoch.c_str(),
                candidate.attempt.c_str(), candidate.globalAttempt.c_str(),
                candidate.word.c_str(), candidate.r.c_str(),
                candidate.u.c_str());
        continue;
      }
      string singletonString;
      for (size_t position = 0;
           position < result.singletonPoints.size(); position++) {
        if (position) singletonString.push_back(',');
        singletonString += to_string(result.singletonPoints[position]);
      }
      fprintf(output,
              "%s %s %s %s %s %s %s %s %s %s %s %s "
              "%d %d %d %d %zu %d %s\n",
              target.tid.c_str(), target.cls.c_str(),
              target.block.c_str(), target.coordinates.c_str(),
              target.target.c_str(), candidate.survivor.c_str(),
              candidate.epoch.c_str(), candidate.attempt.c_str(),
              candidate.globalAttempt.c_str(), candidate.word.c_str(),
              candidate.r.c_str(), candidate.u.c_str(),
              result.colorCounts[0], result.colorCounts[1],
              result.colorCounts[2], result.cells,
              result.singletonPoints.size(), (int)result.strictBase,
              singletonString.c_str());
    }
  }
  fprintf(output, "# done first_target=%zu targets=%zu candidates=%lld "
                  "failures=%lld strict_base=%lld\n",
          firstTarget, targets.size(), (long long)candidatesDone.load(),
          (long long)failures.load(), (long long)strict.load());
  fclose(output);
  fprintf(stderr,
      "H-WL-all done: first=%zu targets=%zu candidates=%lld failures=%lld "
      "strict=%lld t=%llds\n",
      firstTarget, targets.size(), (long long)candidatesDone.load(),
      (long long)failures.load(), (long long)strict.load(),
      (long long)difftime(time(nullptr), started));
}

int main(int argc, char **argv) {
  initializePath();
  loadData();
  if (argc < 3 || argc > 5) {
    fprintf(stderr,
        "usage: %s <sum_candidates.tsv> <output.tsv> "
        "[limit] [first_target]\n", argv[0]);
    return 1;
  }
  size_t limit = argc >= 4 ? (size_t)strtoull(argv[3], nullptr, 10) : 0;
  size_t first = argc >= 5 ? (size_t)strtoull(argv[4], nullptr, 10) : 0;
  refineAllCandidates(argv[1], argv[2], limit, first);
  return 0;
}
