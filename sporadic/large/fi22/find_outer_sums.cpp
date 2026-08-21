// Deterministic CREATE-side witness search for the 743 outer-involution
// fixed-space targets of 3.Fi22.2.
//
// Reuse the reference arithmetic, dual-orbit screen, and random-walk machinery
// from find_inner_sums.cpp.  The only changed input is the outer target
// table.  Target ids are the global ids 5363..6105, so this output can be
// concatenated directly with the existing 5363-row-family inner screen.
//
// Usage:
//   find_outer_sums <root> <outfile> <threads> <max-attempts> <survivors>

#define main fi22_inner_witness_search_main
#include "find_inner_sums.cpp"
#undef main

static constexpr size_t OUTER_TARGETS = 743;
static constexpr size_t GLOBAL_TID_START = 5363;

static vector<string> splitFields(const string &line) {
  vector<string> fields;
  size_t position = 0;
  while (position < line.size()) {
    while (position < line.size() &&
           isspace(static_cast<unsigned char>(line[position])))
      ++position;
    if (position == line.size()) break;
    size_t end = position;
    while (end < line.size() &&
           !isspace(static_cast<unsigned char>(line[end])))
      ++end;
    fields.push_back(line.substr(position, end - position));
    position = end;
  }
  return fields;
}

static bool hasOnlyDigits(const string &text, char first, char last) {
  return !text.empty() &&
      all_of(text.begin(), text.end(), [=](char character) {
        return first <= character && character <= last;
      });
}

static vector<Target> collectOuterTargets(const string &root) {
  const string path =
      root + "/generated/outer_fixed/outer_targets_wl.tsv";
  FILE *file = fopen(path.c_str(), "r");
  if (!file) die("cannot open outer target table " + path);

  vector<Target> targets;
  unordered_set<string> seen;
  for (;;) {
    string line = readLine(file);
    if (line.empty() && feof(file)) break;
    if (line.rfind("TARGET ", 0) != 0) continue;
    vector<string> fields = splitFields(line);
    if (fields.size() != 8 || fields[0] != "TARGET")
      die("malformed outer target row");

    char *end = nullptr;
    unsigned long localId = strtoul(fields[1].c_str(), &end, 10);
    if (!end || *end != '\0' || localId != targets.size())
      die("outer target ids are not consecutive");
    if (fields[2].size() != 4 ||
        fields[2].rfind("OUT", 0) != 0 ||
        fields[2][3] < '1' || fields[2][3] > '3')
      die("invalid outer family");
    if (!hasOnlyDigits(fields[5], '0', '1') ||
        fields[5].size() != 54)
      die("invalid 54-dimensional outer representative");
    if (!hasOnlyDigits(fields[6], '0', '3') ||
        fields[6].size() != DIM)
      die("invalid GF(4) outer target");
    if (fields[7] != "0")
      die("outer target overlaps an inner target");

    Target target;
    target.cls = fields[2];
    target.block = fields[2][3] - '0';
    // Keep the exact 54-bit representative as the origin coordinate.
    target.coordinates = fields[5];
    target.digits = fields[6];
    decodeDigits(target.digits, target.vector, DIM);
    if (target.digits.find_first_not_of('0') == string::npos)
      die("zero outer target");
    if (!seen.insert(target.digits).second)
      die("duplicate outer target");
    targets.push_back(target);
  }
  fclose(file);

  if (targets.size() != OUTER_TARGETS)
    die("outer target count is " + to_string(targets.size()) +
        ", expected " + to_string(OUTER_TARGETS));
  return targets;
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
  vector<Target> targets = collectOuterTargets(root);
  fprintf(stderr,
      "loaded %zu outer targets; global tids=%zu..%zu; "
      "reference seed norm=%d zeros=%d\n",
      targets.size(), GLOBAL_TID_START,
      GLOBAL_TID_START + targets.size() - 1, referenceNorm, referenceZeros);

  vector<SearchResult> results(targets.size());
  atomic<size_t> next{0};
  atomic<size_t> completed{0};
  mutex logMutex;
  vector<thread> threads;
  for (int threadIndex = 0; threadIndex < threadCount; ++threadIndex) {
    threads.emplace_back([&]() {
      ZeroWorkspace workspace;
      for (;;) {
        size_t localIndex = next.fetch_add(1);
        if (localIndex >= targets.size()) break;
        // The random stream is keyed by the final global target id.
        results[localIndex] = searchTarget(
            targets[localIndex], GLOBAL_TID_START + localIndex,
            maxAttempts, want, workspace);
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
  fprintf(output, "# FI22_OUTER_WITNESS_SCREEN_V1\n");
  fprintf(output,
      "# tid cls block coordinates target survivor epoch attempt "
      "global_attempt norm zeros word r u\n");
  size_t failures = 0, survivorCount = 0;
  long long totalAttempts = 0;
  for (size_t localIndex = 0; localIndex < targets.size(); ++localIndex) {
    const Target &target = targets[localIndex];
    const SearchResult &result = results[localIndex];
    totalAttempts += result.attempts;
    if (static_cast<int>(result.survivors.size()) != want) ++failures;
    for (size_t survivorIndex = 0;
         survivorIndex < result.survivors.size(); ++survivorIndex) {
      const Survivor &survivor = result.survivors[survivorIndex];
      fprintf(output,
          "%zu %s %d %s %s %zu %d %d %d %d %d %s %s %s\n",
          GLOBAL_TID_START + localIndex, target.cls.c_str(), target.block,
          target.coordinates.c_str(), target.digits.c_str(), survivorIndex,
          survivor.epoch, survivor.attempt, survivor.globalAttempt,
          survivor.norm, survivor.zeros, survivor.word.c_str(),
          survivor.r.c_str(), survivor.u.c_str());
      ++survivorCount;
    }
  }
  fprintf(output,
      "# done targets=%zu survivors=%zu failures=%zu attempts=%lld "
      "global_tid_start=%zu global_tid_end=%zu\n",
      targets.size(), survivorCount, failures, totalAttempts,
      GLOBAL_TID_START, GLOBAL_TID_START + targets.size() - 1);
  fclose(output);
  fprintf(stderr,
      "DONE targets=%zu survivors=%zu failures=%zu attempts=%lld "
      "global_tids=%zu..%zu\n",
      targets.size(), survivorCount, failures, totalAttempts,
      GLOBAL_TID_START, GLOBAL_TID_START + targets.size() - 1);
  return failures == 0 ? 0 : 2;
}
