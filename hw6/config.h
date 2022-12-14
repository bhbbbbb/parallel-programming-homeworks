#ifndef CONFIG

#define CONFIG

using WeightType = unsigned long;

constexpr int DEFAULT_NUM_THREADS = 2;
constexpr int ALPHA = 1;
constexpr int BETA = 1;
constexpr int LAMBDA = 1;
constexpr int NUM_ANTS = 1 << 7;
constexpr int MU = 8;
constexpr int NUM_COLONIES = 1 << 9;
constexpr int QPOW = (sizeof(int) * 8 + 2);
constexpr WeightType Q = (1UL << QPOW);
constexpr WeightType QQ = Q >> 8;

constexpr int GLOBAL_EXANGE_INTERVAL = NUM_COLONIES >> 4;
constexpr bool PRINT_BEST_IN_THREADS = false;
constexpr bool PRINT_PHEROROME = false;

#endif
