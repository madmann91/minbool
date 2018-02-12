#include <unordered_set>
#include <random>
#include <chrono>
#include <iostream>
#include <fstream>

#include "minbool.h"

int main(int /*argc*/, char** /*argv*/) {
    using namespace std::chrono;
    using namespace minbool;

    auto seed = 1518454078902;//duration_cast<milliseconds>(high_resolution_clock::now().time_since_epoch()).count();
    std::mt19937 gen(seed);
    auto rand1024 = std::uniform_int_distribution<size_t>(0, 1024);
    auto rand8192 = std::uniform_int_distribution<size_t>(0, 8192);
    auto rand65535 = std::uniform_int_distribution<size_t>(0, 65535);

    std::unordered_set<uint16_t> on_set, dc_set;
    for (size_t i = 0, n = rand8192(gen); i < n; ++i)
        on_set.emplace(rand65535(gen));
    for (size_t i = 0, n = rand1024(gen); i < n; ++i) {
        auto value = rand65535(gen);
        if (!on_set.count(value))
            dc_set.emplace(value);
    }
    auto start = high_resolution_clock::now();
    std::vector<uint16_t> on(on_set.begin(), on_set.end());
    std::vector<uint16_t> dc(dc_set.begin(), dc_set.end());
    auto solution = minimize_boolean<16>(on, dc);
    auto end = high_resolution_clock::now();

    std::cout << "seed " << seed << std::endl;
    std::cout << duration_cast<milliseconds>(end - start).count() << " ms" << std::endl;
    std::cout << solution.size() << " terms" << std::endl;
    for (auto& term : solution)
        std::cout << term << std::endl;

    std::ofstream of("test.esp");
    of << ".i 16\n";
    of << ".o 1\n";
    for (auto i : on) of << std::bitset<16>(i) << " 1\n";
    for (auto i : dc) of << std::bitset<16>(i) << " -\n";
    of << ".e\n";
    return 0;
}
