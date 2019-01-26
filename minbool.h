#ifndef MINBOOL_H
#define MINBOOL_H

#include <cstdint>
#include <bitset>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cassert>
#include <ostream>
#include <numeric>

namespace minbool {

inline size_t popcount(uint8_t  n) { return std::bitset<8>(n).count();  }
inline size_t popcount(uint16_t n) { return std::bitset<16>(n).count(); }
inline size_t popcount(uint32_t n) { return std::bitset<32>(n).count(); }
inline size_t popcount(uint64_t n) { return std::bitset<64>(n).count(); }

template <size_t Nbits> struct IntType {};
template <> struct IntType<8>  { using Type = uint8_t;  };
template <> struct IntType<16> { using Type = uint16_t; };
template <> struct IntType<32> { using Type = uint32_t; };
template <> struct IntType<64> { using Type = uint64_t; };

constexpr size_t compute_int_size(size_t Nbits) {
    return Nbits <= 8 ? 8 : 2 * compute_int_size(Nbits / 2);
}

template <size_t Nbits>
struct MinTerm {
    using IntTypeN = typename IntType<compute_int_size(Nbits)>::Type;

    struct Hash {
        size_t operator () (const MinTerm& term) const {
            // Bernstein's hash function
            return 33 * term.value ^ term.dash;
        }
    };

    enum Value {
        Zero = 0,
        One  = 1,
        Dash = 2
    };

    MinTerm(IntTypeN value = 0, IntTypeN dash = 0)
        : value(value), dash(dash)
    {
        assert((value & ~((1 << Nbits) - 1)) == 0);
        assert((dash  & ~((1 << Nbits) - 1)) == 0);
    }

    Value operator [] (size_t i) const {
        assert(i < Nbits);
        return (dash >> i) & 1 ? Dash : Value((value >> i) & 1);
    }

    size_t count_literals() const {
        return Nbits - popcount(dash);
    }

    MinTerm combine(const MinTerm& other) const {
        IntTypeN mask = (value ^ other.value) | (dash ^ other.dash);
        return MinTerm(value & ~mask, dash | mask);
    }

    template <typename F>
    void foreach_value(F f, size_t bit = 0, IntTypeN cur = 0) const {
        if (bit == Nbits) {
            f(cur);
        } else {
            auto value = (*this)[bit];
            if (value == Dash) {
                foreach_value(f, bit + 1, cur);
                foreach_value(f, bit + 1, cur | (1 << bit));
            } else {
                foreach_value(f, bit + 1, cur | (value == Zero ? 0 : (1 << bit)));
            }
        }
    }

    bool operator < (const MinTerm& other) const {
        return value < other.value || (value == other.value && dash < other.dash);
    }

    bool operator == (const MinTerm& other) const {
        return value == other.value && dash == other.dash;
    }

    IntTypeN value;
    IntTypeN dash;
};

template <size_t Nbits>
std::ostream& operator << (std::ostream& os, const MinTerm<Nbits>& term) {
    for (int i = Nbits - 1; i >= 0; --i) {
        auto value = term[i];
        if (value != MinTerm<Nbits>::Dash)
            os << (uint32_t)value;
        else
            os << '-';
    }
    return os;
}

template <size_t Nbits>
struct ImplicantTable {
    using MinTermN = MinTerm<Nbits>;

    size_t groups[Nbits + 2];
    std::vector<bool>     marks;
    std::vector<MinTermN> terms;

    size_t size() const { return terms.size(); }

    void fill(const std::vector<MinTermN>& minterms) {
        std::fill(groups, groups + Nbits + 2, 0);
        for (auto& term : minterms)
            groups[popcount(term.value)]++;
        std::partial_sum(groups, groups + Nbits + 2, groups);
        terms.resize(minterms.size());
        marks.resize(minterms.size());
        for (auto& term : minterms)
            terms[--groups[popcount(term.value)]] = term;
    }

    void combine(std::vector<MinTermN>& res) {
        for (size_t i = 0; i < Nbits; ++i) {
            for (size_t j = groups[i]; j < groups[i + 1]; ++j) {
                for (size_t k = groups[i + 1]; k < groups[i + 2]; ++k) {
                    auto& term_a = terms[j];
                    auto& term_b = terms[k];
                    if ((term_a.value & term_b.value) == term_a.value && (term_a.dash == term_b.dash)) {
                        marks[j] = true;
                        marks[k] = true;
                        res.push_back(term_a.combine(term_b));
                    }
                }
            }
        }
    }

    void primes(std::vector<MinTermN>& res) {
        for (size_t i = 0; i < terms.size(); ++i) {
            if (!marks[i])
                res.push_back(terms[i]);
        }
    }
};

template <size_t Nbits>
std::ostream& operator << (std::ostream& os, const ImplicantTable<Nbits>& table) {
    for (size_t i = 0; i < Nbits; ++i) {
        if (!table.terms[i].empty()) {
            for (size_t j = 0; j < table.terms[i].size(); ++j) {
                os << table.terms[i][j] << ' ' << (table.marks[i][j] ? 'X' : ' ') << '\n';
            }
            for (size_t j = 0; j < Nbits + 2; ++j)
                os << '-';
            os << '\n';
        }
    }
    return os;
}

template <size_t Nbits>
struct PrimeChart {
    using MinTermN = MinTerm<Nbits>;
    using IntTypeN = typename MinTermN::IntTypeN;

    std::unordered_map<IntTypeN, std::vector<MinTermN>> columns;

    size_t size() const {
        return columns.size();
    }

    void fill(const std::vector<MinTermN>& primes) {
        for (auto& prime : primes) {
            prime.foreach_value([&] (IntTypeN value) {
                columns[value].emplace_back(prime);
            });
        }
        for (auto& pair : columns)
            std::sort(pair.second.begin(), pair.second.end());
    }

    void remove_columns(const std::vector<IntTypeN>& values) {
        for (auto value : values)
            columns.erase(value);
    }

    bool remove_essentials(std::vector<MinTermN>& essentials) {
        size_t count = essentials.size();
        for (auto& pair : columns) {
            if (pair.second.size() == 1)
                essentials.push_back(pair.second.front());
        }
        // No essential prime has been found
        if (essentials.size() == count)
            return false;
        std::sort(essentials.begin() + count, essentials.end());
        essentials.erase(std::unique(essentials.begin() + count, essentials.end()), essentials.end());

        std::for_each(essentials.begin() + count, essentials.end(), [&] (const MinTermN& term) {
            term.foreach_value([&] (IntTypeN value) {
                columns.erase(value);
            });
        });
        return true;
    }

    void remove_heuristic(std::vector<MinTermN>& solution) {
        assert(size() > 0);
        std::unordered_map<MinTermN, size_t, typename MinTermN::Hash> covers;
        for (auto& pair : columns) {
            for (auto& term : pair.second)
                covers[term]++;
        }
        // Heuristic: Remove the term that covers the most columns
        size_t max_covers = 0;
        MinTermN term;
        for (auto& pair : covers) {
            if (pair.second > max_covers) {
                max_covers = pair.second;
                term = pair.first;
            }
        }
        solution.emplace_back(term);
        term.foreach_value([&] (IntTypeN value) {
            columns.erase(value);
        });
    }

    bool simplify() {
        bool change = false;

        for (auto& pair1 : columns) {
            for (auto& pair2 : columns) {
                if (pair1.first == pair2.first)
                    continue;
                // Dominating columns are eliminated
                if (std::includes(pair2.second.begin(), pair2.second.end(),
                                  pair1.second.begin(), pair1.second.end())) {
                    columns.erase(pair2.first);
                    change = true;
                    break;
                }
            }
        }

        // Transpose columns => rows
        std::unordered_map<MinTermN, std::vector<IntTypeN>, typename MinTermN::Hash> rows;
        for (auto& pair : columns) {
            for (auto& term : pair.second)
                rows[term].emplace_back(pair.first);
            pair.second.clear();
        }
        for (auto& pair : rows)
            std::sort(pair.second.begin(), pair.second.end());

        for (auto& pair1 : rows) {
            for (auto& pair2 : rows) {
                if (pair1.first == pair2.first)
                    continue;
                // Dominated rows are eliminated
                if (std::includes(pair1.second.begin(), pair1.second.end(),
                                  pair2.second.begin(), pair2.second.end())) {
                    rows.erase(pair2.first);
                    change = true;
                    break;
                }
            }
        }

        // Transpose rows => columns
        for (auto& pair : rows) {
            for (auto& value : pair.second)
                columns[value].emplace_back(pair.first);
        }
        for (auto& pair : columns)
            std::sort(pair.second.begin(), pair.second.end());

        return change;
    }
};

template <size_t Nbits>
std::vector<MinTerm<Nbits>> prime_implicants(std::vector<MinTerm<Nbits>>& terms) {
    std::vector<MinTerm<Nbits>> primes;

    while (!terms.empty()) {
        ImplicantTable<Nbits> table;
        table.fill(terms);
        terms.clear();
        table.combine(terms);
        // Remove duplicates
        std::sort(terms.begin(), terms.end());
        terms.erase(std::unique(terms.begin(), terms.end()), terms.end());
        table.primes(primes);
    }
    return primes;
}

template <size_t Nbits>
bool eval_boolean(const std::vector<MinTerm<Nbits>>& solution, typename MinTerm<Nbits>::IntTypeN v) {
    for (auto& term : solution) {
        bool prod = true;
        for (size_t i = 0; i < Nbits; ++i) {
            bool bit = ((v >> i) & 1) ? true : false;
            if (term[i] == MinTerm<Nbits>::One)
                prod &= bit;
            else if (term[i] == MinTerm<Nbits>::Zero)
                prod &= !bit;
        }
        if (prod) return true;
    }
    return false;
}

template <size_t Nbits>
bool check_solution(const std::vector<MinTerm<Nbits>>& solution,
                    const std::vector<typename MinTerm<Nbits>::IntTypeN>& on_values,
                    const std::vector<typename MinTerm<Nbits>::IntTypeN>& dc_values) {
    using IntTypeN = typename MinTerm<Nbits>::IntTypeN;
    std::unordered_set<IntTypeN> not_off;
    for (auto v : on_values) not_off.emplace(v);
    for (auto v : dc_values) not_off.emplace(v);
    for (auto v : on_values) {
        if (!eval_boolean(solution, v))
            return false;
    }
    for (IntTypeN i = (1 << Nbits) - 1; i > 0; --i) {
        if (not_off.count(i) == 0 && eval_boolean(solution, i))
            return false;
    }
    return not_off.count(0) != 0 || !eval_boolean(solution, 0);
}

template <size_t Nbits>
std::vector<MinTerm<Nbits>> minimize_boolean(
    const std::vector<typename MinTerm<Nbits>::IntTypeN>& on_values,
    const std::vector<typename MinTerm<Nbits>::IntTypeN>& dc_values)
{
    if (on_values.size() <= 1) {
        return on_values.size() == 1
            ? std::vector<MinTerm<Nbits>>{ on_values.front() }
            : std::vector<MinTerm<Nbits>>();
    }

    std::vector<MinTerm<Nbits>> init;
    init.reserve(on_values.size() + dc_values.size());
    for (auto on : on_values)
        init.emplace_back(on);
    for (auto dc : dc_values)
        init.emplace_back(dc);
    auto primes = prime_implicants(init);

    PrimeChart<Nbits> chart;
    chart.fill(primes);
    chart.remove_columns(dc_values);

    std::vector<MinTerm<Nbits>> solution;
    do {
        bool change = chart.remove_essentials(solution);
        change |= chart.simplify();
        if (!change && chart.size() > 0)
            chart.remove_heuristic(solution);
    } while (chart.size() > 0);

    assert(check_solution(solution, on_values, dc_values));
    return solution;
}

} // namespace minbool

#endif // MINBOOL_H
