#ifndef MINBOOL_H
#define MINBOOL_H

#include <cstdint>
#include <bitset>
#include <vector>
#include <algorithm>
#include <cassert>
#include <ostream>

namespace minbool {

template <size_t Nbits>
struct MinTerm {
    enum Value {
        Zero = 0,
        One  = 1,
        Dash = 2
    };

    MinTerm(std::bitset<Nbits> value, std::bitset<Nbits> dash)
        : value(value), dash(dash)
    {}

    Value operator [] (size_t i) const {
        return dash[i] ? Dash : (value[i] ? One : Zero);
    }

    MinTerm combine(const MinTerm& other) const {
        auto mask = (value ^ other.value) | (dash ^ other.dash);
        return MinTerm(value & ~mask, dash | mask);
    }

    template <typename F>
    void foreach_value(F f, size_t bit = 0, std::bitset<Nbits> cur = 0) const {
        if (bit == Nbits) {
            f(cur);
        } else {
            auto value = (*this)[bit];
            if (value == Dash) {
                foreach_value(f, bit + 1, cur);
                foreach_value(f, bit + 1, cur.set(bit));
            } else {
                foreach_value(f, bit + 1, cur.set(bit, value == One));
            }
        }
    }

    bool operator < (const MinTerm& other) const {
        return value < other.value || (value == other.value && dash < other.dash);
    }

    bool operator == (const MinTerm& other) const {
        return value == other.value && dash == other.dash;
    }

    std::bitset<Nbits> value;
    std::bitset<Nbits> dash;
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

    size_t groups[Nbits + 1];
    std::vector<MinTermN> terms;

    ImplicantTable(const std::vector<MinTermN>& init) {
        for (auto& term : init)
            groups[term.value.count()]++;
        std::partial_sum(groups, groups + Nbits + 1, groups);
        terms.resize(init.size());
        for (auto& term : init)
            terms[--groups[term.value.count()]] = term;
    }

    void combine(std::vector<MinTermN>& res, std::vector<bool>& marks) {
        marks.resize(terms.size());
        std::fill(marks.begin(), marks.end(), false);

        for (size_t i = 0; i < Nbits - 1; ++i) {
            for (size_t j = groups[i]; j < groups[i + 1]; ++j) {
                for (size_t k = groups[i + 1]; k < groups[i + 2], ++k) {
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

    void primes(std::vector<MinTermN>& res, const std::vector<bool>& marks) {
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

    std::vector<std::bitset<Nbits>> columns;
    std::vector<MinTermN> rows;

    PrimeChart(const std::vector<MinTermN>& primes)
        : rows(primes)
    {
        for (auto& prime : primes) {
            prime.foreach_value([&] (std::bitset<Nbits> value) {
                columns.emplace_back(value);
            });
        }
        std::sort(columns.begin(), columns.end());
        columns.erase(std::unique(columns.begin(), columns.end()), columns.end());7
    }

    bool remove_essentials(std::vector<MinTermN>& essentials, size_t ) {
        std::vector<size_t> covers;
        size_t count = essentials.size();
        for (auto

        // No essential prime has been found
        if (essentials.size() == count)
            return false;
        std::sort(essentials.begin() + count, essentials.end());
        essentials.erase(std::unique(essentials.begin() + count, essentials.end()), essentials.end());
        std::remove_if(elements.begin(), elements.end(), [] (const Element& elem) {
            return 
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

    std::vector<bool> marks;
    while (!terms.empty()) {
        ImplicantTable<Nbits> table(terms);
        terms.clear();
        table.combine(terms, marks);
        // Remove duplicates
        std::sort(terms.begin(), terms.end());
        terms.erase(std::unique(terms.begin(), terms.end()), terms.end());
        table.primes(primes, marks);
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
