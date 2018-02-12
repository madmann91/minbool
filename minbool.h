#ifndef MINBOOL_H
#define MINBOOL_H

#include <cstdint>
#include <bitset>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <ostream>

namespace minbool {

template <size_t N>
bool increase(std::bitset<N>& bs)
{
    for (size_t i = 0; i != bs.size(); ++i) {
        if (bs.flip(i).test(i) == true)
            return true;
    }
    return false;
}

template<size_t N>
bool operator < (const std::bitset<N>& x, const std::bitset<N>& y) {
    for (int i = N-1; i >= 0; i--) {
        if (x[i] ^ y[i]) return y[i];
    }
    return false;
}

template <size_t Nbits>
struct MinTerm {
    struct Hash {
        size_t operator () (const MinTerm& term) const noexcept {
            std::hash<std::bitset<Nbits>> hash;
            return 33 * hash(term.value) ^ hash(term.dash);
        }
    };

    enum Value {
        Zero = 0,
        One  = 1,
        Dash = 2
    };

    MinTerm(std::bitset<Nbits> value, std::bitset<Nbits> dash = 0)
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

    std::unordered_map<std::bitset<Nbits>, std::vector<MinTermN>> columns;
    std::unordered_map<MinTermN, std::vector<std::bitset<Nbits>>, typename MinTermN::Hash> rows;

    PrimeChart(const std::vector<MinTermN>& primes) {
        for (auto& prime : primes) {
            prime.foreach_value([&] (std::bitset<Nbits> value) {
                columns[value].emplace_back(prime);
                rows[prime].emplace_back(value);
            });
        }
        for (auto& pair : rows)
            std::sort(pair.second.begin(), pair.second.end());
        for (auto& pair : columns)
            std::sort(pair.second.begin(), pair.second.end());
    }

    size_t size() const { return columns.size(); }

    void remove_row(const MinTermN& row) {
        rows.erase(row);
        // Remove columns that are covered by this row
        for (auto it = columns.begin(); it != columns.end(); ++it) {
            if (std::binary_search(it->second.begin(), it->second.end(), row))
                columns.erase(it);
        }
    }

    void remove_column(const std::bitset<Nbits>& column) {
        columns.erase(column);
        // If a row only covers this column, remove it
        for (auto it = rows.begin(); it != rows.end(); ++it) {
            it->second.erase(std::remove(it->second.begin(), it->second.end(), column), it->second.end());
            if (it->second.size() == 0)
                rows.erase(it);
        }
    }

    bool remove_essentials(std::vector<MinTermN>& essentials) {
        size_t count = essentials.size();
        for (auto& pair : columns) {
            if (pair.second.size() == 1)
                essentials.emplace_back(pair.second.front());
        }
        // No essential prime has been found
        if (essentials.size() == count)
            return false;
        std::sort(essentials.begin() + count, essentials.end());
        essentials.erase(std::unique(essentials.begin() + count, essentials.end()), essentials.end());
        for (auto& essential : essentials)
            remove_row(essential);
        return true;
    }

    void remove_heuristic(std::vector<MinTermN>& solution) {
        assert(size() > 0);
        // Heuristic: Remove the term that covers the most columns
        size_t max_covers = 0;
        MinTermN term(0, 0);
        for (auto& pair : rows) {
            if (pair.second.size() > max_covers) {
                max_covers = pair.second.size();
                term = pair.first;
            }
        }
        solution.emplace_back(term);
        remove_row(term);
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
                    remove_column(pair2.first);
                    change = true;
                    break;
                }
            }
        }

        for (auto& pair1 : rows) {
            for (auto& pair2 : rows) {
                if (pair1.first == pair2.first)
                    continue;
                // Dominated rows are eliminated
                if (std::includes(pair1.second.begin(), pair1.second.end(),
                                  pair2.second.begin(), pair2.second.end())) {
                    remove_row(pair2.first);
                    change = true;
                    break;
                }
            }
        }

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
                    const std::vector<std::bitset<Nbits>>& on_values,
                    const std::vector<std::bitset<Nbits>>& dc_values) {
    std::unordered_set<std::bitset<Nbits>> not_off;
    for (auto v : on_values) not_off.emplace(v);
    for (auto v : dc_values) not_off.emplace(v);
    for (auto v : on_values) {
        if (!eval_boolean(solution, v))
            return false;
    }
    auto i = std::bitset<Nbits>(0);
    do {
        if (not_off.count(i) == 0 && eval_boolean(solution, i))
            return false;
    } while (increase(i));
    return true;
}

template <size_t Nbits>
std::vector<MinTerm<Nbits>> minimize_boolean(
    const std::vector<std::bitset<Nbits>>& on_values,
    const std::vector<std::bitset<Nbits>>& dc_values)
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

    PrimeChart<Nbits> chart(primes);
    for (auto dc : dc_values)
        chart.remove_column(dc);

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
