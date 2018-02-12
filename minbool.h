#ifndef MINBOOL_H
#define MINBOOL_H

#include <cstdint>
#include <vector>
#include <array>
#include <bitset>
#include <limits>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <ostream>

namespace minbool {

// Same as std::bitset but supports iteration and comparison
template <size_t Nbits>
struct BitSet {
    struct Hash {
        size_t operator () (const BitSet& bitset) const noexcept {
            return bitset.hash();
        }
    };

    std::array<uint64_t, Nbits / 64 + (Nbits % 64 ? 1 : 0)> bits;

    BitSet() { reset(); }
    BitSet(std::bitset<Nbits> value) {
        for (size_t i = 0; i < Nbits; ++i) {
            if (value.test(i))
                set(i);
        }
    }
    BitSet(uint64_t value)
        : bits { value }
    {}

    size_t count() const {
        size_t sum = 0;
        for (auto& word : bits)
            sum += std::bitset<64>(word).count();
        return sum;
    }

    BitSet<Nbits>& set(size_t i, bool value = true) {
        uint64_t mask = 1 << (i % 64);
        if (value)
            bits[i / 64] |= mask;
        else
            bits[i / 64] &= ~mask;
        return *this;
    }

    BitSet<Nbits>& flip(size_t i) {
        bits[i / 64] ^= 1 << (i % 64);
        return *this;
    }

    BitSet<Nbits>& reset() {
        for (auto& word : bits)
            word = 0;
        return *this;
    }

    BitSet<Nbits>& reset(size_t i) {
        set(i, false);
    }

    bool operator [] (size_t i) const {
        return (bits[i / 64] & (1 << (i % 64))) != 0;
    }

    size_t hash() const {
        size_t h = 0;
        for (auto word : bits)
            h = h * 33 ^ word;
        return h;
    }

    bool next() {
        static constexpr uint64_t max_word = Nbits < 64
            ? (1 << Nbits) - 1
            : std::numeric_limits<uint64_t>::max();
        for (auto& word : bits) {
            if (word != max_word) {
                word++;
                return true;
            }
            word = 0;
        }
        return false;
    }

    bool operator == (const BitSet<Nbits>& other) const {
        return bits == other.bits;
    }

    bool operator < (const BitSet<Nbits>& other) const {
        return std::lexicographical_compare(bits.begin(), bits.end(), other.bits.begin(), other.bits.end());
    }
    
    template <typename Op>
    BitSet<Nbits> map(Op op) const {
        BitSet<Nbits> res;
        for (size_t i = 0; i < bits.size(); ++i)
            res.bits[i] = op(bits[i]);
        return res;
    }

    template <typename Op>
    BitSet<Nbits> zip(const BitSet<Nbits>& other, Op op) const {
        BitSet<Nbits> res;
        for (size_t i = 0; i < bits.size(); ++i)
            res.bits[i] = op(bits[i], other.bits[i]);
        return res;
    }

    BitSet<Nbits> operator & (const BitSet<Nbits>& other) const { return zip(other, [] (uint64_t a, uint64_t b) { return a & b; }); }
    BitSet<Nbits> operator | (const BitSet<Nbits>& other) const { return zip(other, [] (uint64_t a, uint64_t b) { return a | b; }); }
    BitSet<Nbits> operator ^ (const BitSet<Nbits>& other) const { return zip(other, [] (uint64_t a, uint64_t b) { return a ^ b; }); }
    BitSet<Nbits> operator ~ () const { return map([] (uint64_t a) { return ~a; }); }
};

template <size_t Nbits>
std::ostream& operator << (std::ostream& os, const BitSet<Nbits>& bitset) {
    for (size_t i = 0; i < Nbits; ++i)
        os << (bitset[i] ? '1' : '0');
    return os;
}

template <size_t Nbits>
struct MinTerm {
    struct Hash {
        size_t operator () (const MinTerm& term) const noexcept {
            return 33 * term.value.hash() ^ term.dash.hash();
        }
    };

    enum Value {
        Zero = 0,
        One  = 1,
        Dash = 2
    };

    MinTerm() {}

    MinTerm(BitSet<Nbits> value, BitSet<Nbits> dash = 0)
        : value(value), dash(dash)
    {}

    Value operator [] (size_t i) const {
        return dash[i] ? Dash : (value[i] ? One : Zero);
    }

    bool cover(const BitSet<Nbits>& bitset) {
        return (value ^ bitset) & ~dash == BitSet<Nbits>();
    }

    MinTerm combine(const MinTerm& other) const {
        auto mask = (value ^ other.value) | (dash ^ other.dash);
        return MinTerm(value & ~mask, dash | mask);
    }

    template <typename F>
    void foreach_value(F f, size_t bit = 0, BitSet<Nbits> cur = 0) const {
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

    BitSet<Nbits> value;
    BitSet<Nbits> dash;
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
        std::fill(groups, groups + Nbits + 1, 0);
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
    using BitSetN  = BitSet<Nbits>;

    std::unordered_map<BitSetN, std::vector<MinTermN>, typename BitSetN::Hash> columns;
    std::unordered_map<MinTermN, std::vector<BitSetN>, typename MinTermN::Hash> rows;

    PrimeChart(const std::vector<MinTermN>& primes) {
        for (auto& prime : primes) {
            prime.foreach_value([&] (BitSetN value) {
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
        for (auto it = columns.begin(); it != columns.end();) {
            it->second.erase(std::remove(it->second.begin(), it->second.end(), row), it->second.end());
            if (it->second.empty())
                it = columns.erase(it);
            else
                ++it;
        }
    }

    void remove_column(const BitSetN& column) {
        columns.erase(column);
        for (auto it = rows.begin(); it != rows.end();) {
            it->second.erase(std::remove(it->second.begin(), it->second.end(), column), it->second.end());
            if (it->second.empty())
                it = rows.erase(it);
            else
                ++it;
        }
    }

    void eliminate_prime(const MinTermN& prime) {
        auto it = rows.find(prime);
        assert(it != rows.end());
        for (auto& value : it->second)
            columns.erase(value);
        rows.erase(it);
    }

    bool eliminate_essentials(std::vector<MinTermN>& essentials) {
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
        for (size_t i = count; i < essentials.size(); ++i)
            eliminate_prime(essentials[i]);
        return true;
    }

    void eliminate_with_heuristic(std::vector<MinTermN>& solution) {
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
        eliminate_prime(term);
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
bool eval_boolean(const std::vector<MinTerm<Nbits>>& solution, const BitSet<Nbits>& v) {
    for (auto& term : solution) {
        bool prod = true;
        for (size_t i = 0; i < Nbits; ++i) {
            bool bit = v[i];
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
                    const std::vector<BitSet<Nbits>>& on_values,
                    const std::vector<BitSet<Nbits>>& dc_values) {
    std::unordered_set<BitSet<Nbits>, typename BitSet<Nbits>::Hash> not_off;
    for (auto v : on_values) not_off.emplace(v);
    for (auto v : dc_values) not_off.emplace(v);
    for (auto v : on_values) {
        if (!eval_boolean(solution, v))
            return false;
    }
    auto i = BitSet<Nbits>(0);
    do {
        if (not_off.count(i) == 0 && eval_boolean(solution, i))
            return false;
    } while (i.next());
    return true;
}

template <size_t Nbits>
std::vector<MinTerm<Nbits>> minimize_boolean(
    const std::vector<BitSet<Nbits>>& on_values,
    const std::vector<BitSet<Nbits>>& dc_values)
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
        bool change = chart.eliminate_essentials(solution);
        change |= chart.simplify();
        if (!change && chart.size() > 0)
            chart.eliminate_with_heuristic(solution);
    } while (chart.size() > 0);

    assert(check_solution(solution, on_values, dc_values));
    return solution;
}

} // namespace minbool

#endif // MINBOOL_H
