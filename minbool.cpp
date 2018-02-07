#include <cassert>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <chrono>
#include <iostream>

inline size_t popcount(uint8_t  n) { return __builtin_popcount(n);   }
inline size_t popcount(uint16_t n) { return __builtin_popcount(n);   }
inline size_t popcount(uint32_t n) { return __builtin_popcount(n);   }
inline size_t popcount(uint64_t n) { return __builtin_popcountll(n); }

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

    std::vector<bool>     marks[Nbits + 1];
    std::vector<MinTermN> terms[Nbits + 1];

    size_t size() const {
        size_t sum = 0;
        for (size_t i = 0; i <= Nbits; ++i)
            sum += terms[i].size();
        return sum;
    }

    void fill(const std::vector<MinTermN>& minterms) {
        for (auto& term : minterms) {
            auto i = popcount(term.value);
            terms[i].push_back(term);
            marks[i].push_back(false);
        }
    }

    void combine(std::vector<MinTermN>& res) {
        for (size_t i = 0; i <= Nbits - 1; ++i) {
            for (size_t j = 0; j < terms[i].size(); ++j) {
                for (size_t k = 0; k < terms[i + 1].size(); ++k) {
                    auto& term_a = terms[i][j];
                    auto& term_b = terms[i + 1][k];
                    if ((term_a.value & term_b.value) == term_a.value && (term_a.dash == term_b.dash)) {
                        marks[i][j] = true;
                        marks[i + 1][k] = true;
                        res.push_back(term_a.combine(term_b));
                    }
                }
            }
        }
    }

    void primes(std::vector<MinTermN>& res) {
        for (size_t i = 0; i <= Nbits; ++i) {
            for (size_t j = 0; j < terms[i].size(); ++j) {
                if (!marks[i][j])
                    res.push_back(terms[i][j]);
            }
        }
    }
};

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
    }

    void remove_columns(const std::vector<IntTypeN>& values) {
        for (auto value : values)
            columns.erase(value);
    }

    void remove_essentials(std::vector<MinTermN>& res) {
        for (auto& pair : columns) {
            if (pair.second.size() == 1)
                res.push_back(pair.second[0]);
        }
        std::sort(res.begin(), res.end());
        res.erase(std::unique(res.begin(), res.end()), res.end());
        for (auto& essential : res) {
            essential.foreach_value([&] (IntTypeN value) {
                columns.erase(value);
            });
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
struct Prod {
    using MinTermN = MinTerm<Nbits>;

    std::vector<MinTermN> terms;

    size_t size() const { return terms.size(); }

    size_t count_literals() const {
        size_t literals = 0;
        for (auto& term : terms)
            literals += term.count_literals();
        return literals;
    }

    bool contains(const Prod& other) const {
        for (auto& term : other.terms) {
            if (std::find(terms.begin(), terms.end(), term) == terms.end())
                return false;
        }
        return true;
    }

    void mul(const MinTermN& term) {
        // XX = X
        if (std::find(terms.begin(), terms.end(), term) != terms.end())
            return;
        terms.push_back(term);
    }

    void mul(const Prod& other) {
        for (auto& term : other.terms)
            mul(term);
    }
};

template <size_t Nbits>
struct Sum {
    using ProdN = Prod<Nbits>;

    std::vector<ProdN> prods;

    size_t size() const { return prods.size(); }

    void add(const ProdN& prod) {
        for (auto& p : prods) {
            // X + XY = X
            // X + X = X
            if (prod.contains(p)) return;
        }
        // XY + X = X
        prods.erase(std::remove_if(prods.begin(), prods.end(), [&] (const ProdN& p) {
            return p.contains(prod);
        }), prods.end());

        prods.push_back(prod);
    }

    void add(const Sum& other) {
        for (auto& prod : other.prods)
            add(prod);
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
std::vector<MinTerm<Nbits>> simplify_chart(const PrimeChart<Nbits>& chart) {
    // Petrick's method
    Sum<Nbits> sum;
    for (auto& pair : chart.columns) {
        if (sum.size() == 0) {
            for (auto& term : pair.second) {
                Prod<Nbits> prod;
                prod.mul(term);
                sum.add(prod);
            }
        } else {
            // Distribute
            Sum<Nbits> new_sum;
            for (auto& term : pair.second) {
                for (auto& prod : sum.prods) {
                    Prod<Nbits> new_prod(prod);
                    new_prod.mul(term);
                    new_sum.add(new_prod);
                }
            }
            sum = std::move(new_sum);
        }
    }
    // Find best product
    size_t min_terms = std::numeric_limits<size_t>::max();
    size_t min_literals = std::numeric_limits<size_t>::max();
    const Prod<Nbits>* found = nullptr;
    for (auto& prod : sum.prods) {
        size_t terms = prod.size();
        if (terms < min_terms) {
            size_t literals = prod.count_literals();
            if (literals < min_literals) {
                found = &prod;
                min_terms = terms;
                min_literals = literals;
            }
        }
    }
    assert(found);
    return found->terms;
}

template <size_t Nbits>
std::vector<MinTerm<Nbits>> minimize_boolean(
    const std::vector<typename MinTerm<Nbits>::IntTypeN>& on_values,
    const std::vector<typename MinTerm<Nbits>::IntTypeN>& dc_values)
{
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
    chart.remove_essentials(solution);
    if (chart.size() == 0)
        return solution;

    auto simplified = simplify_chart(chart);
    solution.insert(solution.end(), simplified.begin(), simplified.end());
    return solution;
}

int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();
    //std::vector<uint8_t> on { 9, 12, 13, 15 };
    //std::vector<uint8_t> dc { 1, 4, 5, 7, 8, 11, 14 };
    std::vector<uint8_t> on { 4, 8, 10, 11, 12, 15 };
    std::vector<uint8_t> dc { 9, 14 };
    auto solution = minimize_boolean<4>(on, dc);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " us" << std::endl;
    for (auto& term : solution)
        std::cout << term << std::endl;
    return 0;
}
