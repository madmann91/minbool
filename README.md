# MinBool

This project is a boolean expression simplifier in C++. It uses the Quine–McCluskey algorithm to compute the set of prime implicants, and then iteratively extracts prime essentials and simplifies the implicant chart. When no chart simplification can be applied and no prime essential can be removed, it relies on heuristics to find a good prime to extract from the chart.

This project does not depend on any other library, except of course the standard C++ library. It is distributed under the MIT license.

## Building and Running

    g++ minbool.cpp -pedantic -Wall -Wextra -O3 -DNDEBUG -march=native -std=c++11 -o minbool
    ./minbool

## Using the Simplifier

The boolean function to simplify is represented as two vectors containing the set of values for which the function evaluates to 1, and the set of values for which the values of the function does not matter (it can be either 0 or 1).

Consider the following example for a function of 4 variables (X mark "don't care" values):

|    | A | B | C | D | f(A, B, C, D) |
|----|---|---|---|---|---------------|
| 0  | 0 | 0 | 0 | 0 |       0       |
| 1  | 0 | 0 | 0 | 1 |       0       |
| 2  | 0 | 0 | 1 | 0 |       0       |
| 3  | 0 | 0 | 1 | 1 |       0       |
| 4  | 0 | 1 | 0 | 0 |       1       |
| 5  | 0 | 1 | 0 | 1 |       0       |
| 6  | 0 | 1 | 1 | 0 |       0       |
| 7  | 0 | 1 | 1 | 1 |       0       |
| 8  | 1 | 0 | 0 | 0 |       1       |
| 9  | 1 | 0 | 0 | 1 |       X       |
| 10 | 1 | 0 | 1 | 0 |       1       |
| 11 | 1 | 0 | 1 | 1 |       1       |
| 12 | 1 | 1 | 0 | 0 |       1       |
| 13 | 1 | 1 | 0 | 1 |       0       |
| 14 | 1 | 1 | 1 | 0 |       X       |
| 15 | 1 | 1 | 1 | 1 |       1       |

The code to run the simplifier on this function is:

```cpp
std::vector<uint8_t> on { 4, 8, 10, 11, 12, 15 };
std::vector<uint8_t> dc { 9, 14 };
std::vector<MinTerm<4>> solution = minimize_boolean<4>(on, dc);
```

The function `minimize_boolean` takes the number of bits to consider as a template argument, and the two vectors containing the definition of the function as described above. The result is a vector of `MinTerm` objects representing the sum of products found during minimization. For the example above, the result will be (one `MinTerm` per line):

    -100
    1-1-
    10--

This is equivalent to the minimal boolean expression:

    BC'D' + AC + AB'

As demonstrated in this example, a `MinTerm` indicates, for each variable, if it is present in the product (`1`), present but negated (`0`), or not present at all (`-`). The state of every variable in the product is read from left to right (i.e. in the example, the leftmost symbol represents the state of variable A, and the right most symbol the state of variable D).
