# MinBool

This project is a boolean expression simplifier in C++. It uses the Quine–McCluskey algorithm to compute the set of prime implicants, and then runs Petrick's method to find the minimum sum of products solution.

This project does not depend on any other library, except of course the standard C++ library.

## Building and Running

    g++ min_bool.cpp -pedantic -Wall -Wextra -O3 -DNDEBUG -march=native -std=c++11 -o min_bool
    ./min_bool

## Using the Simplifier

The boolean function to simplify is represented as two vectors containing the set of values for which the function evaluates to `true`, and the set of values for which the values of the function does not matter (it can be either 0 or 1).

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
auto solution = minimize_boolean<4>(on, dc);
```
