# Root Separation and Refinement Using Newton's Method

This repository contains a program for solving equations with one variable using **Newton's Method**. The program identifies roots of a given function on a specified interval, refines them using Newton's Method, and provides performance metrics including convergence behavior, number of iterations, and computation time.

---

## Problem Statement

Solving equations with one variable involves finding the values of the variable that satisfy the equation $( f(x) = 0 \)$. This program implements **Newton's Method** to refine roots identified within user-defined intervals.

### Input:
- A function $f(x) = 0.5x^2 - 10 + 2^{-x}$, its first and second derivatives.
- Interval boundaries $[a, b]$.
- Step size $( h \)$ for root separation.
- Tolerances:
  - $( \epsilon_1 \)$: precision for arguments.
  - $( \epsilon_2 \)$: precision for function values.

### Output:
- Root $( x \)$.
- Function value at the root $f(x)$.
- Number of iterations $( n \)$.
- Number of function evaluations.
- Convergence parameter: $\alpha = \frac{|x_{n+1} - x_n|}{|x_n - x_{n-1}|^2}$
  
- Computation time for each root.

---

## Implementation Details

### Algorithm Description: Newton's Method

Newton's Method refines a root approximation iteratively using the formula:

$$
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
$$

**Steps:**

1. **Root Separation:**
   - Divide the interval $[a, b]$ into subintervals of size $( h \)$.
   - Identify subintervals where $f(x)$ changes sign, indicating a root.

2. **Initial Guess Selection:**
   - Choose an initial guess $( x_0 \)$ based on the sign of $f''(x)$ at the subinterval boundaries.

3. **Newton's Iterative Process:**
   - Apply the Newton iteration until:
     $|x_{n+1} - x_n| < \epsilon_1$, and
     $|f(x_n)| < \epsilon_2$.

4. **Convergence Analysis:**
   - Calculate the convergence parameter to evaluate the quadratic convergence property.

5. **Performance Metrics:**
   - Count function evaluations for $f(x)$, $f'(x)$, and $f''(x)$.
   - Measure the time taken for each root refinement.

---

## Example

Given the function:

$$
f(x) = 0.5x^2 - 10 + 2^{-x}
$$

For the interval $[-5, 5]$ with tolerances $\epsilon_1 = 10^{-6}$ and $\epsilon_2 = 10^{-6}$, and step size $( h = 7)$, the program computes the roots and outputs the following:

### Sample Output:
```
Enter the interval boundaries [a, b]: -5 5
Enter the required accuracy for root determination by argument (eps1): 1e-6
Enter the required accuracy for root determination by function (eps2): 1e-6

Root: -2.68003
Accuracy: 2.31997e+00
Function value f(Xi):  1.82712e-08
Number of iterations: 5
Total number of function evaluations and their derivatives: 15
Execution time: 20 microseconds
Convergence parameter: 0.28818

Root: 4.46198
Accuracy: 5.38022e-01
Function value f(Xi):  5.12005e-09
Number of iterations: 3
Total number of function evaluations and their derivatives: 25
Execution time: 2 microseconds
Convergence parameter: 0.11528

Press Enter to exit...
```

---

## Code Structure

### Main Function
Handles user input and initiates root separation and refinement:
```cpp
int main() {
    setlocale(LC_ALL, "");
    double a, b, eps1, eps2, h = 7.0;
    cout << "Enter the interval boundaries [a, b]: ";
    cin >> a >> b;
    cout << "Enter the required accuracy for root determination by argument (eps1): ";
    cin >> eps1;
    cout << "Enter the required accuracy for root determination by function (eps2): ";
    cin >> eps2;
    double root, convergence;
    double f_root;
    int n;
    auto start = chrono::high_resolution_clock::now();
    double prev_f = f(a);
    for (double x = a; x <= b; x += h) {
        double x2 = min(x + h, b);
        double curr_f = f(x2);
        if (prev_f * curr_f <= 0) {
            double x0 = (prev_f * f_double_prime(x) > 0) ? x : x2;
            root = newton_method(x0, eps1, eps2, n, convergence);
            f_root = f(root);
            auto stop = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
            cout << "\nRoot: " << std::fixed << std::setprecision(5) << root << endl;
            cout << "Function value f(Xi): " << f_root << endl;
            cout << "Number of iterations: " << n << endl;
            cout << "Execution time: " << duration.count() << " microseconds" << endl;
            cout << "Convergence parameter: " << std::fixed << std::setprecision(5) << convergence << endl;
            start = chrono::high_resolution_clock::now();
        }
        prev_f = curr_f;
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout << "Press Enter to exit...";
    std::cin.get();
    return 0;
}
```

### Newton's Method Implementation
Refines root approximations:
```cpp
double newton_method(double x, double eps1, double eps2, int& n, double& convergence) {
    double x_prev = x;
    double x_next = x - f(x) / f_prime(x);
    double x_prev_prev = x;
    n = 1;
    double fx_next = f(x_next);
    while (abs(x_next - x_prev) >= eps1 && abs(fx_next) >= eps2) {
        x_prev_prev = x_prev;
        x_prev = x_next;
        x_next = x_prev - f(x_prev) / f_prime(x_prev);
        fx_next = f(x_next);
        n++;
    }
    convergence = abs((x_next - x_prev) / pow(abs(x_prev - x_prev_prev), 2));
    return x_next;
}
```

### Function Definitions
Defines \( f(x) \), \( f'(x) \), and \( f''(x) \):
```cpp
double f(double x) {
    return 0.5 * pow(x, 2) - 10 + pow(2, -x);
}

double f_prime(double x) {
    return x - pow(2, -x) * log(2);
}

double f_double_prime(double x) {
    return 1 + pow(2, -x) * log(2) * log(2);
}
```

---

## Requirements
- **Language**: C++
- **Libraries**: `<iostream>`, `<cmath>`, `<chrono>`, `<iomanip>`
- **Compiler**: Any standard-compliant C++ compiler

---

## Usage
Clone the repository and execute the program:
```bash
git clone https://github.com/Nissmoline/Root-Separation-and-Refinement-Program-using-Newton-s-Method.git
cd Root-Separation-and-Refinement-Program-using-Newton-s-Method
g++ -o root_refinement main.cpp
./root_refinement
```

---

## License
This project is open-source and available under the MIT License.
