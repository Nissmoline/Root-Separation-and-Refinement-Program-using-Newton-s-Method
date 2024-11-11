#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>

using namespace std;

int calc_count = 0; 

double f(double x) {
    calc_count++;
    return 0.5 * pow(x, 2) - 10 + pow(2, -x);
}

double f_prime(double x) {
    calc_count++;
    return x - pow(2, -x) * log(2);
}

double f_double_prime(double x) {
    calc_count++;
    return 1 + pow(2, -x) * log(2) * log(2);
}

double newton_method(double x, double eps1, double eps2, int& n, double& convergence) {
    double x_prev = x;
    double x_next = x - f(x) / f_prime(x);
    double x_prev_prev = x;
    n = 1;
    convergence = 0;
    double fx_next = f(x_next);
    while (abs(x_next - x_prev) >= eps1 && abs(fx_next) >= eps2) {
        x_prev_prev = x_prev;
        x_prev = x_next;
        double fx_prev = fx_next;
        double f_prime_x_prev = f_prime(x_prev);
        x_next = x_prev - fx_prev / f_prime_x_prev;
        fx_next = f(x_next);
        n++;
    }

    convergence = abs((x_next - x_prev) / pow(abs(x_prev - x_prev_prev), 2));
    return x_next;
}

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
            double x0;
            double curr_f_double_prime = f_double_prime(x);
            if (prev_f * curr_f_double_prime > 0) {
                x0 = x;
            }
            else {
                x0 = x2;
            }
            root = newton_method(x0, eps1, eps2, n, convergence);
            f_root = f(root);
            auto stop = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
            cout << "\nRoot: " << std::fixed << std::setprecision(5) << root << endl;
            cout << std::scientific << "Accuracy: " << abs(root - x0) << endl;
            cout << "Function value f(Xi): " << f_root << endl;
            cout << "Number of iterations: " << n << endl;
            cout << "Total number of function evaluations and their derivatives: " << calc_count << endl;
            cout << "Execution time: " << duration.count() << " microseconds" << endl;
            cout << "Convergence parameter: " << std::fixed << std::setprecision(5) << convergence << endl << endl;
            start = chrono::high_resolution_clock::now();
        }
        prev_f = curr_f;
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout << "Press Enter to exit...";
    std::cin.get();
    return 0;
}