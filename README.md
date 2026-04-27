// 3 задание.
#include <iostream>
#include <cmath>

using namespace std;

// Изменяем уравнение. 
double f(double x) {
    return pow(x, 3) - 0.2 * pow(x, 2) - 0.2 * x - 1.2;
}

// Метод Вегстейна
void wegstein(double x0, double x1, double eps, int maxIter) {
    double x2, f0, f1, f2;
    int iter = 0;

    f0 = f(x0);
    f1 = f(x1);

    cout << "Начальное x0 = " << x0 << ", f(x0) = " << f0 << endl;
    cout << "Начальное x1 = " << x1 << ", f(x1) = " << f1 << endl;

    while (iter < maxIter) {
        if (fabs(f1 - f0) < 1e-15) {
            cout << "Знаменатель близок к нулю. Останов.\n";
            return;
        }
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        f2 = f(x2);

        cout << "Итерация " << iter + 1 << ": x = " << x2 << ", f(x) = " << f2 << endl;

        if (fabs(f2) < eps) {
            cout << "Корень найден: x = " << x2 << ", f(x) = " << f2 << endl;
            return;
        }

        x0 = x1; f0 = f1;
        x1 = x2; f1 = f2;
        iter++;
    }
    cout << "Метод не сошёлся за " << maxIter << " итераций.\n";
}

int main() {
    setlocale(LC_ALL, "Russian");

    double a, b, eps;
    int maxIter;

    cout << "Решение уравнения f(x)=0 методом Вегстейна\n";
    cout << "Уравнение задано в функции f(x) (смотрите код)\n";

    cout << "Введите отрезок [a, b]:\n";
    cout << "a = "; cin >> a;
    cout << "b = "; cin >> b;

    if (f(a) * f(b) >= 0) {
        cout << "На отрезке нет корня или корней несколько.\n";
        return 1;
    }

    cout << "Введите точность (eps): "; cin >> eps;
    cout << "Введите максимальное число итераций: "; cin >> maxIter;

    wegstein(a, b, eps, maxIter);

    return 0;
}

//примеры:
1 пример return pow(x,4) - 3*x*x + 1; 
a = 0
b = 1
eps = 0.000001
maxIter = 30
2 пример  return x * x * x - 2 * x - 5;
a = 2
b = 3
eps = 0.000001
maxIter = 30
3 пример return x*x - 4;
[1, 1.5]
