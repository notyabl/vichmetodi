2 задание код:

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <sstream>
#include <iomanip>

using namespace std;

// Решение СЛАУ методом Гаусса (для малых систем) с выбором гл. элемента по столбцу
// Решение СЛАУ методом Гаусса (для малых систем) – возвращает false, если система вырождена, и ture тогда система решена успешно.
bool solveLinearSystem(vector<vector<double>> A, vector<double> b, vector<double>& x) {
    int n = A.size();
    // Прямой ход
    for (int i = 0; i < n; ++i) {
        // Поиск главного элемента
        int maxRow = i;
        for (int k = i + 1; k < n; ++k)
            if (fabs(A[k][i]) > fabs(A[maxRow][i]))
                maxRow = k;
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        double pivot = A[i][i];
        if (fabs(pivot) < 1e-12) return false;   // вырожденность

        for (int j = i; j < n; ++j) A[i][j] /= pivot;
        b[i] /= pivot;

        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i];
            for (int j = i; j < n; ++j)
                A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }
    // Обратный ход: начиная с последнего уравнения, вычисляются неизвестные.
    x.assign(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j)
            x[i] -= A[i][j] * x[j];
    }
    return true;
}
// построить аппроксимирующий полином степени deg по точкам (xw[i], yw[i]) 
// методом наименьших квадратов и вычислить его значение в точке x.
double localFit(double x, const vector<double>& xw, const vector<double>& yw, int deg) {
    int m = xw.size();
    if (m <= deg) deg = m - 1;
    if (deg < 0) return yw[0];

    int ncoeff = deg + 1;
    vector<vector<double>> A(ncoeff, vector<double>(ncoeff, 0.0));
    vector<double> B(ncoeff, 0.0);

    for (int i = 0; i < m; ++i) {
        double xi = xw[i], yi = yw[i];
        double power = 1.0;
        for (int p = 0; p < ncoeff; ++p) {
            B[p] += yi * power;
            for (int q = 0; q < ncoeff; ++q) {
                A[p][q] += power * pow(xi, q);
            }
            power *= xi;
        }
    }

    vector<double> coeff;
    if (!solveLinearSystem(A, B, coeff)) {
        // Если система вырождена, вернуть среднее арифметическое yw
        double sum = 0.0;
        for (double v : yw) sum += v;
        return sum / m;
    }

    double res = 0.0, xpow = 1.0;
    for (int p = 0; p < ncoeff; ++p) {
        res += coeff[p] * xpow;
        xpow *= x;
    }
    return res;
}

// основная функция сглаживания –
// для каждой точки исходных данных строит скользящее окно и вычисляет сглаженное значение.
vector<double> smooth(const vector<double>& x, const vector<double>& y, int radius, int deg) {
    int n = x.size();
    vector<double> ys(n);
    for (int i = 0; i < n; ++i) {
        int left = max(0, i - radius);
        int right = min(n - 1, i + radius);
        vector<double> xw, yw;
        for (int j = left; j <= right; ++j) {
            xw.push_back(x[j]);
            yw.push_back(y[j]);
        }
        ys[i] = localFit(x[i], xw, yw, deg);
    }
    return ys;
}
// Генерация тестовых данных: y = sin(x) + шум
void genTestData(vector<double>& x, vector<double>& y, int n, double xmin, double xmax) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> noise(-0.2, 0.2);
    x.resize(n);
    y.resize(n);
    double step = (xmax - xmin) / (n - 1);
    for (int i = 0; i < n; ++i) {
        x[i] = xmin + i * step;
        y[i] = sin(x[i]) + noise(gen);
    }
}

// Ввод данных с клавиатуры
void inputData(vector<double>& x, vector<double>& y) {
    int n;
    cout << "Количество точек: ";
    cin >> n;
    x.resize(n);
    y.resize(n);
    cout << "Введите x y (через пробел, каждая точка с новой строки):\n";
    for (int i = 0; i < n; ++i)
        cin >> x[i] >> y[i];
}

// Форматирование числа: точка -> запятая
string fmt(double val) {
    ostringstream oss;
    oss << fixed << setprecision(6) << val;
    string s = oss.str();
    size_t pos = s.find('.');
    if (pos != string::npos) s[pos] = ',';
    return s;
}
int main() {
    setlocale(LC_ALL, "Russian");

    cout << "=== Локальное сглаживание ===\n";
    cout << "1 - тестовые данные (sin(x) + шум)\n";
    cout << "2 - ввести свои данные\n";
    int choice;
    cin >> choice;

    vector<double> x, y;
    if (choice == 1) {
        int n; double xmin, xmax;
        cout << "n = "; cin >> n;
        cout << "xmin xmax = "; cin >> xmin >> xmax;
        genTestData(x, y, n, xmin, xmax);
    }
    else {
        inputData(x, y);
    }

    int rad, deg;
    cout << "Радиус окна (соседей слева и справа): "; cin >> rad;
    cout << "Степень полинома: "; cin >> deg;

    // Проверка радиуса
    if (rad < 1) {
        cout << "Ошибка: радиус не может быть меньше 1. Установлено значение 1.\n";
        rad = 1;
    }
    cout << "Параметры сглаживания: радиус = " << rad << ", степень = " << deg << endl;

    vector<double> ys = smooth(x, y, rad, deg);

    // Вывод таблицы с тремя столбцами (x, y, ys) через табуляцию
    cout << "\n=== Таблица для копирования в Excel (x, исходное y, сглаженное y) ===\n";
    cout << "x\tисходное y\tсглаженное y\n";
    for (size_t i = 0; i < x.size(); ++i) {
        cout << fmt(x[i]) << "\t" << fmt(y[i]) << "\t" << fmt(ys[i]) << "\n";
    }

    // Дополнительно: вывод вертикальных списков (на случай, если нужно по отдельности)
    cout << "\n=== Вертикальные списки (для копирования каждого столбца отдельно) ===\n";
    cout << "x:\n";
    for (size_t i = 0; i < x.size(); ++i) cout << fmt(x[i]) << "\n";
    cout << "\nисходное y:\n";
    for (size_t i = 0; i < y.size(); ++i) cout << fmt(y[i]) << "\n";
    cout << "\nсглаженное y:\n";
    for (size_t i = 0; i < ys.size(); ++i) cout << fmt(ys[i]) << "\n";

    // Ошибка
    double rmse = 0;
    for (size_t i = 0; i < y.size(); ++i) {
        double d = ys[i] - y[i];
        rmse += d * d;
    }
    rmse = sqrt(rmse / y.size());
    cout << "\nСреднеквадратичная ошибка: " << fmt(rmse) << "\n";

    return 0;
}
