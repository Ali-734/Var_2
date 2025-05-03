#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <iomanip>
#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;

class FFT {
public:
    using Complex = complex<double>;

    // Выполняю прямое (inverse=false) или обратное (inverse=true) FFT
    static bool transform(vector<Complex>& data, bool inverse) {
        size_t n = data.size();
        if (!is_valid_length(n)) {
            cerr << "Ошибка: Длина данных " << n << " не кратна 2, 3 или 5." << endl;
            return false;
        }
        fft_recursive(data, inverse);
        if (inverse) {
            // Нормировка обратного преобразования
            for (auto& x : data) x /= static_cast<double>(n);
        }
        return true;
    }


    static bool is_valid_length(size_t n) {
        if (n == 0) return false;
        for (int f : {2, 3, 5}) {
            while (n % f == 0) n /= f;
        }
        return n == 1;
    }

private:
    // Рекурсивная реализация алгоритма Кули-Тьюки
    static void fft_recursive(vector<Complex>& a, bool inverse) {
        size_t n = a.size();
        if (n <= 1) return;

        size_t d = find_divisor(n);
        size_t m = n / d;

        vector<vector<Complex>> groups(d, vector<Complex>(m));
        for (size_t i = 0; i < d; ++i)
            for (size_t j = 0; j < m; ++j)
                groups[i][j] = a[j * d + i];

        for (auto& grp : groups)
            fft_recursive(grp, inverse);

        double angle = (inverse ? 2 : -2) * M_PI / static_cast<double>(n);

        for (size_t k = 0; k < m; ++k) {
            for (size_t j = 0; j < d; ++j) {
                Complex sum(0.0, 0.0);
                for (size_t r = 0; r < d; ++r) {
                    double phi = angle * (j * r * m);
                    sum += groups[r][k] * Complex(cos(phi), sin(phi));
                }
                a[k * d + j] = sum;
            }
        }
    }

    // Возвращаю первый делитель из {2,3,5}
    static size_t find_divisor(size_t n) {
        for (int d : {2, 3, 5}) {
            if (n % d == 0) return d;
        }
        return 1;
    }
};

// Генерация случайных комплексных данных
vector<complex<double>> generate_random_data(size_t n) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(-1.0, 1.0);

    vector<complex<double>> data(n);
    for (size_t i = 0; i < n; ++i)
        data[i] = complex<double>(dis(gen), dis(gen));
    return data;
}

// Вычисление относительной среднеквадратической ошибки
double calculate_error(const vector<complex<double>>& original,
                       const vector<complex<double>>& reconstructed) {
    double squared_error = 0.0;
    double signal_energy = 0.0;
    for (size_t i = 0; i < original.size(); ++i) {
        squared_error += norm(original[i] - reconstructed[i]);
        signal_energy += norm(original[i]);
    }
    if (signal_energy == 0.0) return 0.0;
    return sqrt(squared_error) / sqrt(signal_energy);
}


void print_complex_vector(const vector<complex<double>>& vec,
                         const string& label,
                         size_t max_elements = 100) {  // Значение по умолчанию
    cout << label << endl;
    cout << fixed << setprecision(4);

    const size_t total = vec.size();
    const bool need_truncate = (total > max_elements);
    const size_t to_show = need_truncate ? max_elements : total;

    for (size_t i = 0; i < to_show; ++i) {
        double re = vec[i].real();
        double im = vec[i].imag();
        re = (re == -0.0) ? 0.0 : re;
        im = (im == -0.0) ? 0.0 : im;

        cout << "[" << setw(5) << i << "] "
             << setw(8) << re << " + "
             << setw(8) << im << "j" << endl;
    }

    if (need_truncate) {
        cout << "\n... [пропущено " << (total - max_elements)
             << " элементов]\n" << endl;
    }
    else {
        cout << endl;
    }
}

int main() {
    SetConsoleOutputCP(CP_UTF8);

    size_t n;
    cout << "Введите длину данных (кратную 2, 3, 5): ";
    cin >> n;

    if (!FFT::is_valid_length(n)) {
        cerr << "Ошибка: Недопустимая длина данных." << endl;
        return 1;
    }

    auto data = generate_random_data(n);
    auto original = data;

    print_complex_vector(original, "Исходные данные:");

    if (!FFT::transform(data, false)) return 1;
    print_complex_vector(data, "Прямое преобразование Фурье:");

    if (!FFT::transform(data, true)) return 1;
    print_complex_vector(data, "Обратное преобразование Фурье:");

    double error = calculate_error(original, data);
    cout << scientific << "Среднеквадратичная ошибка: " << error << endl;

    return 0;
}
