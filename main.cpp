#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>
#include <numeric>

#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;

class FFT {
public:
    using Complex = complex<double>;

    static bool transform(vector<Complex>& data, bool inverse) {
        // Проверяет длину данных и запускает рекурсивное преобразование
        size_t n = data.size();
        if (!is_valid_length(n)) {
            cerr << "Ошибка: Длина данных " << n << " не является произведением степеней 2, 3 или 5." << endl;
            return false;
        }
        fft_recursive(data, inverse);
        return true;
    }

    // Проверяет, является ли длина данных n произведением степеней 2, 3 и 5.
    static bool is_valid_length(size_t n) {
        if (n == 0) return false;
        // Удаляем все множители 2, 3 и 5 из n
        for (int factor : {2, 3, 5}) {
            while (n > 0 && n % factor == 0) {
                n /= factor;
            }
        }
        // Если после удаления всех множителей 2, 3, 5 осталось 1, то длина допустима
        return n == 1;
    }

private:
    static void fft_recursive(vector<Complex>& a, bool inverse) {
        // Базовый случай рекурсии: массив из 0 или 1 элемента уже преобразован
        size_t n = a.size();
        if (n <= 1) return;
        // Находим наименьший простой множитель (2, 3 или 5)
        size_t d = find_divisor(n);
        // Количество подпоследовательностей (ДПФ меньшего размера)
        size_t m = n / d;

        // 1. Распределение данных по d группам (столбцам матрицы d x m)
        vector<vector<Complex>> groups(d, vector<Complex>(m));
        for (size_t i = 0; i < d; ++i) {
            for (size_t j = 0; j < m; ++j) {
                groups[i][j] = a[j * d + i];
            }
        }

        // 2. Рекурсивные вызовы БПФ для каждой из d групп (для каждого столбца)
        for (auto& group_vector : groups) {
            fft_recursive(group_vector, inverse);
        }

        // 3. Объединение результатов: применение поворотных множителей и d-точечных ДПФ
        const double PI = std::acos(-1.0);
        double base_angle_N = (inverse ? 2.0 : -2.0) * PI / static_cast<double>(n);

        for (size_t k_m = 0; k_m < m; ++k_m) {
            for (size_t k_d = 0; k_d < d; ++k_d) {
                Complex current_sum(0.0, 0.0);
                for (size_t r_idx = 0; r_idx < d; ++r_idx) {

                    // Вычисляем угол для поворотного множителя
                    double phi_twiddle_N = base_angle_N * static_cast<double>(k_m * r_idx);
                    Complex twiddle_factor_N(cos(phi_twiddle_N), sin(phi_twiddle_N));

                    // Вычисляем угол для ядра ДПФ
                    double phi_dft_d = base_angle_N * static_cast<double>(m * k_d * r_idx);
                    Complex dft_kernel_d(cos(phi_dft_d), sin(phi_dft_d));

                    // Суммируем с учетом поворотных множителей
                    current_sum += groups[r_idx][k_m] * twiddle_factor_N * dft_kernel_d;
                }
                // Записываем результат в массив
                a[k_d * m + k_m] = current_sum;
            }
        }

        // 4. Нормировка для обратного преобразования
        if (inverse) {
            for (auto& x : a) {
                x /= static_cast<double>(d);
            }
        }
    }

    static size_t find_divisor(size_t n) {
        // Ищем первый подходящий делитель (2, 3 или 5)
        if (n % 2 == 0) return 2;
        if (n % 3 == 0) return 3;
        if (n % 5 == 0) return 5;
        return n;
    }
};

vector<FFT::Complex> generate_random_data(size_t n) {
    // Создаем вектор случайных комплексных чисел
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(-1.0, 1.0);

    vector<FFT::Complex> data(n);
    for (size_t i = 0; i < n; ++i) {
        // Заполняем каждый элемент случайными значениями
        data[i] = FFT::Complex(dis(gen), dis(gen));
    }
    return data;
}

double calculate_error(const vector<FFT::Complex>& original,
                       const vector<FFT::Complex>& reconstructed) {
    // Проверяем размеры векторов и их пустоту
    if (original.size() != reconstructed.size() || original.empty()) {
        return (original.empty() && reconstructed.empty()) ? 0.0 : 1.0;
    }

    double squared_error_sum = 0.0;
    double original_signal_energy = 0.0;

    for (size_t i = 0; i < original.size(); ++i) {
        // Суммируем квадраты разностей и исходного сигнала
        squared_error_sum += norm(original[i] - reconstructed[i]);
        original_signal_energy += norm(original[i]);
    }

    // Проверяем нулевой случай
    if (original_signal_energy == 0.0) {
        return (squared_error_sum == 0.0) ? 0.0 : 1.0;
    }

    // Возвращаем относительную среднеквадратичную ошибку
    return sqrt(squared_error_sum) / sqrt(original_signal_energy);
}

void print_complex_vector(const vector<FFT::Complex>& vec, const string& label, size_t max_elements = 10) {
    cout << label << endl;
    cout << fixed << setprecision(4);
    size_t elements_to_print = min(vec.size(), max_elements);

    for (size_t i = 0; i < elements_to_print; ++i) {
        double re = vec[i].real();
        double im = vec[i].imag();

        if (std::abs(re) < 1e-9 && std::signbit(re)) re = 0.0;
        if (std::abs(im) < 1e-9 && std::signbit(im)) im = 0.0;

        // Выводим каждый элемент
        cout << "[" << setw(5) << i << "] "
             << setw(8) << re
             << (im < 0 ? " - " : " + ")
             << setw(8) << abs(im) << "j" << endl;
    }


    if (vec.size() > elements_to_print) {
        cout << "\n... [пропущено " << (vec.size() - elements_to_print) << " элементов]\n" << endl;
    } else if (!vec.empty()) {
        cout << endl;
    } else {
        cout << "[пустой вектор]" << endl << endl;
    }
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    size_t n_data;
    // Запрашиваем длину данных у пользователя
    cout << "Введите длину данных (число, являющееся произведением степеней 2, 3, 5, например, 12, 30, 60, 100): ";
    cin >> n_data;

    // Проверяем допустимость длины
    if (!FFT::is_valid_length(n_data)) {
        if (n_data > 0) {
             cerr << "Ошибка: Недопустимая длина данных (" << n_data <<"). Завершение программы." << endl;
        } else if (cin.fail()){
             cerr << "Ошибка: Некорректный ввод. Пожалуйста, введите число." << endl;
        } else {
             cerr << "Ошибка: Длина данных должна быть положительным числом. Завершение программы." << endl;
        }
        return 1;
    }

    // Генерируем и сохраняем исходные данные
    vector<FFT::Complex> data_vector = generate_random_data(n_data);
    vector<FFT::Complex> original_data_vector = data_vector;

    // Выводим исходные данные
    print_complex_vector(original_data_vector, "Исходные данные:");

    // Выполняем прямое БПФ
    if (!FFT::transform(data_vector, false)) {
        cerr << "Ошибка во время прямого преобразования Фурье. Завершение программы." << endl;
        return 1;
    }
    print_complex_vector(data_vector, "Результат прямого преобразования Фурье:");

    // Выполняем обратное БПФ
    if (!FFT::transform(data_vector, true)) {
        cerr << "Ошибка во время обратного преобразования Фурье. Завершение программы." << endl;
        return 1;
    }
    print_complex_vector(data_vector, "Результат обратного преобразования Фурье (восстановленные данные):");

    // Вычисляем и выводим ошибку
    double error = calculate_error(original_data_vector, data_vector);
    cout << scientific << setprecision(8)
         << "Относительная среднеквадратичная ошибка между исходными и восстановленными данными: " << error << endl;

    return 0;
}