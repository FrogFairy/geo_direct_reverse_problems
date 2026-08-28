#include "primitives.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
#include <vector>
#include <limits>
#include <windows.h>

double to_radians(Degrees degrees)
{
    double deg = std::abs(degrees.degrees) + degrees.minutes / 60.0 + degrees.seconds / 3600.0;
    deg = degrees.degrees < 0 ? -deg : deg;
    return deg / 180 * M_PI;
}

std::vector<double> to_seconds(const std::vector<double>& data)
{
    std::vector<double> result;
    for (double val : data)
    {
        result.push_back(to_seconds(val));
    }
    return result;
}

Degrees to_degrees(double radians)
{
    double degrees;
    double minutes;
    double seconds;

    minutes = std::modf(std::abs(radians) / M_PI * 180, &degrees) * 60;
    seconds = std::modf(minutes, &minutes) * 60;
    degrees = radians < 0 ? -degrees : degrees;
    return Degrees{(int) degrees, (int) minutes, seconds};    
}

std::vector<Degrees> to_degrees(std::vector<double>& data)
{
    std::vector<Degrees> result;
    for (double val : data)
    {
        result.push_back(to_degrees(val));
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, const Degrees& data)
{
    os << data.degrees << "\370"
        << data.minutes << "' "
        << data.seconds << "''";
    return os;
}

std::ostream& operator<<(std::ostream& os, const PointPrinter& pp)
{
    if (pp.show_radians)
    {
        os << "B" << pp.index << " = " << to_degrees(pp.point->B) << " = " << pp.point->B << " рад\n";
        os << "L" << pp.index << " = " << to_degrees(pp.point->L) << " = " << pp.point->L << " рад\n";
        os << "A" << pp.index << " = " << to_degrees(pp.point->A) << " = " << pp.point->A << " рад";
    } else
    {
        os << "B" << pp.index << " = " << to_degrees(pp.point->B) << "\n";
        os << "L" << pp.index << " = " << to_degrees(pp.point->L) << "\n";
        os << "A" << pp.index << " = " << to_degrees(pp.point->A);
    }
    return os;
}

std::istream& operator>>(std::istream& is, Degrees& data)
{
    is >> data.degrees >> data.minutes >> data.seconds;
    return is;
}

std::string alignment(const std::string& data, int width, char align, char filling)
{
    int left = 0;
    int right = 0;
    switch (align)
    {
        case '^':
            left = (width - data.size()) / 2;
            right = width - data.size() - left;
            break;
        case '<':
            right = width - data.size();
            break;
        case '>':
            left = width - data.size();
            break;
    }
    left = left >= 0 ? left : 0;
    right = right >= 0 ? right : 0;

    return std::string(left, filling) + data + std::string(right, filling);
}

std::string to_string_with_precision(double value, int precision)
{
    std::ostringstream out;
    out.precision(precision);
    out << std::fixed << value;
    return out.str();
}

std::vector<std::string> data_to_string(const std::vector<double>& data)
{
    std::vector<std::string> result;
    for (double val : data)
    {
        result.push_back(to_string_with_precision(val, std::numeric_limits<double>::max_digits10));
    }
    return result;
}

std::vector<std::string> data_to_string(const std::vector<Degrees>& data)
{
    std::vector<std::string> result;
    for (Degrees val : data)
    {
        std::ostringstream out;
        out << val;
        result.push_back(out.str());
    }
    return result;
}

void Problem::input_ellipsoid()
{
    double b, e2, e_2;
    out << "Параметры эллипсоида." << std::endl;
    out << "Введите малую полуось эллипсоида (м): b = ";
    std::cin >> b;
    out << StreamType::FileOnly << b << std::endl << StreamType::Both;
    out << "Введите квадрат эксцентриситета: e^2 = ";
    std::cin >> e2;
    out << StreamType::FileOnly << e2 << std::endl << StreamType::Both;
    out << "Введите квадрат второго эксцентриситета: (e')^2 = ";
    std::cin >> e_2;
    out << StreamType::FileOnly << e_2 << std::endl << StreamType::Both;

    this->ellipsoid = Ellipsoid{b, e2, e_2};
}

void Problem::set_output()
{
    std::cout << "Хотите записать результаты в файл? y - да, n - нет: ";
    std::string ans;
    std::cin >> ans;

    while (ans != "y" && ans != "n")
    {
        std::cout << "Хотите записать результаты в файл? y - да, n - нет: ";
        std::cin >> ans;
    }

    if (ans == "y")
    {
        std::string filename;
        std::cout << "Введите название файла: ";
        std::cin >> filename;
        out.enableFileLogging(filename);
    }
}

void Problem::init()
{
    set_output();
    input_ellipsoid();
}

Coefficients ProblemMath::calc_coefficients_simplified(const Ellipsoid& ellipsoid, double cos2_A0)
{
    double A = ellipsoid.b * (1 + 1.0 / 32 * 5 / 256 * std::pow(ellipsoid.e_2, 3)) +
                ellipsoid.b * ellipsoid.e_2 / 4 *
                (1 - 9.0 / 4 * 5 / 256 * std::pow(ellipsoid.e_2, 2)) * cos2_A0 +
                ellipsoid.b * std::pow(ellipsoid.e_2, 2) * 3 / 64 *
                (-1 + 5.0 / 8 * ellipsoid.e_2) * std::pow(cos2_A0, 2);
    double B = ellipsoid.b * ellipsoid.e_2 / 8 *
                (1 - 9.0 / 2 * 15 / 1024 * std::pow(ellipsoid.e_2, 2)) * cos2_A0 +
                ellipsoid.b * std::pow(ellipsoid.e_2, 2) / 32 *
                (-1 + 3.0 / 64 * 15 * ellipsoid.e_2) * std::pow(cos2_A0, 2);
    double C = ellipsoid.b * std::pow(ellipsoid.e_2, 3) / 16 / 512 * 9 * 3 +
                ellipsoid.b * std::pow(ellipsoid.e_2, 2) / 128 *
                (1 - 3.0 / 8 * 3 * ellipsoid.e_2) * std::pow(cos2_A0, 2);
    return Coefficients{A, B, C};
}

Coefficients ProblemMath::calc_coefficients_standart(const Ellipsoid& ellipsoid, double cos2_A0)
{
    double A = ellipsoid.b * (1 + 1.0 / 32 * 5 / 256 * std::pow(ellipsoid.e_2, 3)) +
                ellipsoid.b * ellipsoid.e_2 / 4 *
                (1 - 9.0 / 4 * 5 / 256 * std::pow(ellipsoid.e_2, 2)) * cos2_A0 +
                ellipsoid.b * std::pow(ellipsoid.e_2, 2) * 3 / 64 *
                (-1 + 5.0 / 8 * ellipsoid.e_2) * std::pow(cos2_A0, 2);
    double B = ellipsoid.b * 1.0 / 32 * 15 / 1024 * std::pow(ellipsoid.e_2, 3) +
                ellipsoid.b * ellipsoid.e_2 / 8 *
                (1 - 9.0 / 2 * 15 / 1024 * std::pow(ellipsoid.e_2, 2)) * cos2_A0 +
                ellipsoid.b * std::pow(ellipsoid.e_2, 2) / 32 *
                (-1 + 3.0 / 64 * 15 * ellipsoid.e_2) * std::pow(cos2_A0, 2);
    double C = -ellipsoid.b * 1.0 / 32 * 3 / 512 * std::pow(ellipsoid.e_2, 3) +
                ellipsoid.b * std::pow(ellipsoid.e_2, 3) / 16 / 512 * 9 * 3 * cos2_A0 +
                ellipsoid.b * std::pow(ellipsoid.e_2, 2) / 128 *
                (1 - 3.0 / 8 * 3 * ellipsoid.e_2) * std::pow(cos2_A0, 2);
    return Coefficients{A, B, C};
}

Coefficients ProblemMath::calc_coefficients_accurate(const Ellipsoid& ellipsoid, double cos2_A0)
{
    double A = ellipsoid.b * (1 + 1.0 / 32 * 5 / 256 * std::pow(ellipsoid.e_2, 3)) +
                ellipsoid.b * ellipsoid.e_2 / 4 *
                (1 - 9.0 / 4 * 5 / 256 * std::pow(ellipsoid.e_2, 2)) * cos2_A0 +
                ellipsoid.b * std::pow(ellipsoid.e_2, 2) * 3 / 64 *
                (-1 + 5.0 / 8 * ellipsoid.e_2) * std::pow(cos2_A0, 2);
    double B = ellipsoid.b * 1.0 / 32 * 15 / 1024 * std::pow(ellipsoid.e_2, 3) +
                ellipsoid.b * ellipsoid.e_2 / 8 *
                (1 - 9.0 / 2 * 15 / 1024 * std::pow(ellipsoid.e_2, 2)) * cos2_A0 +
                ellipsoid.b * std::pow(ellipsoid.e_2, 2) / 32 *
                (-1 + 3.0 / 64 * 15 * ellipsoid.e_2) * std::pow(cos2_A0, 2);
    double C = -ellipsoid.b * 1.0 / 32 * 3 / 512 * std::pow(ellipsoid.e_2, 3) +
                ellipsoid.b * std::pow(ellipsoid.e_2, 3) / 16 / 512 * 9 * 3 * cos2_A0 +
                ellipsoid.b * std::pow(ellipsoid.e_2, 2) / 128 *
                (1 - 3.0 / 8 * 3 * ellipsoid.e_2) * std::pow(cos2_A0, 2);
    return Coefficients{A, B, C};
}

Coefficients ProblemMath::calc_coefficients(const Way& way, const Ellipsoid& ellipsoid, double cos2_A0)
{
    switch (way)
    {
        case Simplified:
            return ProblemMath::calc_coefficients_simplified(ellipsoid, cos2_A0);
        case Standart:
            return ProblemMath::calc_coefficients_standart(ellipsoid, cos2_A0);
        case Accurate:
            return ProblemMath::calc_coefficients_accurate(ellipsoid, cos2_A0);
    }
    return Coefficients{0, 0, 0};
}

// КОСТЫЛЬ ДЛЯ ПЕРЕВОДА В 1251
std::string StreamSplitter::Cp1125ToCp1251(const std::string& cp1125Text) {
    int wideSize = MultiByteToWideChar(866, 0, cp1125Text.c_str(), -1, nullptr, 0);
    if (wideSize <= 0) return "";
    
    std::wstring wide(wideSize, L'\0');
    MultiByteToWideChar(866, 0, cp1125Text.c_str(), -1, wide.data(), wideSize);
    wide.pop_back();
    
    int cp1251Size = WideCharToMultiByte(1251, 0, wide.c_str(), -1, nullptr, 0, nullptr, nullptr);
    if (cp1251Size <= 0) return "";
    
    std::string result(cp1251Size, '\0');
    WideCharToMultiByte(1251, 0, wide.c_str(), -1, result.data(), cp1251Size, nullptr, nullptr);
    result.pop_back();
    
    return result;
}

// КОСТЫЛЬ ДЛЯ ПЕРЕВОДА В 1251
void StreamSplitter::SaveAsWindows1251(const std::stringstream& cp1125Text) {
    std::string s = cp1125Text.str();
    std::string cp1251Text = Cp1125ToCp1251(s);
    file->write(cp1251Text.c_str(), cp1251Text.size());
}

void StreamSplitter::enableFileLogging(const std::string& filename) {
    file = std::make_unique<std::ofstream>(filename, std::ios::app | std::ios::binary);
    if (file->is_open()) {
        writeToFile = true;
        this->filename = filename;
    } else {
        std::cout << "Ошибка при открытии файла: " << filename << std::endl;
    }
}

void StreamSplitter::disableFileLogging() {
    if (file && file->is_open()) {
        file->flush();
        file->close();
    }
    file.reset();
    writeToFile = false;
}

StreamSplitter& StreamSplitter::operator<<(StreamType stream_type) {
    writeToConsole = stream_type == StreamType::FileOnly ? false : true;
    writeToFile = stream_type == StreamType::ConsoleOnly ? false : true;
    return *this;
}

StreamSplitter& StreamSplitter::operator<<(std::ostream& (*manip)(std::ostream&)) {
    if (writeToConsole) {
        console << manip;
    }
    if (writeToFile && file && file->is_open()) {
        std::stringstream ss;
        ss << manip;
        SaveAsWindows1251(ss);
    }
    return *this;
}