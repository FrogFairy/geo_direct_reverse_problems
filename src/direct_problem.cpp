#include "direct_problem.h"

#include <iostream>
#include <cmath>

void DirectProblem::input_ellipsoid()
{
    out << "Прямая задача." << std::endl;
    Problem::input_ellipsoid();
}

void DirectProblem::input_data(DirectState& state) 
{
    Degrees B1, L1, A1;
    double s;
    
    out << "Входные данные." << std::endl;
    out << "Введите широту (\"ГРАДУСЫ МИНУТЫ СЕКУНДЫ\"): B1 = ";
    std::cin >> B1;
    out << StreamType::FileOnly << B1 << std::endl << StreamType::Both;
    out << "Введите долготу (\"ГРАДУСЫ МИНУТЫ СЕКУНДЫ\"): L1 = ";
    std::cin >> L1;
    out << StreamType::FileOnly << L1 << std::endl << StreamType::Both;
    out << "Введите азимут (\"ГРАДУСЫ МИНУТЫ СЕКУНДЫ\"): A1 = ";
    std::cin >> A1;
    out << StreamType::FileOnly << A1 << std::endl << StreamType::Both;
    out << "Введите длину (м): s = ";
    std::cin >> s;
    out << StreamType::FileOnly << s << std::endl << StreamType::Both;

    out << "Выберите способ подсчета коэффициентов A, B и C." << std::endl;
    out << "0 - экономизация и отбрасывание незначительных коэффициентов (как у Морозова)" << std::endl;
    out << "1 - экономизация без отбрасывания" << std::endl;
    out << "2 - сумма ряда до 100 члена" << std::endl;

    int w;
    out << "Введите цифру от 0 до 2: ";
    std::cin >> w;
    out << StreamType::FileOnly << w << std::endl << StreamType::Both;
    Way way = static_cast<Way>(w);

    Point point{to_radians(B1), to_radians(L1), to_radians(A1)};

    state.point1 = point;
    state.s = s;
    this->way = way;
}

void DirectProblem::calc_problem(DirectState& state)
{
    state.W1 = std::sqrt(1 - this->ellipsoid.e2 * std::pow(std::sin(state.point1.B), 2));
    state.sin_u1 = std::sin(state.point1.B) / state.W1 * std::sqrt(1 - this->ellipsoid.e2);
    state.cos_u1 = std::cos(state.point1.B) / state.W1;

    state.sin_A1 = std::sin(state.point1.A);
    state.cos_A1 = std::cos(state.point1.A);

    state.sin_A0 = state.cos_u1 * state.sin_A1;
    state.ctg_sigma1 = state.cos_u1 * state.cos_A1 / state.sin_u1;
    state.sin_2sigma1 = 2 / (std::pow(state.ctg_sigma1, 2) + 1) * state.ctg_sigma1;
    state.cos_2sigma1 = (std::pow(state.ctg_sigma1, 2) - 1) / (std::pow(state.ctg_sigma1, 2) + 1);

    state.cos2_A0 = 1 - std::pow(state.sin_A0, 2);

    state.coeffs = ProblemMath::calc_coefficients(this->way, this->ellipsoid, state.cos2_A0);

    state.alpha = (this->ellipsoid.e2 / 2 + std::pow(this->ellipsoid.e2, 2) / 8 +
                    std::pow(this->ellipsoid.e2, 3) / 16) -
                    (std::pow(this->ellipsoid.e2, 2) / 16 + std::pow(this->ellipsoid.e2, 3) / 16) * state.cos2_A0 +
                    3.0 / 128 * std::pow(this->ellipsoid.e2, 3) * std::pow(state.cos2_A0, 2);
    state.beta = (std::pow(this->ellipsoid.e2, 2) / 32 + std::pow(this->ellipsoid.e2, 3) / 32) * state.cos2_A0 -
                    std::pow(this->ellipsoid.e2, 3) / 64 * std::pow(state.cos2_A0, 2);

    state.sigma0 = (state.s - (state.coeffs.B + state.coeffs.C * state.cos_2sigma1) * state.sin_2sigma1) / state.coeffs.A;
    state.sin_2sigma10 = state.sin_2sigma1 * std::cos(2 * state.sigma0) + state.cos_2sigma1 * std::sin(2 * state.sigma0);
    state.cos_2sigma10 = state.cos_2sigma1 * std::cos(2 * state.sigma0) - state.sin_2sigma1 * std::sin(2 * state.sigma0);
    state.sigma = state.sigma0 + (state.coeffs.B + 5 * state.coeffs.C * state.cos_2sigma10)
                * state.sin_2sigma10 / state.coeffs.A;

    state.delta = (state.alpha * state.sigma + state.beta * (state.sin_2sigma10 - state.sin_2sigma1)) * state.sin_A0;

    state.sin_sigma = std::sin(state.sigma);
    state.cos_sigma = std::cos(state.sigma);
    double sin_u2 = state.sin_u1 * state.cos_sigma + state.cos_u1 * state.cos_A1 * state.sin_sigma;
    state.point2.B = std::atan(sin_u2 / std::sqrt(1 - this->ellipsoid.e2) / std::sqrt(1 - std::pow(sin_u2, 2)));

    state.lambda = std::atan(state.sin_A1 * state.sin_sigma / 
                    (state.cos_u1 * state.cos_sigma - state.sin_u1 * state.sin_sigma * state.cos_A1));

    if (state.sin_A1 >= 0)
    {
        if (std::tan(state.lambda) >= 0)
        {
            state.lambda = std::abs(state.lambda);
        }
        else
        {
            state.lambda = M_PI - std::abs(state.lambda);
        }
    } else
    {
        if (std::tan(state.lambda) >= 0)
        {
            state.lambda = std::abs(state.lambda) - M_PI;
        }
        else
        {
            state.lambda = -std::abs(state.lambda);
        }
    }

    state.point2.L = state.point1.L + state.lambda - state.delta;

    state.point2.A = std::atan(state.cos_u1 * state.sin_A1 / 
        (state.cos_u1 * state.cos_sigma * state.cos_A1 - state.sin_u1 * state.sin_sigma));

    if (state.sin_A1 < 0)
    {
        if (std::tan(state.point2.A) >= 0)
        {
            state.point2.A= std::abs(state.point2.A);
        }
        else
        {
            state.point2.A = M_PI - std::abs(state.point2.A);
        }
    } else
    {
        if (std::tan(state.point2.A) >= 0)
        {
            state.point2.A = std::abs(state.point2.A) + M_PI;
        }
        else
        {
            state.point2.A = 2 * M_PI - std::abs(state.point2.A);
        }
    }
}

void DirectProblem::output_result(DirectState& state)
{
    out << "Выходные данные." << std::endl;
    out << "W1 = " << state.W1 << std::endl;
    out << "sin(u1) = " << state.sin_u1 << std::endl;
    out << "cos(u1) = " << state.cos_u1 << std::endl;
    out << "sin(A1) = " << state.sin_A1 << std::endl;
    out << "cos(A1) = " << state.cos_A1 << std::endl;
    out << "sin(A0) = " << state.sin_A0 << std::endl;
    out << "(cos(A0))^2 = " << state.cos2_A0 << std::endl;
    out << "ctg(sigma1) = " << state.ctg_sigma1 << std::endl;
    out << "sin(2sigma1) = " << state.sin_2sigma1 << std::endl;
    out << "cos(2sigma1) = " << state.cos_2sigma1 << std::endl;
    out << "A = " << state.coeffs.A << std::endl;
    out << "B = " << state.coeffs.B << std::endl;
    out << "C = " << state.coeffs.C << std::endl;

    out << "alpha = " << to_seconds(state.alpha) << "''";
#ifdef DEBUG
    out << " = " << state.alpha << " рад ";
#endif
    out << std::endl;

    out << "beta = " << to_seconds(state.beta) << "''";
#ifdef DEBUG
    out << " = " << state.beta << " рад ";
#endif
    out << std::endl;

    out << "sigma0 = " << state.sigma0 << std::endl;
    out << "sin(2(sigma1 + sigma0)) = " << state.sin_2sigma10 << std::endl;
    out << "cos(2(sigma1 + sigma0)) = " << state.cos_2sigma10 << std::endl;
    out << "sigma = " << state.sigma << std::endl;

    out << "delta = " << to_seconds(state.delta) << "''";
#ifdef DEBUG
    out << " = " << state.delta << " рад ";
#endif
    out << std::endl;

    out << "sin(sigma) = " << state.sin_sigma << std::endl;
    out << "cos(sigma) = " << state.cos_sigma << std::endl;

    out << "lambda = " << to_degrees(state.lambda);
#ifdef DEBUG
    out << " = " << state.lambda << " рад ";
#endif
    out << std::endl;

#ifdef DEBUG
    out << PointPrinter(state.point2, 2, true) << std::endl;
#else
    out << PointPrinter(state.point2, 2) << std::endl;
#endif
    out << StreamType::FileOnly << std::endl << StreamType::Both;
}

void DirectProblem::solving_problem()
{
    DirectState state;
    input_data(state);
    calc_problem(state);
    output_result(state);
}