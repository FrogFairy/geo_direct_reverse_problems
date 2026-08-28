#include "reverse_problem.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <utility>

void ReverseProblem::input_ellipsoid()
{
    out << "Обратная задача." << std::endl;
    Problem::input_ellipsoid();
}

void ReverseProblem::input_data(ReverseState& state) 
{
    Degrees B1, B2, L1, L2;
    
    out << "Обратная задача." << std::endl;
    out << "Входные данные." << std::endl;
    out << "Введите широту начальной точки (\"ГРАДУСЫ МИНУТЫ СЕКУНДЫ\"): B1 = ";
    std::cin >> B1;
    out << StreamType::FileOnly << B1 << std::endl << StreamType::Both;
    out << "Введите широту конечной точки (\"ГРАДУСЫ МИНУТЫ СЕКУНДЫ\"): B2 = ";
    std::cin >> B2;
    out << StreamType::FileOnly << B2 << std::endl << StreamType::Both;
    out << "Введите долготу начальной точки (\"ГРАДУСЫ МИНУТЫ СЕКУНДЫ\"): L1 = ";
    std::cin >> L1;
    out << StreamType::FileOnly << L1 << std::endl << StreamType::Both;
    out << "Введите долготу конечной точки (\"ГРАДУСЫ МИНУТЫ СЕКУНДЫ\"): L2 = ";
    std::cin >> L2;
    out << StreamType::FileOnly << L2 << std::endl << StreamType::Both;

    out << "Выберите способ подсчета коэффициентов A, B и C." << std::endl;
    out << "0 - экономизация и отбрасывание незначительных коэффициентов (как у Морозова)" << std::endl;
    out << "1 - экономизация без отбрасывания" << std::endl;
    out << "2 - сумма ряда до 100 члена" << std::endl;

    int w;
    out << "Введите цифру от 0 до 2: ";
    std::cin >> w;
    out << StreamType::FileOnly << w << std::endl << StreamType::Both;
    Way way = static_cast<Way>(w);

    Point point1{to_radians(B1), to_radians(L1), 0};
    Point point2{to_radians(B2), to_radians(L2), 0};

    state.point1 = point1;
    state.point2 = point2;
    this->way = way;
}

void ReverseProblem::calc_problem(ReverseState& state)
{
    double sin_B1 = std::sin(state.point1.B);
    double sin_B2 = std::sin(state.point2.B);
    double cos_B1 = std::cos(state.point1.B);
    double cos_B2 = std::cos(state.point2.B);

    double W1 = std::sqrt(1 - this->ellipsoid.e2 * std::pow(sin_B1, 2));
    double W2 = std::sqrt(1 - this->ellipsoid.e2 * std::pow(sin_B2, 2));

    state.sin_u1 = sin_B1 / W1 * std::sqrt(1 - this->ellipsoid.e2);
    state.sin_u2 = sin_B2 / W2 * std::sqrt(1 - this->ellipsoid.e2);

    state.cos_u1 = cos_B1 / W1;
    state.cos_u2 = cos_B2 / W2;

    state.l = state.point2.L - state.point1.L;
    state.a1 = state.sin_u1 * state.sin_u2;
    state.a2 = state.cos_u1 * state.cos_u2;
    state.b1 = state.cos_u1 * state.sin_u2;
    state.b2 = state.sin_u1 * state.cos_u2;

    double delta = 0;
    double eps = 1e-4;
    double lambda;
    do
    {
        if (!state.delta.empty())
        {
            delta = state.delta.back();
        }

        lambda = state.l + delta;
        state.sin_lambda.push_back(std::sin(lambda));
        state.cos_lambda.push_back(std::cos(lambda));
        state.p.push_back(state.cos_u2 * state.sin_lambda.back());
        state.q.push_back(state.b1 - state.b2 * state.cos_lambda.back());

        double A1 = std::atan(state.p.back() / state.q.back());
        if (state.p.back() > 0)
        {
            if (state.q.back() > 0)
            {
                state.A1.push_back(std::abs(A1));
            } else
            {
                state.A1.push_back(M_PI - std::abs(A1));
            }
        } else
        {
            if (state.q.back() > 0)
            {
                state.A1.push_back(2 * M_PI - std::abs(A1));
            } else
            {
                state.A1.push_back(M_PI + std::abs(A1));
            }
        }

        state.sin_A1.push_back(std::sin(state.A1.back()));
        state.cos_A1.push_back(std::cos(state.A1.back()));

        state.sin_sigma.push_back(state.p.back() * state.sin_A1.back() + 
                                    state.q.back() * state.cos_A1.back());
        state.cos_sigma.push_back(state.a1 + state.a2 * state.cos_lambda.back());

        double sigma = std::atan(state.sin_sigma.back() / state.cos_sigma.back());
        if (state.cos_sigma.back() > 0)
        {
            state.sigma.push_back(std::abs(sigma));
        } else
        {
            state.sigma.push_back(M_PI - std::abs(sigma));
        }

        state.sin_A0.push_back(state.cos_u1 * state.sin_A1.back());
        double cos2_A0 = 1 - std::pow(state.sin_A0.back(), 2);

        state.x.push_back(2 * state.a1 - cos2_A0 * state.cos_sigma.back());

        state.alpha.push_back((this->ellipsoid.e2 / 2 + std::pow(this->ellipsoid.e2, 2) / 8 +
                    std::pow(this->ellipsoid.e2, 3) / 16) -
                    (std::pow(this->ellipsoid.e2, 2) / 16 + std::pow(this->ellipsoid.e2, 3) / 16) * cos2_A0 +
                    3.0 / 128 * std::pow(this->ellipsoid.e2, 3) * std::pow(cos2_A0, 2));
        
        double beta = (std::pow(this->ellipsoid.e2, 2) / 32 + std::pow(this->ellipsoid.e2, 3) / 32) * cos2_A0 -
                    std::pow(this->ellipsoid.e2, 3) / 64 * std::pow(cos2_A0, 2);
        state.beta_.push_back(2 * beta / cos2_A0);

        state.delta.push_back((state.alpha.back() * state.sigma.back() -
                                state.beta_.back() * state.x.back() * state.sin_sigma.back())
                                * state.sin_A0.back());
    } while (std::abs(to_seconds(state.delta.back()) - to_seconds(delta)) >= eps);

    double cos2_A0 = 1 - std::pow(-0.705597732, 2);

    Coefficients coeffs = ProblemMath::calc_coefficients(this->way, this->ellipsoid, cos2_A0);
    state.A = coeffs.A;
    state.B_ = 2 / cos2_A0 * coeffs.B;
    state.C_ = 2 / std::pow(cos2_A0, 2) * coeffs.C;
    out << cos2_A0 << " " << coeffs.C << std::endl;

    state.y = (std::pow(cos2_A0, 2) - 2 * std::pow(state.x.back(), 2)) * state.cos_sigma.back();
    state.s = state.A * state.sigma.back() + (state.B_ * state.x.back() + state.C_ * state.y) * state.sin_sigma.back();

    double p = state.cos_u1 * state.sin_lambda.back();
    double q = state.b1 * state.cos_lambda.back() - state.b2;
    double A2 = std::atan(p / q);
    if (p > 0)
    {
        if (q > 0)
        {
            state.A2 = M_PI + std::abs(A2);
        } else
        {
            state.A2 = 2 * M_PI - std::abs(A2);
        }
    } else if (p < 0)
    {
        if (q > 0)
        {
            state.A2 = M_PI - std::abs(A2);
        } else
        {
            state.A2 = std::abs(A2);
        }
    }
}

void ReverseProblem::output_header(int size_col, int num_cols)
{
    out << alignment("Величина", size_col, '<') << "|"
                << alignment("Приближения", size_col * num_cols + num_cols - 1) << "|" << std::endl;
    out << std::string(size_col + 1 + size_col * num_cols + num_cols, '-') << std::endl;
    out << std::string(size_col, ' ') << "|";
    for (int i = 1; i <= num_cols; ++i)
    {
        out << alignment(std::to_string(i), size_col) << "|";
    }
    out << std::endl;
    out << std::string(size_col + 1 + size_col * num_cols + num_cols, '-') << std::endl;
}

void ReverseProblem::output_row(const std::string& name, const std::vector<std::string>& data, 
                                int size_col, int num_cols)
{
    out << alignment(name, size_col, '<') << "|";
    for (int i = 0; i < num_cols; ++i)
    {
        out << alignment(data[i], size_col) << "|";
    }
    out << std::endl;
}

void ReverseProblem::output_result(ReverseState& state)
{
    out << "Выходные данные." << std::endl;
    out << "l = " << state.l << std::endl;
    out << "sin(u1) = " << state.sin_u1 << std::endl;
    out << "cos(u1) = " << state.cos_u1 << std::endl;
    out << "sin(u2) = " << state.sin_u2 << std::endl;
    out << "cos(u2) = " << state.cos_u2 << std::endl;
    out << "a1 = " << state.a1 << std::endl;
    out << "a2 = " << state.a2 << std::endl;
    out << "b1 = " << state.b1 << std::endl;
    out << "b2 = " << state.b2 << std::endl;

    int size_col = 25;
    int num_cols = state.delta.size();
    
    output_header(size_col, num_cols);
    output_row("sin(lambda)", data_to_string(state.sin_lambda), size_col, num_cols);
    output_row("cos(lambda)", data_to_string(state.cos_lambda), size_col, num_cols);
    output_row("p", data_to_string(state.p), size_col, num_cols);
    output_row("q", data_to_string(state.q), size_col, num_cols);
    output_row("A1", data_to_string(to_degrees(state.A1)), size_col, num_cols);
#ifdef DEBUG
    output_row("A1", data_to_string(state.A1), size_col, num_cols);
#endif
    output_row("sin(A1)", data_to_string(state.sin_A1), size_col, num_cols);
    output_row("cos(A1)", data_to_string(state.cos_A1), size_col, num_cols);
    output_row("sin(sigma)", data_to_string(state.sin_sigma), size_col, num_cols);
    output_row("cos(sigma)", data_to_string(state.cos_sigma), size_col, num_cols);
    output_row("sigma", data_to_string(state.sigma), size_col, num_cols);
    output_row("sin(A0)", data_to_string(state.sin_A0), size_col, num_cols);
    output_row("x", data_to_string(state.x), size_col, num_cols);
    output_row("alpha", data_to_string(state.alpha), size_col, num_cols);
    output_row("beta'", data_to_string(state.beta_), size_col, num_cols);

    std::vector<std::string> delta = data_to_string(to_seconds(state.delta));
    for (std::string& s : delta)
    {
        s += "''";
    }
    output_row("delta", delta, size_col, num_cols);
#ifdef DEBUG
    output_row("delta", data_to_string(state.delta), size_col, num_cols);
#endif

    out << "A = " << state.A << std::endl;
    out << "B' = " << state.B_ << std::endl;
    out << "C' = " << state.C_ << std::endl;
    out << "y = " << state.y << std::endl;
    out << "s = " << state.s << std::endl;
    out << "A2 = " << to_degrees(state.A2);
#ifdef DEBUG
    out << " = " << state.A2 << " рад ";
#endif
    out << std::endl;
    out << StreamType::FileOnly << std::endl << StreamType::Both;
}


void ReverseProblem::solving_problem()
{
    ReverseState state;

    input_data(state);
    calc_problem(state);
    output_result(state);
}