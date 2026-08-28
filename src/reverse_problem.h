#ifndef REVERSE_PROBLEM_H
#define REVERSE_PROBLEM_H

#include "primitives.h"

struct ReverseState
{
    //input
    Point point1; // !without A
    Point point2; // !without A

    // intermediate
    double l;
    double sin_u1;
    double cos_u1;
    double sin_u2;
    double cos_u2;
    double a1;
    double a2;
    double b1;
    double b2;

    std::vector<double> sin_lambda;
    std::vector<double> cos_lambda;
    std::vector<double> p;
    std::vector<double> q;
    std::vector<double> sin_A1;
    std::vector<double> cos_A1;
    std::vector<double> sin_sigma;
    std::vector<double> cos_sigma;
    std::vector<double> sigma;
    std::vector<double> sin_A0;
    std::vector<double> x;
    std::vector<double> alpha;
    std::vector<double> beta_;
    std::vector<double> delta;

    double A;
    double B_;
    double C_;
    double y;

    // output
    double s;
    std::vector<double> A1;
    double A2;
};

struct ReverseProblem : Problem
{
public:
    void solving_problem() override;

protected:
    void input_ellipsoid() override;

private:
    void input_data(ReverseState&);
    void calc_problem(ReverseState&);
    void output_result(ReverseState&);

    void output_header(int size_col, int num_cols);
    void output_row(const std::string& name, const std::vector<std::string>& data, 
                    int size_col, int num_cols);
};

#endif // REVERSE_PROBLEM_H