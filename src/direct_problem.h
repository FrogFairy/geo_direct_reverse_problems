#ifndef DIRECT_PROBLEM_H
#define DIRECT_PROBLEM_H

#include "primitives.h"

#include <utility>

struct DirectState
{
    //input
    Point point1;
    double s;

    // intermediate
    double W1;
    double sin_u1;
    double cos_u1;
    double sin_A1;
    double cos_A1;
    double sin_A0;
    double cos2_A0;
    double ctg_sigma1;
    double sin_2sigma1;
    double cos_2sigma1;
    Coefficients coeffs;
    double alpha;
    double beta;
    double sigma0;
    double sin_2sigma10;
    double cos_2sigma10;
    double sigma;
    double delta;
    double sin_sigma;
    double cos_sigma;
    double lambda;

    // output
    Point point2;
};

struct DirectProblem : Problem
{
public:
    void solving_problem() override;
private:
    void input_data(DirectState&);
    void calc_problem(DirectState&);
    void output_result(DirectState&);
};

#endif // DIRECT_PROBLEM_H