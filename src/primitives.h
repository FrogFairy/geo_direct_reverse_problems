#ifndef PRIMITIVES_H
#define PRIMITIVES_H

#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
#include <vector>
#include <limits>

struct Ellipsoid
{
    double b;
    double e2;
    double e_2;
};

struct Point
{
    double B; // radians
    double L; // radians
    double A; // radians
};

struct PointPrinter
{
    const Point* point;
    int index;
    bool show_radians;

    PointPrinter(const Point& p, int ind = 1, bool show_rad = false) 
        : point(&p), index(ind), show_radians(show_rad) {}
};

struct Degrees
{
    int degrees;
    int minutes;
    double seconds;
};

enum Way
{
    Simplified = 0,
    Standart = 1,
    // Accurate = 2
};

struct Coefficients
{
    double A;
    double B;
    double C;
};

enum class StreamType
{
    FileOnly,
    ConsoleOnly,
    Both
};

struct StreamSplitter {
private:
    std::ostream& console;
    std::unique_ptr<std::ofstream> file;
    std::string filename;
    bool writeToFile = false;
    bool writeToConsole = true;

    std::string Cp1125ToCp1251(const std::string& cp1125Text);
    void SaveAsWindows1251(const std::stringstream& cp1125Text);
    
public:
    StreamSplitter() : console(std::cout) {}
    
    void enableFileLogging(const std::string& filename);
    
    // void disableFileLogging();

    StreamSplitter& operator<<(StreamType);
    
    template<typename T>
    StreamSplitter& operator<<(const T& value) {
        if (writeToConsole)
        {
            console << value;
        }
        if (writeToFile && file && file->is_open()) {
            std::stringstream ss;
            ss << value;
            SaveAsWindows1251(ss);
        }
        return *this;
    }

    StreamSplitter& operator<<(std::ostream& (*manip)(std::ostream&));
};

struct Problem
{
public:
    Problem() = default;
    virtual ~Problem() = default;

    virtual void solving_problem() = 0;

    void set_output();
    void init();
protected:
    Ellipsoid ellipsoid;
    Way way;
    StreamSplitter out;

    virtual void input_ellipsoid();
};

namespace ProblemMath
{
    Coefficients calc_coefficients(const Way& way, const Ellipsoid& ellipsoid, double cos2_A0);
    Coefficients calc_coefficients_simplified(const Ellipsoid& ellipsoid, double cos2_A0);
    Coefficients calc_coefficients_standart(const Ellipsoid& ellipsoid, double cos2_A0);
    // Coefficients calc_coefficients_accurate(const Ellipsoid& ellipsoid, double cos2_A0);
};

double to_radians(Degrees);

double inline to_seconds(double radians)
{
    return radians / M_PI * 180 * 3600;
}
std::vector<double> to_seconds(const std::vector<double>&);

Degrees to_degrees(double);
std::vector<Degrees> to_degrees(std::vector<double>&);

std::ostream& operator<<(std::ostream&, const Degrees&);
std::ostream& operator<<(std::ostream&, const PointPrinter&);

std::istream& operator>>(std::istream&, Degrees&);

std::string alignment(const std::string& data, int width, char align='^', char filling=' ');
std::string to_string_with_precision(double value, int precision);

std::vector<std::string> data_to_string(const std::vector<double>&);
std::vector<std::string> data_to_string(const std::vector<Degrees>&);

#endif // PRIMITIVES_H