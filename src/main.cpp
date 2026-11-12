#include <iostream>
#include <limits>

#include "direct_problem.h"
#include "reverse_problem.h"
#include "primitives.h"

void help()
{
    std::cout << "Доступные команды: \\help" << std::endl;
    std::cout << "Прямая задача: \\direct" << std::endl;
    std::cout << "Обратная задача: \\reverse" << std::endl;
    std::cout << "Выход: \\quit" << std::endl;
}

int main()
{
    std::cout.precision(std::numeric_limits<double>::max_digits10);

    help();

    std::string command;
    while (std::cin >> command && command != "\\quit")
    {
        if (command == "\\direct")
        {
            DirectProblem problem;
            problem.init();
            problem.solving_problem();
        } else if (command == "\\reverse")
        {
            ReverseProblem problem;
            problem.init();
            problem.solving_problem();
        } else
        {
            help();
        }
    }
    
    return 0;
}