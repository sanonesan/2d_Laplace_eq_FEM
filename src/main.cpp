#include <iostream>
#include <fstream>
#include <vector>
#include "./include/Solver_2d_Laplace_eq.hpp"

int main(int args, char **argv){

    typedef double T;
    Class_2d_Laplace_equation<T> laplace_eq;
    laplace_eq.DEFAULT_TEST();


    
    return 0;
}
