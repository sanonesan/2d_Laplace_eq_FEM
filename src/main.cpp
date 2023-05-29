#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "./include/Solver_2d_Laplace_eq.hpp"
#include "./include/Tests_for_domains/Test_domain_1_sin.hpp"
#include "./include/Tests_for_domains/Test_domain_2_rectangle.hpp"


int main(int args, char **argv){

    typedef double T;
    
    Solver_2d_Laplace_eq<T> solver;
    solver.tol = 1e-16;
    solver.notifications = true;

    Test_domain_1_sin<T> test1;
    // solver.output_folder = "../output/domain_1/mesh001/";
    // solver.file_name = test1._name;
    // solver.solve_eq(test1);

    // test1.Set_mesh_01();
    // solver.output_folder = "../output/domain_1/mesh01/";
    // solver.file_name = test1._name;
    // solver.solve_eq(test1);

    // test1.Set_mesh_005();
    // solver.output_folder = "../output/domain_1/mesh005/";
    // solver.file_name = test1._name;
    // solver.solve_eq(test1);

    test1.Set_mesh_0001();
    solver.output_folder = "../output/domain_1/mesh0001/";
    solver.file_name = test1._name;
    solver.solve_eq(test1);

    // Test_domain_2_rectangle<T> test2;
    // solver.output_folder = "../output/domain_2/mesh001/";
    // solver.file_name = test2._name;
    // solver.solve_eq(test2);
    
    return 0;
}
