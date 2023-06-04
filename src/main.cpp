#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "./include/Solver_2d_Laplace_eq.hpp"
#include "./include/Tests_for_domains/Test_domain_1_sin.hpp"
#include "./include/Tests_for_domains/Test_domain_1_1_sin.hpp"
#include "./include/Tests_for_domains/Test_domain_1_1_1_sin.hpp"
#include "./include/Tests_for_domains/Test_domain_2_rectangle.hpp"
#include "./include/Tests_for_domains/Test_domain_2_rectangle_dirichlet_only.hpp"
#include "./include/Tests_for_domains/Test_domain_3_sin_sin.hpp"
#include "./include/Tests_for_domains/Test_domain_4_abs_x.hpp"

#include <time.h>

int main(int args, char **argv){

    clock_t t_start, t_end;


    typedef double T;
    
    Solver_2d_Laplace_eq<T> solver;
    solver.tol = 1e-16;
    solver.notifications = true;




    t_start  = clock();

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

    test1.Set_mesh_001();
    solver.output_folder = "../output/domain_1/mesh001/";
    solver.file_name = test1._name;
    solver.solve_eq(test1);

    solver.output_folder = "../output/domain_1/mesh001_1/";
    solver.file_name = test1._name;
    solver.solve_eq_mod(test1);

    

    // Test_domain_2_rectangle<T> test2;
    // solver.output_folder = "../output/domain_2/mesh001/";
    // solver.file_name = test2._name;
    // solver.solve_eq(test2);

    // Test_domain_3_sin_sin<T> test3;
    // solver.output_folder = "../output/domain_3/mesh001/";
    // solver.file_name = test3._name;
    // solver.solve_eq(test3);

    // Test_domain_2_rectangle_dirichlet_only<T> test4;
    // test4.Set_mesh_001();
    // solver.output_folder = "../output/domain_2_extra/mesh001/";
    // solver.file_name = test4._name;
    // solver.solve_eq_testing(test4);

    // test4.Set_mesh_0001();
    // solver.output_folder = "../output/domain_2_extra/mesh0001/";
    // solver.file_name = test4._name;
    // solver.solve_eq_mod_THREAD(test4);

    // Test_domain_1_1_sin<T> test1_1;
    // solver.output_folder = "../output/domain_1_1/mesh001/";
    // solver.file_name = test1_1._name;
    // solver.solve_eq(test1_1);

    // Test_domain_1_1_1_sin<T> test1_1_1;
    // solver.output_folder = "../output/domain_1_1_1/mesh001/";
    // solver.file_name = test1_1_1._name;
    // solver.solve_eq(test1_1_1);


    Test_domain_4_abs_x<T> test5;
    solver.output_folder = "../output/domain_4/mesh001_1/";
    solver.file_name = test5._name;
    solver.solve_eq_mod(test5);
    solver.output_folder = "../output/domain_4/mesh001/";
    solver.solve_eq(test5);

    
    t_end = clock();

    std::cout << ("Time taken: %.2fs\n", (double)(t_end - t_start)/CLOCKS_PER_SEC) << "\n";

    return 0;
}
