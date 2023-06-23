#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <time.h>

//Import Solver
#include "./include/Solver_2d_Laplace_eq.hpp"

//Import tests
#include "./include/Tests_for_domains/Test_domain_1_sin.hpp"
#include "./include/Tests_for_domains/Test_domain_1_1_sin.hpp"
#include "./include/Tests_for_domains/Test_domain_1_1_1_sin.hpp"
#include "./include/Tests_for_domains/Test_domain_2_rectangle.hpp"
#include "./include/Tests_for_domains/Test_domain_2_rectangle_dirichlet_only.hpp"
#include "./include/Tests_for_domains/Test_domain_3_sin_sin.hpp"
#include "./include/Tests_for_domains/Test_domain_4_abs_x.hpp"


int main(int args, char **argv){

    clock_t t_start, t_end;


    typedef double T;
    
    Solver_2d_Laplace_eq<T> solver;
    solver.tol = 1e-16;
    solver.notifications = true;


    t_start  = clock();

    // Test_domain_1_sin<T> test1;

    // test1.Set_mesh_01();
    // solver.output_folder = "../output/domain_1/mesh01/";
    // solver.file_name = test1._name;
    // solver.solve_eq(test1);

    // test1.Set_mesh_005();
    // solver.output_folder = "../output/domain_1/mesh005/";
    // solver.file_name = test1._name;
    // solver.solve_eq(test1);

    // test1.Set_mesh_001();
    // solver.output_folder = "../output/domain_1/mesh001/";
    // solver.file_name = test1._name;
    // solver.solve_eq(test1);       

    // test1.Set_mesh_0001();
    // solver.output_folder = "../output/domain_1/mesh0001/";
    // solver.file_name = test1._name;
    // solver.solve_eq(test1);


    /*
        Rectangle periodic
    */

    // Test_domain_2_rectangle<T> test2;
    // test2.Set_mesh_005();    
    // solver.output_folder = "../output/domain_2/mesh005/";
    // solver.file_name = test2._name;
    // solver.solve_eq(test2);

    // test2.Set_mesh_005_calfem();    
    // solver.output_folder = "../output/domain_2/mesh005_calfem/";
    // solver.file_name = test2._name;
    // solver.solve_eq(test2);

    // test2.Set_mesh_001();    
    // solver.output_folder = "../output/domain_2/mesh001_calfem/";
    // solver.file_name = test2._name;
    // solver.solve_eq(test2);

    // test2.Set_mesh_0001();    
    // solver.output_folder = "../output/domain_2/mesh0001/";
    // solver.file_name = test2._name;
    // solver.solve_eq(test2);

    // test2.Set_mesh_00008();    
    // solver.output_folder = "../output/domain_2/mesh0008/";
    // solver.file_name = test2._name;
    // solver.solve_eq(test2);

    // test2.Set_mesh_00008_s();    
    // solver.output_folder = "../output/domain_2/mesh0008_s/";
    // solver.file_name = test2._name;
    // solver.solve_eq(test2);



    // Test_domain_3_sin_sin<T> test3;
    // solver.output_folder = "../output/domain_3/mesh001/";
    // solver.file_name = test3._name;
    // solver.solve_eq(test3);



    /*
        Dirichlet boundary tests
    */



    Test_domain_2_rectangle_dirichlet_only<T> test4;

    // test4.Set_mesh_005();
    // solver.output_folder = "../output/domain_2_extra/mesh005/";
    // solver.file_name = test4._name;
    // solver.solve_eq_dirichlet_only(test4);

    // test4.Set_mesh_005_rev_calfem();
    // solver.output_folder = "../output/domain_2_extra/mesh005_calfem/";
    // solver.file_name = test4._name;
    // solver.solve_eq_dirichlet_only(test4);

    // test4.Set_mesh_005_calfem();
    // solver.output_folder = "../output/domain_2_extra/mesh005_calfem/";
    // solver.file_name = test4._name;
    // solver.solve_eq_dirichlet_only(test4);

    // test4.Set_mesh_001();
    // solver.output_folder = "../output/domain_2_extra/mesh001/";
    // solver.file_name = test4._name;
    // solver.solve_eq_dirichlet_only(test4);

    // test4.Set_mesh_001_calfem();
    // solver.output_folder = "../output/domain_2_extra/mesh001_calfem/";
    // solver.file_name = test4._name;
    // solver.solve_eq_dirichlet_only(test4);

    // test4.Set_mesh_0005_calfem();
    // solver.output_folder = "../output/domain_2_extra/mesh0005_calfem/";
    // solver.file_name = test4._name;
    // solver.solve_eq_dirichlet_only(test4);

    // test4.Set_mesh_0001();
    // solver.output_folder = "../output/domain_2_extra/mesh0001/";
    // solver.file_name = test4._name;
    // solver.solve_eq_dirichlet_only(test4);

    // test4.Set_mesh_0001_calfem();
    // solver.output_folder = "../output/domain_2_extra/mesh0001_calfem/";
    // solver.file_name = test4._name;
    // solver.solve_eq_dirichlet_only(test4);

    // test4.Set_mesh_00005_calfem();
    // solver.output_folder = "../output/domain_2_extra/mesh00005_calfem/";
    // solver.file_name = test4._name;
    // solver.solve_eq_dirichlet_only(test4);

    // test4.Set_mesh_00008();
    // solver.output_folder = "../output/domain_2_extra/mesh00008/";
    // solver.file_name = test4._name;
    // solver.solve_eq_dirichlet_only(test4);

    test4.Set_mesh_00008_calfem();
    solver.output_folder = "../output/domain_2_extra/mesh00008_calfem/";
    solver.file_name = test4._name;
    solver.solve_eq_dirichlet_only(test4);

    // test4.Set_mesh_00008_s();
    // solver.output_folder = "../output/domain_2_extra/mesh00008_s/";
    // solver.file_name = test4._name;
    // solver.solve_eq_dirichlet_only(test4);

    




    // Test_domain_1_1_sin<T> test1_1;
    // test1_1.Set_mesh_0005_calfem();
    // solver.output_folder = "../output/domain_1_1/mesh0005_calfem/";
    // solver.file_name = test1_1._name;
    // solver.solve_eq(test1_1);

    // test1_1.Set_mesh_0005_3_in_row_calfem();
    // solver.output_folder = "../output/domain_1_1/mesh0005_3_in_row_calfem/";
    // solver.file_name = test1_1._name;
    // solver.solve_eq(test1_1);

    // test1_1.Set_mesh_0001_calfem();
    // solver.output_folder = "../output/domain_1_1/mesh0001_calfem/";
    // solver.file_name = test1_1._name;
    // solver.solve_eq(test1_1);

    // test1_1.Set_mesh_0001_3_in_row_calfem();
    // solver.output_folder = "../output/domain_1_1/mesh0001_3_in_row_calfem/";
    // solver.file_name = test1_1._name;
    // solver.solve_eq(test1_1);


    // Test_domain_1_1_1_sin<T> test1_1_1;
    // solver.output_folder = "../output/domain_1_1_1/mesh001/";
    // solver.file_name = test1_1_1._name;
    // solver.solve_eq(test1_1_1);


    // Test_domain_4_abs_x<T> test5;

    // test5.Set_mesh001_calfem();
    // solver.output_folder = "../output/domain_4/mesh001_calfem/";
    // solver.file_name = test5._name;
    // solver.solve_eq(test5);


    // solver.output_folder = "../output/domain_4/mesh001_3_in_row_calfem/";
    // solver.file_name = test5._name;
    // solver.solve_eq(test5);
    

    
    t_end = clock();

    std::cout << ("Time taken: %.2fs\n", (double)(t_end - t_start)/CLOCKS_PER_SEC) << "\n";

    return 0;
}
