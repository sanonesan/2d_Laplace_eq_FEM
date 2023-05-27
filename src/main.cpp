#include <iostream>
#include <fstream>
#include <vector>
#include "./include/Solver_2d_Laplace_eq.hpp"

#include "../../computational_methods_labs/structures/linalg/Matrix_n_Vector.hpp"

template<typename T>
class Element{
    public:

    std::vector<T> x;
    std::vector<T> y;
    std::vector<std::size_t> poly;

};


template<typename T> 
void stiffness_matrix(Class_2d_Laplace_equation<T> laplace_eq){


    auto ind_shift = [](std::size_t index){
        if(index == 3){
            return (std::size_t)0;
        }
        if(index == 4){
            return (std::size_t)1;
        }
        
        return index;
    };


    std::ofstream fout;
    fout.open("./stiffness_matrix_elements_of_each_triangle.txt");
    fout << std::scientific;
    fout << std::setprecision(8);


    Matrix<T> full_matrix(laplace_eq._nodes.size(), laplace_eq._nodes.size());
    Vector<T> full_b(laplace_eq._nodes.size() - laplace_eq._boundary_nodes.size());
    
    std::size_t reduced_sys_size = 0;
    reduced_sys_size = laplace_eq._nodes.size() - laplace_eq._dirichlet_upper_boundary_nodes.size() - laplace_eq._dirichlet_lower_boundary_nodes.size();

    Matrix<T> reduced_matrix(reduced_sys_size, reduced_sys_size);
    Vector<T> reduced_b(reduced_sys_size);

    //check for first element
    Matrix<T> A(3, 3);    
    Element<T> el;
    T S_el = 0.;

    for(std::size_t k = 0; k < laplace_eq._polygons.size(); ++k){


        for(std::size_t i = 0; i < 3; ++i){
            el.x.push_back( laplace_eq._nodes[laplace_eq._polygons[k][i]][0] );
            el.y.push_back( laplace_eq._nodes[laplace_eq._polygons[k][i]][1] );
        }
        el.poly = laplace_eq._polygons[k];

        //std::cout << el.poly[0] << "\t" << el.poly[1] << "\t" << el.poly[2] << "\n";
        // el.print();
        S_el = ((el.x[1] - el.x[0]) * (el.y[2] - el.y[0]) - (el.x[2] - el.x[0]) * (el.y[1] - el.y[0])) / 2;
        // std:: cout << S_el << "\n";
        S_el = 1. / 4. / S_el;
        // std:: cout << S_el << "\n";

        

        for(std::size_t i = 0; i < 3; ++i){
            for(std::size_t j = 0; j < 3; ++j){            
                A[i][j] = S_el * ( (el.y[ind_shift(i+1)] - el.y[ind_shift(i+2)]) * (el.y[ind_shift(j+1)] - el.y[ind_shift(j+2)]) + (el.x[ind_shift(i+2)] - el.x[ind_shift(i+1)]) * (el.x[ind_shift(j+2)] - el.x[ind_shift(j+1)]));
                fout << A[i][j] << "\t";
                full_matrix[laplace_eq._polygons[k][i]][laplace_eq._polygons[k][j]] += A[i][j];
            }
            fout << "\n";
        }
        fout << "\n";

    }
    fout.close();

    fout.open("./full_matrix.txt");

    for(std::size_t i = 0; i < full_matrix.get_rows(); ++i){
        for(std::size_t j = 0; j < full_matrix.get_cols(); ++j){            
            fout << full_matrix[i][j] << "\t";
        }
        fout << "\n";
    }
    fout << "\n";
    T sum = 0.;
    for(std::size_t i = 0; i < full_matrix.get_rows(); ++i){
        for(auto it = laplace_eq._ind_dirichlet_lower_boundary_nodes.begin(); it != laplace_eq._ind_dirichlet_lower_boundary_nodes.end(); ++it){
            full_b[i] -= full_matrix[i][*it] * laplace_eq._dirichlet_lower_boundary_condition(laplace_eq._nodes[*it][0], laplace_eq._nodes[*it][1]);
        }
    }
    for(std::size_t i = 0; i < full_matrix.get_rows(); ++i){
        for(auto it = laplace_eq._ind_dirichlet_upper_boundary_nodes.begin(); it != laplace_eq._ind_dirichlet_upper_boundary_nodes.end(); ++it){
            full_b[i] -= full_matrix[i][*it] * laplace_eq._dirichlet_upper_boundary_condition(laplace_eq._nodes[*it][0], laplace_eq._nodes[*it][1]);
        }
    }

    // std::cout << "dagdasdg";
    
    // for(auto j = 0; j < full_b.size(); ++j)
    //     std::cout << full_b[j] << "\t";
    
    //std::cout << laplace_eq._polygons.size();
    //std::cout << laplace_eq._nodes.size() - laplace_eq._dirichlet_upper_boundary_nodes.size() - laplace_eq._dirichlet_lower_boundary_nodes.size();
    //A.print();
}

int main(int args, char **argv){

    typedef double T;
    Class_2d_Laplace_equation<T> laplace_eq;
    laplace_eq.DEFAULT_TEST();

    stiffness_matrix(laplace_eq);
    
    return 0;
}
