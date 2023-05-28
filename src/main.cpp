#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "./include/Solver_2d_Laplace_eq.hpp"
#include <cmath>
#include "../../computational_methods_labs/structures/linalg/Matrix_n_Vector.hpp"
#include "./include/eigen-3.4.0/Eigen/Sparse"




template<typename T>
class Element{
    public:

    std::vector<T> x;
    std::vector<T> y;
    std::vector<std::size_t> poly;

};


template<typename T> 
void stiffness_matrix(Class_2d_Laplace_equation<T> laplace_eq){

    typedef Eigen::SparseMatrix<T> SpMat; // declares a column-major sparse matrix type of T
    typedef Eigen::Triplet<T> Tr;



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
    reduced_sys_size = laplace_eq._nodes.size() - laplace_eq._dirichlet_upper_boundary_nodes.size() - laplace_eq._dirichlet_lower_boundary_nodes.size() - laplace_eq._left_boundary_nodes.size();

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
    fout.close();

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

    // Vector<T> tmp;
    // tmp = laplace_eq._ind_left_boundary_nodes;
    // std::cout << Vector<std::size_t> (laplace_eq._ind_left_boundary_nodes) << "\n";
    // std::cout << Vector<std::size_t> (laplace_eq._ind_right_boundary_nodes) << "\n";

    std::vector<std::pair<std::size_t, std::size_t>> ind_left_n_right_boundary_correlation;
    
    for(std::size_t i = 0; i < laplace_eq._left_boundary_nodes.size(); ++i){
        for(std::size_t j = 0; j < laplace_eq._right_boundary_nodes.size(); ++j){
            if(laplace_eq._left_boundary_nodes[i][1] == laplace_eq._right_boundary_nodes[j][1]){
                ind_left_n_right_boundary_correlation.push_back(std::pair<std::size_t, std::size_t> (laplace_eq._ind_left_boundary_nodes[i], laplace_eq._ind_right_boundary_nodes[j]));
                //std:: cout << laplace_eq._ind_left_boundary_nodes[i] << "\t" << laplace_eq._ind_right_boundary_nodes[j] << "\n";
            }
        }
    }

    for(std::size_t i = 0; i < full_matrix.get_rows(); ++i){
        for(auto it = ind_left_n_right_boundary_correlation.begin(); it != ind_left_n_right_boundary_correlation.end(); ++it){
            full_matrix[i][(*it).first] += full_matrix[i][(*it).second];
        }
    }

    std::vector<std::size_t> nodes_ind_4_reduced_matrix;
    
    auto check_it_ind = std::find(laplace_eq._ind_dirichlet_lower_boundary_nodes.begin(), laplace_eq._ind_dirichlet_lower_boundary_nodes.end(), 0);
    
    for(std::size_t i = 0; i < laplace_eq._nodes.size(); ++i){
        check_it_ind = std::find(laplace_eq._ind_dirichlet_lower_boundary_nodes.begin(), laplace_eq._ind_dirichlet_lower_boundary_nodes.end(), i);
        
        if (check_it_ind == laplace_eq._ind_dirichlet_lower_boundary_nodes.end()){

            check_it_ind = std::find(laplace_eq._ind_dirichlet_upper_boundary_nodes.begin(), laplace_eq._ind_dirichlet_upper_boundary_nodes.end(), i);
            if (check_it_ind == laplace_eq._ind_dirichlet_upper_boundary_nodes.end()){

                check_it_ind = std::find(laplace_eq._ind_right_boundary_nodes.begin(), laplace_eq._ind_right_boundary_nodes.end(), i);
                if (check_it_ind == laplace_eq._ind_right_boundary_nodes.end()){
                    nodes_ind_4_reduced_matrix.push_back(i);
                }
            }
        }
    }  

    // std::cout << Vector<std::size_t> (laplace_eq._ind_left_boundary_nodes) << "\n";
    // std::cout << Vector<std::size_t> (nodes_ind_4_reduced_matrix);
    
    for(std::size_t i = 0; i < reduced_sys_size; ++i){
        for(std::size_t j = 0; j < reduced_sys_size; ++j){
            reduced_matrix[i][j] = full_matrix[nodes_ind_4_reduced_matrix[i]][nodes_ind_4_reduced_matrix[j]];
        }
        reduced_b[i] = full_b[nodes_ind_4_reduced_matrix[i]];
    }

    fout.open("./reduced_matrix.txt");

    for(std::size_t i = 0; i < reduced_matrix.get_rows(); ++i){
        for(std::size_t j = 0; j < reduced_matrix.get_cols(); ++j){            
            fout << reduced_matrix[i][j] << "\t";
        }
        fout << "\n";
    }
    fout << "\n";
    fout.close();


    std::vector<Tr> tripletList;
    tripletList.reserve(reduced_sys_size);
    Eigen::VectorX<T> sparse_solution(reduced_sys_size), b_sp(reduced_sys_size);

    
    T mat_el = 0.;
    T vec_el = 0.;

    T tol = 1e-9;

    for(std::size_t i = 0; i < reduced_matrix.get_rows(); ++i){
        for(std::size_t j = 0; j < reduced_matrix.get_cols(); ++j){ 
            mat_el = reduced_matrix[i][j];
            if (fabs(mat_el) > tol)
                tripletList.push_back(Tr(i, j, mat_el));
        }
        vec_el = reduced_b[i];
        if (fabs(vec_el) > tol)
            b_sp[i] = vec_el;
    }

    SpMat sparse_matrix(reduced_sys_size, reduced_sys_size);
    sparse_matrix.setFromTriplets(tripletList.begin(), tripletList.end());

    Eigen::BiCGSTAB<SpMat> solver;
    solver.compute(sparse_matrix);
    sparse_solution = solver.solve(b_sp);

    for(std::size_t i = 0; i < reduced_sys_size; ++i){
        std::cout << i << ": " << sparse_solution[i] << "\n";
    }

}

int main(int args, char **argv){

    typedef double T;
    Class_2d_Laplace_equation<T> laplace_eq;
    laplace_eq.DEFAULT_TEST();

    stiffness_matrix(laplace_eq);
    
    return 0;
}
