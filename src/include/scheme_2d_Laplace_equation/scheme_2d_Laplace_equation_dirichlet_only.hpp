#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../Class_2d_Laplace_equation.hpp"
#include "../Eigen/Sparse"
#include "../structures/Matrix_n_Vector.hpp"
#include "../structures/struct_triag_element.hpp"




template<typename T>
void scheme_2d_Laplace_equation_dirichlet_only(
    const Class_2d_Laplace_equation<T>& laplace_eq, 
    const T& tol, const 
    std::string& out_path, 
    bool output_stiffness_matricies_bool = true){

    typedef Eigen::SparseMatrix<T> SpMat; // declares a column-major sparse matrix type of T
    typedef Eigen::Triplet<T> Tr; // declares a triplet (i, j, el_val) for Eigen::SparseMatrix
    
    std::ofstream fout;
    std::vector<T> solution(laplace_eq._nodes.size());


    fout << std::scientific;
    fout << std::setprecision(8);

    if (output_stiffness_matricies_bool)
        fout.open(out_path + "_stiffness_matrix_elements_for_each_triangle.txt");

    /**
     * Assembling full stiffness matrix
    */

    Matrix<T> full_stiffness_matrix(laplace_eq._nodes.size(), laplace_eq._nodes.size());
    Vector<T> full_b(laplace_eq._nodes.size());

    Matrix<T> local_stiffness_matrix(3, 3); // stiffness matrix for single triangle element    
    Triag_Element<T> el;
    T S_el = 0.;

    auto ind_shift = [](std::size_t index){
        if(index == 3){
            return (std::size_t)0;
        }
        if(index == 4){
            return (std::size_t)1;
        }
        
        return index;
    };

    for(std::size_t k = 0; k < laplace_eq._polygons.size(); ++k){

        if (output_stiffness_matricies_bool) fout << k << ":\n";

        el.x.clear();
        el.y.clear();
        el.poly.clear();
        // Reading triangle element (its coordinates && polygons in global nodes vector)
        for(std::size_t i = 0; i < 3; ++i){
            el.x.push_back( laplace_eq._nodes[laplace_eq._polygons[k][i]][0] );
            el.y.push_back( laplace_eq._nodes[laplace_eq._polygons[k][i]][1] );
        }
        el.poly = laplace_eq._polygons[k];

        // Find square of triangle element
        S_el = ((el.x[1] - el.x[0]) * (el.y[2] - el.y[0]) - (el.x[2] - el.x[0]) * (el.y[1] - el.y[0])) / 2;
        // Transformation for the next formula (for numerial integral over the triangle)
        S_el = 1 / ( 4. * S_el );
        
        for(std::size_t i = 0; i < 3; ++i){
            for(std::size_t j = i; j < 3; ++j){     
                // Integral[ \nabla \phi_i * \nabla \phi_j, over triangle]
                local_stiffness_matrix[i][j] = S_el * ( (el.y[ind_shift(i+1)] - el.y[ind_shift(i+2)]) * (el.y[ind_shift(j+1)] - el.y[ind_shift(j+2)]) + (el.x[ind_shift(i+2)] - el.x[ind_shift(i+1)]) * (el.x[ind_shift(j+2)] - el.x[ind_shift(j+1)]));
                
                //output if needed
                if (output_stiffness_matricies_bool) fout << local_stiffness_matrix[i][j] << "\t";
                // Parallel assemble of full_stiffness matrix
                full_stiffness_matrix[laplace_eq._polygons[k][i]][laplace_eq._polygons[k][j]] += local_stiffness_matrix[i][j];
                if (i != j) {
                    local_stiffness_matrix[j][i] = local_stiffness_matrix[i][j];
                    full_stiffness_matrix[laplace_eq._polygons[k][j]][laplace_eq._polygons[k][i]] += local_stiffness_matrix[i][j];
                }
            }
            if (output_stiffness_matricies_bool) fout << "\n";
        }

        if (output_stiffness_matricies_bool) fout << "\n";
    }
    fout.close();

    // output full stiffness matrix in txt
    if (output_stiffness_matricies_bool){
        fout.open(out_path + "_full_matrix.txt");
        for(std::size_t i = 0; i < full_stiffness_matrix.get_rows(); ++i){
            for(std::size_t j = 0; j < full_stiffness_matrix.get_cols(); ++j){            
                fout << full_stiffness_matrix[i][j] << "\t";
            }
            fout << "\n";
        }
        fout << "\n";
        fout.close();
    }


    /**
     * Applying Dirichlet boundary conditions
     * 
     * Adding -(value_of_know_node * element_of_full_stiffnes_matrix) 
     * to vector full_b
     * 
     * Later the rows and columns of full_stiffness_matrix 
     * corresponding to the known values in the nodes will be removed 
     * from full_stiffness_matrix (new reduced_matrix will be assembled)
    */

    for(std::size_t i = 0; i < full_b.size(); ++i){
        for(auto it = laplace_eq._ind_dirichlet_lower_boundary_nodes.begin(); it != laplace_eq._ind_dirichlet_lower_boundary_nodes.end(); ++it){
            solution[*it] = laplace_eq._dirichlet_lower_boundary_condition(laplace_eq._nodes[*it][0], laplace_eq._nodes[*it][1]);
            full_b[i] -= full_stiffness_matrix[i][*it] * solution[*it];
        }
    }
    for(std::size_t i = 0; i < full_b.size(); ++i){
        for(auto it = laplace_eq._ind_dirichlet_upper_boundary_nodes.begin(); it != laplace_eq._ind_dirichlet_upper_boundary_nodes.end(); ++it){
            solution[*it] = laplace_eq._dirichlet_upper_boundary_condition(laplace_eq._nodes[*it][0], laplace_eq._nodes[*it][1]);
            full_b[i] -= full_stiffness_matrix[i][*it] * solution[*it];
        }
    }
    for(std::size_t i = 0; i < full_b.size(); ++i){
        for(auto it = laplace_eq._ind_left_boundary_nodes.begin(); it != laplace_eq._ind_left_boundary_nodes.end(); ++it){
            solution[*it] = laplace_eq._dirichlet_left_boundary_condition(laplace_eq._nodes[*it][0], laplace_eq._nodes[*it][1]);
            full_b[i] -= full_stiffness_matrix[i][*it] * solution[*it];
        }
    }
    for(std::size_t i = 0; i < full_b.size(); ++i){
        for(auto it = laplace_eq._ind_right_boundary_nodes.begin(); it != laplace_eq._ind_right_boundary_nodes.end(); ++it){
            solution[*it] = laplace_eq._dirichlet_right_boundary_condition(laplace_eq._nodes[*it][0], laplace_eq._nodes[*it][1]);
            full_b[i] -= full_stiffness_matrix[i][*it] * solution[*it];
        }
    }

    /**
     * Applying Periodic boundary conditions
     * 
     * Summing values of full_stiffness matrix of nodes 
     * that correspond periodicity 
     * 
     * Later the rows and columns of full_stiffness_matrix 
     * corresponding to these nodes will be removed 
     * from full_stiffness_matrix (new reduced_matrix will be assembled)
    */

    // std::vector<std::pair<std::size_t, std::size_t>> ind_left_n_right_boundary_correlation;

    // for(std::size_t i = 0; i < laplace_eq._left_boundary_nodes.size(); ++i){
    //     for(std::size_t j = 0; j < laplace_eq._right_boundary_nodes.size(); ++j){
    //         if(laplace_eq._left_boundary_nodes[i][1] == laplace_eq._right_boundary_nodes[j][1]){
    //             ind_left_n_right_boundary_correlation.push_back(std::pair<std::size_t, std::size_t> (laplace_eq._ind_left_boundary_nodes[i], laplace_eq._ind_right_boundary_nodes[j]));
    //             // // For debugging purposes
    //             //std:: cout << laplace_eq._ind_left_boundary_nodes[i] << "\t" << laplace_eq._ind_right_boundary_nodes[j] << "\n";
    //         }
    //     }
    // }

    // for(std::size_t i = 0; i < full_stiffness_matrix.get_rows(); ++i){
    //     for(auto it = ind_left_n_right_boundary_correlation.begin(); it != ind_left_n_right_boundary_correlation.end(); ++it){
    //         full_stiffness_matrix[i][(*it).first] += full_stiffness_matrix[i][(*it).second];
    //     }
    // }


    /**
     * Making vector with indexes (in all nodes vector) of nodes 
     * which values should be calculated ( nodes that DON'T lie on upper (dirichlete condition), 
     * lower (dirichlete condition), and right (periodic condition) boundaries )
    */
    std::vector<std::size_t> nodes_ind_4_reduced_matrix;
    
    auto check_it_ind = std::find(laplace_eq._ind_dirichlet_lower_boundary_nodes.begin(), laplace_eq._ind_dirichlet_lower_boundary_nodes.end(), 0);
    
    for(std::size_t i = 0; i < laplace_eq._nodes.size(); ++i)   {
        check_it_ind = std::find(laplace_eq._ind_dirichlet_lower_boundary_nodes.begin(), laplace_eq._ind_dirichlet_lower_boundary_nodes.end(), i);
        
        if (check_it_ind == laplace_eq._ind_dirichlet_lower_boundary_nodes.end()){

            check_it_ind = std::find(laplace_eq._ind_dirichlet_upper_boundary_nodes.begin(), laplace_eq._ind_dirichlet_upper_boundary_nodes.end(), i);
            if (check_it_ind == laplace_eq._ind_dirichlet_upper_boundary_nodes.end()){

                check_it_ind = std::find(laplace_eq._ind_right_boundary_nodes.begin(), laplace_eq._ind_right_boundary_nodes.end(), i);
                if (check_it_ind == laplace_eq._ind_right_boundary_nodes.end()){

                    check_it_ind = std::find(laplace_eq._ind_left_boundary_nodes.begin(), laplace_eq._ind_left_boundary_nodes.end(), i);
                    if  (check_it_ind == laplace_eq._ind_left_boundary_nodes.end()){
                        nodes_ind_4_reduced_matrix.push_back(i);
                    }
                }
            }
        }
    }


    /**
     * Assembling reduced_stiffness_matrix
     * (matrix without unnecessary) nodes
    */

    std::size_t reduced_sys_size = 0;
    reduced_sys_size = laplace_eq._nodes.size() - laplace_eq._dirichlet_upper_boundary_nodes.size() - laplace_eq._dirichlet_lower_boundary_nodes.size() - laplace_eq._right_boundary_nodes.size() - laplace_eq._left_boundary_nodes.size();

    Matrix<T> reduced_matrix(reduced_sys_size, reduced_sys_size);
    Vector<T> reduced_b(reduced_sys_size);

    for(std::size_t i = 0; i < reduced_sys_size; ++i){
        for(std::size_t j = 0; j < reduced_sys_size; ++j){
            reduced_matrix[i][j] = full_stiffness_matrix[nodes_ind_4_reduced_matrix[i]][nodes_ind_4_reduced_matrix[j]];
        }
        reduced_b[i] = full_b[nodes_ind_4_reduced_matrix[i]];
    }
    
    
    // // For debugging purposes
    //std::cout << reduced_b << "\n";
    
    // output reduced stiffness matrix if needed
    if (output_stiffness_matricies_bool){
        fout.open(out_path + "_reduced_matrix.txt");
        for(std::size_t i = 0; i < reduced_matrix.get_rows(); ++i){
            for(std::size_t j = 0; j < reduced_matrix.get_cols(); ++j){            
                fout << reduced_matrix[i][j] << "\t";
            }
            fout << "\n";
        }
        fout << "\n";
        fout.close();  
    }


    /**
     * Assembling sparse matrix of Eigen::Sparce lib
     * 
     * And solving system
     * A.dot(x) = b
     * 
     * x = sparse_solution
     * A = sparse_matrix
     * b = RHS_vec
    */

    std::vector<Tr> tripletList(reduced_sys_size); // triplet list for filling sparse_matrix
    tripletList.reserve(reduced_sys_size);
    Eigen::VectorX<T> sparse_solution(reduced_sys_size), RHS_vec(reduced_sys_size);

    

    /**
     * Filling in triplets for sparse matrix of Eigen::Sparce lib
     * And RHS vector
    */    
    
    T tmp_el = 0.; //tmp_el to avoid multiple gettings of matrix & vector elements
    for(std::size_t i = 0; i < reduced_matrix.get_rows(); ++i){
        for(std::size_t j = 0; j < reduced_matrix.get_cols(); ++j){ 
            tmp_el = reduced_matrix[i][j];
            if (fabs(tmp_el) > tol){
                tripletList.push_back(Tr(i, j, tmp_el));
            }            
        }
        tmp_el = reduced_b[i];
        if (fabs(tmp_el) > tol)
            RHS_vec[i] = tmp_el;
        else RHS_vec[i] = 0.;
    }

    // Declare sparse matrix
    SpMat sparse_matrix(reduced_sys_size, reduced_sys_size);
    sparse_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    
    // Importing solver bi-conjugate
    Eigen::BiCGSTAB<SpMat> solver;
    //Eigen::LeastSquaresConjugateGradient<SpMat> solver;
    //Eigen::SparseLU<SpMat> solver;
    //Eigen::SimplicialLDLT<SpMat> solver;

    // Compute and solve (Eigen logic)
    solver.compute(sparse_matrix);
    sparse_solution = solver.solve(RHS_vec);
    
    /**
     * Assembling full solution (in each node)
     * 
    */
    //setting nodes from sparse solution
    for(std::size_t i = 0; i < reduced_sys_size; ++i){
        solution[nodes_ind_4_reduced_matrix[i]] = sparse_solution[i];
    }
    
    //setting nodes with periodic boundaries
    // for(auto it = ind_left_n_right_boundary_correlation.begin(); it != ind_left_n_right_boundary_correlation.end(); ++it){
    //     solution[(*it).second] = solution[(*it).first];
    // }
    // // For debugging purposes
    // for(std::size_t i = 0; i < laplace_eq._nodes.size(); ++i)
    //     std::cout << i << ": \t" << solution[i] << "\n";

    fout.open(out_path + "_solution.csv");
    fout << "sol\n";
    for(auto it = solution.begin(); it != solution.end(); ++it){
        fout << *it << "\n";
    }
    fout.close();
};