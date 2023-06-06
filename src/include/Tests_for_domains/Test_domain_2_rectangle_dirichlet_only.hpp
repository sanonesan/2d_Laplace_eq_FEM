#pragma once

#include<vector>
#include<cmath>
#include<functional> 

#include "../Class_2d_Laplace_equation.hpp"


/**
 * Test in domain: 
 * {(x, y) | 0 <= y <= 2, x \in [0, 4]}
 * 
 * V(x, 2) = 0.
 * V(x, 0) = 10.
 * V(0., y) = V(4., y)
 * 
*/
template<class T>
class Test_domain_2_rectangle_dirichlet_only: virtual public Class_2d_Laplace_equation<T>{

    public:

        Test_domain_2_rectangle_dirichlet_only(){            
            //test name
            this->_name = std::string (__func__);

            this->_dirichlet_lower_boundary_condition = [](const T x, const T y){
                return 0.;
            };

            this->_dirichlet_upper_boundary_condition = [](const T x, const T y){
                return 10.;
            };

            this->_dirichlet_left_boundary_condition = [](const T x, const T y){
                return 0.;
            };

            this->_dirichlet_right_boundary_condition = [](const T x, const T y){
                return 0.;
            };

            this->Set_domain(*this, "../domains/domain_2/mesh001");
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_domain_dir(Class_2d_Laplace_equation<T>& lap_eq, std::string path){

            lap_eq._nodes.clear();
            lap_eq._polygons.clear();

            lap_eq._boundary_nodes.clear();
            lap_eq._left_boundary_nodes.clear();
            lap_eq._right_boundary_nodes.clear();
            lap_eq._dirichlet_lower_boundary_nodes.clear();
            lap_eq._dirichlet_upper_boundary_nodes.clear();

            lap_eq._ind_left_boundary_nodes.clear();
            lap_eq._ind_right_boundary_nodes.clear();
            lap_eq._ind_dirichlet_lower_boundary_nodes.clear();
            lap_eq._ind_dirichlet_upper_boundary_nodes.clear();


            std::ifstream fin;

            //Reading mesh nodes (coordinates)
            std::vector<T> x_y = {0., 0.};
            x_y.shrink_to_fit();


            fin.open(path + "/mesh_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                lap_eq._nodes.push_back(x_y);
            }            
            fin.close();

            //std::cout << lap_eq._nodes.size();

            //Reading mesh polygons (triangles)
            //2D array contains positions of elements in vector<T> _nodes (mesh_nodes.txt)
            //which form a triangulars

            fin.open(path + "/mesh_polygons.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            

            std::vector<std::size_t> p1_p2_p3 = {0, 0, 0};
            p1_p2_p3.shrink_to_fit();
            
            while (fin >> p1_p2_p3[0] >> p1_p2_p3[1] >> p1_p2_p3[2]) {
                
                //minus 1 as Wolfram Mathematica numerates nodes from 1 to n
                //we need from 0 to n-1

                p1_p2_p3[0] -= 1;
                p1_p2_p3[1] -= 1;
                p1_p2_p3[2] -= 1;

                lap_eq._polygons.push_back(p1_p2_p3);
                //std::cout << p1_p2_p3[0] << "\t" << p1_p2_p3[1] << "\t" << p1_p2_p3[2] << "\n";
            }            
            fin.close();
            
            //Reading boundary nodes (coordinates)

            fin.open(path + "/boundaries/mesh_boundary_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                lap_eq._boundary_nodes.push_back(x_y);
            }            
            fin.close();


            fin.open(path + "/boundaries/mesh_left_boundary_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                lap_eq._left_boundary_nodes.push_back(x_y);
            }            
            fin.close();


            fin.open(path + "/boundaries/mesh_right_boundary_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                lap_eq._right_boundary_nodes.push_back(x_y);
            }            
            fin.close();


            fin.open(path + "/boundaries/mesh_lower_boundary_nodes_dir.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                lap_eq._dirichlet_lower_boundary_nodes.push_back(x_y);
            }            
            fin.close();

            fin.open(path + "/boundaries/mesh_upper_boundary_nodes_dir.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                lap_eq._dirichlet_upper_boundary_nodes.push_back(x_y);
            }            
            fin.close();

            //find indexes of boundary nodes in _nodes vector
            auto it = std::find(lap_eq._nodes.begin(), lap_eq._nodes.end(), lap_eq._left_boundary_nodes[0]);

            for(std::size_t i = 0; i < lap_eq._left_boundary_nodes.size(); ++i){
                it = std::find(lap_eq._nodes.begin(), lap_eq._nodes.end(), lap_eq._left_boundary_nodes[i]);
                if (it  != lap_eq._nodes.end()){
                    lap_eq._ind_left_boundary_nodes.push_back(std::distance(lap_eq._nodes.begin(), it));
                }
            }
            // for debuging
            // for(std::size_t i = 0; i < lap_eq._left_boundary_nodes.size(); ++i){
            //     std::cout << lap_eq._ind_left_boundary_nodes[i] << "\t";
            // }
            // std::cout << "\n";


            for(std::size_t i = 0; i < lap_eq._right_boundary_nodes.size(); ++i){
                it = std::find(lap_eq._nodes.begin(), lap_eq._nodes.end(), lap_eq._right_boundary_nodes[i]);
                if (it  != lap_eq._nodes.end()){
                    lap_eq._ind_right_boundary_nodes.push_back(std::distance(lap_eq._nodes.begin(), it));
                }
            }
            
            // for debuging
            // for(std::size_t i = 0; i < lap_eq._right_boundary_nodes.size(); ++i){
            //     std::cout << lap_eq._ind_right_boundary_nodes[i] << "\t";
            // }
            // std::cout << "\n";

            for(std::size_t i = 0; i < lap_eq._dirichlet_lower_boundary_nodes.size(); ++i){
                it = std::find(lap_eq._nodes.begin(), lap_eq._nodes.end(), lap_eq._dirichlet_lower_boundary_nodes[i]);
                if (it  != lap_eq._nodes.end()){
                    lap_eq._ind_dirichlet_lower_boundary_nodes.push_back(std::distance(lap_eq._nodes.begin(), it));
                }
            }

            // // for debuging
            // for(std::size_t i = 0; i < lap_eq._dirichlet_lower_boundary_nodes.size(); ++i){
            //     std::cout << lap_eq._ind_dirichlet_lower_boundary_nodes[i] << "\t";
            // }
            // std::cout << "\n";

            for(std::size_t i = 0; i < lap_eq._dirichlet_upper_boundary_nodes.size(); ++i){
                it = std::find(lap_eq._nodes.begin(), lap_eq._nodes.end(), lap_eq._dirichlet_upper_boundary_nodes[i]);
                if (it  != lap_eq._nodes.end()){
                    lap_eq._ind_dirichlet_upper_boundary_nodes.push_back(std::distance(lap_eq._nodes.begin(), it));
                }
            }

            return *this; 
        }


        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_005(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_005");
            this->Set_domain_dir(*this, "../domains/domain_2/mesh005");
            return *this;
        }


        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_001(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_001");
            this->Set_domain_dir(*this, "../domains/domain_2/mesh001");
            return *this;
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_0001(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_0001");
            this->Set_domain_dir(*this, "../domains/domain_2/mesh0001");
            return *this;
        }


};