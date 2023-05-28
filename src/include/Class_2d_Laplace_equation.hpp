#pragma once

#include <vector>
#include <cmath>
#include <functional> 
#include <fstream>
#include <iostream>
#include <algorithm>

/*
// ..........Class_2d_Laplace_equation............... //
*/
template<class T>
class Class_2d_Laplace_equation{

    public:

        std::string _name;

        //Space
        T _x0 = 0.;
        T _xL = 4.;

        std::vector<std::vector<T>> _nodes;
        std::vector<std::vector<std::size_t>> _polygons;


        std::vector<std::vector<T>> _boundary_nodes;
        std::vector<std::vector<T>> _left_boundary_nodes;
        std::vector<std::vector<T>> _right_boundary_nodes;
        std::vector<std::vector<T>> _dirichlet_lower_boundary_nodes;
        std::vector<std::vector<T>> _dirichlet_upper_boundary_nodes;

        //indexes of boundary nodes in _nodes vector
        std::vector<std::size_t> _ind_left_boundary_nodes;
        std::vector<std::size_t> _ind_right_boundary_nodes;
        std::vector<std::size_t> _ind_dirichlet_lower_boundary_nodes;
        std::vector<std::size_t> _ind_dirichlet_upper_boundary_nodes;

        std::function<T (const T x, const T y)> _dirichlet_lower_boundary_condition;
        std::function<T (const T x, const T y)> _dirichlet_upper_boundary_condition;

        Class_2d_Laplace_equation<T> DEFAULT_TEST(){

            std::string path = "/home/san/Code/2d_Laplace_eq_FEM/meshes/mesh01_rect";

            this->_dirichlet_lower_boundary_condition = [](const T x, const T y){
                return 10.;
            };

            this->_dirichlet_upper_boundary_condition = [](const T x, const T y){
                return 0.;
            };



            std::ifstream fin;

            //Reading mesh nodes (coordinates)
            std::vector<T> x_y = {0., 0.};
            x_y.shrink_to_fit();


            fin.open(path + "/mesh_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                this->_nodes.push_back(x_y);
            }            
            fin.close();

            //std::cout << this->_nodes.size();

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

                this->_polygons.push_back(p1_p2_p3);
                //std::cout << p1_p2_p3[0] << "\t" << p1_p2_p3[1] << "\t" << p1_p2_p3[2] << "\n";
            }            
            fin.close();
            //check for 1 triangle element
            // std::cout << this->_polygons[0][0] << "\t" << this->_polygons[0][1] << "\t" << this->_polygons[0][2] << "\n";
            // std::cout << this->_nodes[this->_polygons[0][0]][0] << "\t" << this->_nodes[this->_polygons[0][0]][1] << "\n";
            // std::cout << this->_nodes[this->_polygons[0][1]][0] << "\t" << this->_nodes[this->_polygons[0][1]][1] << "\n";
            // std::cout << this->_nodes[this->_polygons[0][2]][0] << "\t" << this->_nodes[this->_polygons[0][2]][1] << "\n";



            //Reading boundary nodes (coordinates)

            fin.open(path + "/boundaries/mesh_boundary_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                this->_boundary_nodes.push_back(x_y);
            }            
            fin.close();


            fin.open(path + "/boundaries/mesh_left_boundary_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                this->_left_boundary_nodes.push_back(x_y);
            }            
            fin.close();


            fin.open(path + "/boundaries/mesh_right_boundary_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                this->_right_boundary_nodes.push_back(x_y);
            }            
            fin.close();

            // for(std::size_t i = 0; i < this->_right_boundary_nodes.size(); ++i){
            //     for(std::size_t j = 0; j < 2; ++j){
            //         std::cout << this->_right_boundary_nodes[i][j] << "\t";
            //     }
            //     std::cout << "\n";
            // }


            fin.open(path + "/boundaries/mesh_lower_boundary_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                this->_dirichlet_lower_boundary_nodes.push_back(x_y);
            }            
            fin.close();

            fin.open(path + "/boundaries/mesh_upper_boundary_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                this->_dirichlet_upper_boundary_nodes.push_back(x_y);
            }            
            fin.close();

            //find indexes of boundary nodes in _nodes vector
            auto it = std::find(this->_nodes.begin(), this->_nodes.end(), this->_left_boundary_nodes[0]);

            for(std::size_t i = 0; i < this->_left_boundary_nodes.size(); ++i){
                it = std::find(this->_nodes.begin(), this->_nodes.end(), this->_left_boundary_nodes[i]);
                if (it  != this->_nodes.end()){
                    this->_ind_left_boundary_nodes.push_back(std::distance(this->_nodes.begin(), it));
                }
            }
            // for debuging
            // for(std::size_t i = 0; i < this->_left_boundary_nodes.size(); ++i){
            //     std::cout << this->_ind_left_boundary_nodes[i] << "\t";
            // }
            // std::cout << "\n";


            for(std::size_t i = 0; i < this->_right_boundary_nodes.size(); ++i){
                it = std::find(this->_nodes.begin(), this->_nodes.end(), this->_right_boundary_nodes[i]);
                if (it  != this->_nodes.end()){
                    this->_ind_right_boundary_nodes.push_back(std::distance(this->_nodes.begin(), it));
                }
            }
            
            // for debuging
            // for(std::size_t i = 0; i < this->_right_boundary_nodes.size(); ++i){
            //     std::cout << this->_ind_right_boundary_nodes[i] << "\t";
            // }
            // std::cout << "\n";

            for(std::size_t i = 0; i < this->_dirichlet_lower_boundary_nodes.size(); ++i){
                it = std::find(this->_nodes.begin(), this->_nodes.end(), this->_dirichlet_lower_boundary_nodes[i]);
                if (it  != this->_nodes.end()){
                    this->_ind_dirichlet_lower_boundary_nodes.push_back(std::distance(this->_nodes.begin(), it));
                }
            }

            // // for debuging
            // for(std::size_t i = 0; i < this->_dirichlet_lower_boundary_nodes.size(); ++i){
            //     std::cout << this->_ind_dirichlet_lower_boundary_nodes[i] << "\t";
            // }
            // std::cout << "\n";

            for(std::size_t i = 0; i < this->_dirichlet_upper_boundary_nodes.size(); ++i){
                it = std::find(this->_nodes.begin(), this->_nodes.end(), this->_dirichlet_upper_boundary_nodes[i]);
                if (it  != this->_nodes.end()){
                    this->_ind_dirichlet_upper_boundary_nodes.push_back(std::distance(this->_nodes.begin(), it));
                }
            }

            return *this; 
        }

};