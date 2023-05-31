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
        std::function<T (const T x, const T y)> _dirichlet_left_boundary_condition;
        std::function<T (const T x, const T y)> _dirichlet_right_boundary_condition;

        Class_2d_Laplace_equation<T> Set_domain(Class_2d_Laplace_equation<T>& lap_eq, std::string path){

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

            // lap_eq._nodes.shrink_to_fit();
            // lap_eq._polygons.shrink_to_fit();

            // lap_eq._boundary_nodes.shrink_to_fit();
            // lap_eq._left_boundary_nodes.shrink_to_fit();
            // lap_eq._right_boundary_nodes.shrink_to_fit();
            // lap_eq._dirichlet_lower_boundary_nodes.shrink_to_fit();
            // lap_eq._dirichlet_upper_boundary_nodes.shrink_to_fit();

            // lap_eq._ind_left_boundary_nodes.shrink_to_fit();
            // lap_eq._ind_right_boundary_nodes.shrink_to_fit();
            // lap_eq._ind_dirichlet_lower_boundary_nodes.shrink_to_fit();
            // lap_eq._ind_dirichlet_upper_boundary_nodes.shrink_to_fit();
            
            //std::string path = "/home/san/Code/2d_Laplace_eq_FEM/meshes/mesh001_rect";
            // std::string path = "/home/san/Code/2d_Laplace_eq_FEM/meshes/mesh01";
            // std::string path = "/home/san/Code/2d_Laplace_eq_FEM/meshes/mesh005";
            // std::string path = "/home/san/Code/2d_Laplace_eq_FEM/meshes/mesh001";

            // lap_eq._dirichlet_lower_boundary_condition = [](const T x, const T y){
            //     return 10.;
            // };

            // lap_eq._dirichlet_upper_boundary_condition = [](const T x, const T y){
            //     return 0.;
            // };


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
            //check for 1 triangle element
            // std::cout << lap_eq._polygons[0][0] << "\t" << lap_eq._polygons[0][1] << "\t" << lap_eq._polygons[0][2] << "\n";
            // std::cout << lap_eq._nodes[lap_eq._polygons[0][0]][0] << "\t" << lap_eq._nodes[lap_eq._polygons[0][0]][1] << "\n";
            // std::cout << lap_eq._nodes[lap_eq._polygons[0][1]][0] << "\t" << lap_eq._nodes[lap_eq._polygons[0][1]][1] << "\n";
            // std::cout << lap_eq._nodes[lap_eq._polygons[0][2]][0] << "\t" << lap_eq._nodes[lap_eq._polygons[0][2]][1] << "\n";



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

            // for(std::size_t i = 0; i < lap_eq._right_boundary_nodes.size(); ++i){
            //     for(std::size_t j = 0; j < 2; ++j){
            //         std::cout << lap_eq._right_boundary_nodes[i][j] << "\t";
            //     }
            //     std::cout << "\n";
            // }


            fin.open(path + "/boundaries/mesh_lower_boundary_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                lap_eq._dirichlet_lower_boundary_nodes.push_back(x_y);
            }            
            fin.close();

            fin.open(path + "/boundaries/mesh_upper_boundary_nodes.txt");
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

};