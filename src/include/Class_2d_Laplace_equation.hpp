#pragma once

#include <vector>
#include <cmath>
#include <functional> 
#include <fstream>
#include <iostream>

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
        
        // std::vector<std::vector<T>> _dirichlet_left_boundary_nodes;
        // std::vector<std::vector<T>> _dirichlet_right_boundary_nodes;
        // std::vector<std::vector<T>> _dirichlet_upper_boundary_nodes;
        // std::vector<std::vector<T>> _dirichlet_lower_boundary_nodes;
        std::vector<std::vector<T>> _boundary_nodes;
        std::vector<std::vector<std::size_t>> _dirichlet_left_boundary_edges;
        std::vector<std::vector<std::size_t>> _dirichlet_right_boundary_edges;
        std::vector<std::vector<std::size_t>> _dirichlet_upper_boundary_edges;
        std::vector<std::vector<std::size_t>> _dirichlet_lower_boundary_edges;

        std::function<T (const T x, const T y)> _dirichlet_lower_boundary_condition;
        std::function<T (const T x, const T y)> _dirichlet_upper_boundary_condition;

        Class_2d_Laplace_equation<T> DEFAULT_TEST(){

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


            fin.open("/home/san/Code/Laplace_eq/meshes/mesh/mesh_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                this->_nodes.push_back(x_y);
            }            
            fin.close();

            //Reading mesh polygons (triangles)
            //2D array contains positions of elements in vector<T> _nodes (mesh_nodes.txt)
            //which form a triangulars

            fin.open("/home/san/Code/Laplace_eq/meshes/mesh/mesh_polygons.txt");
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

            fin.open("/home/san/Code/Laplace_eq/meshes/mesh/boundaries/mesh_boundary_nodes.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> x_y[0] >> x_y[1]) {
                this->_boundary_nodes.push_back(x_y);
            }            
            fin.close();


            //Reading boundaty edges
            //2D array contains positions of elements in vector<T> _nodes (mesh_nodes.txt)
            //which form a triangulars
            std::vector<std::size_t> p1_p2 = {0, 0};
            p1_p2.shrink_to_fit();
            
            //left boundary

            fin.open("/home/san/Code/Laplace_eq/meshes/mesh/boundaries/edges/mesh_left_boundary_edges.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> p1_p2[0] >> p1_p2[1]) {
                
                //minus 1 as Wolfram Mathematica numerates nodes from 1 to n
                //we need from 0 to n-1

                p1_p2[0] -= 1;
                p1_p2[1] -= 1;

                this->_dirichlet_left_boundary_edges.push_back(p1_p2);
                //std::cout << p1_p2_p3[0] << "\t" << p1_p2_p3[1] << "\t" << p1_p2_p3[2] << "\n";
            }            
            fin.close();


            //right boundary

            fin.open("/home/san/Code/Laplace_eq/meshes/mesh/boundaries/edges/mesh_right_boundary_edges.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> p1_p2[0] >> p1_p2[1]) {
                
                //minus 1 as Wolfram Mathematica numerates nodes from 1 to n
                //we need from 0 to n-1

                p1_p2[0] -= 1;
                p1_p2[1] -= 1;

                this->_dirichlet_left_boundary_edges.push_back(p1_p2);
                //std::cout << p1_p2_p3[0] << "\t" << p1_p2_p3[1] << "\t" << p1_p2_p3[2] << "\n";
            }            
            fin.close();

            
            //lower

            fin.open("/home/san/Code/Laplace_eq/meshes/mesh/boundaries/edges/mesh_lower_boundary_edges.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> p1_p2[0] >> p1_p2[1]) {
                
                //minus 1 as Wolfram Mathematica numerates nodes from 1 to n
                //we need from 0 to n-1

                p1_p2[0] -= 1;
                p1_p2[1] -= 1;

                this->_dirichlet_lower_boundary_edges.push_back(p1_p2);
                //std::cout << p1_p2_p3[0] << "\t" << p1_p2_p3[1] << "\t" << p1_p2_p3[2] << "\n";
            }            
            fin.close();


            //upper

            fin.open("/home/san/Code/Laplace_eq/meshes/mesh/boundaries/edges/mesh_upper_boundary_edges.txt");
            if (!fin.is_open()) {
                throw std::invalid_argument("Smth's wrong with path");
            }            
            
            while (fin >> p1_p2[0] >> p1_p2[1]) {
                
                //minus 1 as Wolfram Mathematica numerates nodes from 1 to n
                //we need from 0 to n-1

                p1_p2[0] -= 1;
                p1_p2[1] -= 1;

                this->_dirichlet_upper_boundary_edges.push_back(p1_p2);
                //std::cout << p1_p2_p3[0] << "\t" << p1_p2_p3[1] << "\t" << p1_p2_p3[2] << "\n";
            }            
            fin.close();


            return *this; 
        }

};