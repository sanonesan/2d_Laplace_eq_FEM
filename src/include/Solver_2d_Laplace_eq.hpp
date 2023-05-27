#pragma once

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <vector>
#include "Class_2d_Laplace_equation.hpp"


template<class T>
class Solver_2d_Laplace_eq{

    private:

        void check_folder(const std::string& str){

            std::ifstream dir_stream(str.c_str());

            if (!dir_stream) {
                std::cout << "Folder created!\n" << "path:\t" << str << "\n\n";
                const int dir_err = mkdir(str.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                if (dir_err == -1)
                {
                    std::cout << ("Error creating directory!\n");
                    exit(1);
                }               
                
            }

        }

    public:

        T tol = 1e-6;
        bool notifications = false;
        std::string file_name = "";
        std::string output_folder = "../output/";

        //Default constructor
        Solver_2d_Laplace_eq(){
            this->tol = 1e-9;
            this->output_folder = "../output/";
            this->file_name = "";
        }

        //Alt constructor
        Solver_2d_Laplace_eq(T tol, std::string output_folder, std::string file_name){
            this->tol = tol;
            this->output_folder = "../output/";
            this->file_name = file_name;
        }

        void solve_eq(Class_2d_Laplace_equation<T>& Laplace_equation){

            this->check_folder(this->output_folder);
            std::string out_path;

            out_path = this->output_folder + this->file_name;

            if (this->notifications){
                std::cout << file_name << ": \t";
            }

            out_path += "_2d_Laplace_eq_output";                

        }


};