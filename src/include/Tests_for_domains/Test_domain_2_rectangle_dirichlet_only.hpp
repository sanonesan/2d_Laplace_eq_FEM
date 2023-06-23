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


            auto exact_solution = [](T x, T y, std::size_t row_num_el = 350) -> T {
                T res = 0.;
                for (std::size_t n = 1; n < row_num_el + 1; ++n) {
                    if ((n % 2) != 0 ) {
                        res += 20 * (1 - cos(M_PI * n)) / M_PI / n / (exp(- M_PI * n / 2) - exp( M_PI * n / 2)) 
                                * (exp(- M_PI * n / 4 * y) - exp( M_PI * n / 4 * y)) * sin(M_PI * n / 4 * x);
                        }
                    }                        
                if (fabs(res) < 1e-15) return 0;
                                
                return res;
                
            };

            
            //test name
            this->_name = std::string (__func__);

            this->_dirichlet_lower_boundary_condition = [exact_solution](const T x, const T y){
                return exact_solution(x, y);
            };

            this->_dirichlet_upper_boundary_condition = [exact_solution](const T x, const T y){
                return exact_solution(x, y);
            };

            this->_dirichlet_left_boundary_condition = [exact_solution](const T x, const T y){
                return exact_solution(x, y);
            };

            this->_dirichlet_right_boundary_condition = [exact_solution](const T x, const T y){
                return exact_solution(x, y);
            };



            // Other boundary conditions

            // this->_dirichlet_lower_boundary_condition = [exact_solution](const T x, const T y){
            //     return 0.;
            // };

            // this->_dirichlet_upper_boundary_condition = [exact_solution](const T x, const T y){
            //     return 10.;
            // };

            // this->_dirichlet_left_boundary_condition = [exact_solution](const T x, const T y){
            //     return 0.;
            // };

            // this->_dirichlet_right_boundary_condition = [exact_solution](const T x, const T y){
            //     return 0.;
            // };

            this->Set_domain(*this, "../domains/domain_2/mesh001");
        }

        
        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_005(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_005");
            this->Set_domain(*this, "../domains/domain_2/mesh005");
            return *this;
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_005_rev_calfem(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_005_rev_calfem");

            auto exact_solution_rev = [](T x, T y, std::size_t row_num_el = 350) -> T {
                T res = 0.;
                for (std::size_t n = 1; n < row_num_el + 1; ++n) {
                    if ((n % 2) != 0 ) {
                        res += 20 * (1 - cos(M_PI * n)) / M_PI / n / (exp(- M_PI * n / 2) - exp( M_PI * n / 2)) 
                                * (exp(- M_PI * n / 4 * (2 - y)) - exp( M_PI * n / 4 * (2 - y))) * sin(M_PI * n / 4 * x);
                        }
                    }                        
                if (fabs(res) < 1e-15) return 0;
                                
                return res;
                
            };

            this->_dirichlet_lower_boundary_condition = [exact_solution_rev](const T x, const T y){
                return exact_solution_rev(x, y);
            };

            this->_dirichlet_upper_boundary_condition = [exact_solution_rev](const T x, const T y){
                return exact_solution_rev(x, y);
            };

            this->_dirichlet_left_boundary_condition = [exact_solution_rev](const T x, const T y){
                return exact_solution_rev(x, y);
            };

            this->_dirichlet_right_boundary_condition = [exact_solution_rev](const T x, const T y){
                return exact_solution_rev(x, y);
            };


            this->Set_domain(*this, "../domains/domain_2/mesh005_calfem");
            return *this;
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_005_calfem(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_005_calfem");
            this->Set_domain(*this, "../domains/domain_2/mesh005_calfem");
            return *this;
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_001(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_001");
            this->Set_domain(*this, "../domains/domain_2/mesh001");
            return *this;
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_001_calfem(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_001_calfem");
            this->Set_domain(*this, "../domains/domain_2/mesh001_calfem");
            return *this;
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_0005_calfem(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_0005_calfem");
            this->Set_domain(*this, "../domains/domain_2/mesh0005_calfem");
            return *this;
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_0001(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_0001");
            this->Set_domain(*this, "../domains/domain_2/mesh0001");
            return *this;
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_0001_calfem(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_0001_calfem");
            this->Set_domain(*this, "../domains/domain_2/mesh0001_calfem");
            return *this;
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_00005_calfem(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_00005_calfem");
            this->Set_domain(*this, "../domains/domain_2/mesh00005_calfem");
            return *this;
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_00008(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_00008");
            this->Set_domain(*this, "../domains/domain_2/mesh00008");
            return *this;
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_00008_calfem(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_00008_calfem");
            this->Set_domain(*this, "../domains/domain_2/mesh00008_calfem");
            return *this;
        }

        Test_domain_2_rectangle_dirichlet_only<T> Set_mesh_00008_s(){

            this->_name = std::string ("Test_domain_2_rectangle_dirichlet_only_00008_s");
            this->Set_domain(*this, "../domains/domain_2/mesh00008_s");
            return *this;
        }


};