#pragma once

#include<vector>
#include<cmath>
#include<functional> 

#include "../Class_2d_Laplace_equation.hpp"


/**
 * Test in domain: 
 * {(x, y) | sin( pi / 2 * x) <= y <= 2, x \in [0, 4]}
 * 
 * V(x, 2) = 0.
 * V(x, sin( pi / 2 * x )) = 10.
 * V(0., y) = V(4., y)
 * 
*/
template<class T>
class Test_domain_1_1_sin: virtual public Class_2d_Laplace_equation<T>{

    public:

        Test_domain_1_1_sin(){            
            //test name
            //this->_name = std::string (__func__);
            this->_name = std::string ("Test_domain_1_1_sin_mesh_001");


            this->_dirichlet_lower_boundary_condition = [](const T x, const T y){
                return 10.;
            };

            this->_dirichlet_upper_boundary_condition = [](const T x, const T y){
                return 0.;
            };

            this->Set_domain(*this, "../domains/domain_1_1/mesh001");
        }


        // Test_domain_1_1_sin<T> Set_mesh_01(){

        //     this->_name = std::string ("Test_domain_1_1_sin_mesh_01");

        //     this->_dirichlet_lower_boundary_condition = [](const T x, const T y){
        //         return 10.;
        //     };

        //     this->_dirichlet_upper_boundary_condition = [](const T x, const T y){
        //         return 0.;
        //     };

        //     this->Set_domain(*this, "../domains/domain_1_1/mesh01");
        //     return *this;
        // }


        


        Test_domain_1_1_sin<T> Set_mesh_001(){

            this->_name = std::string ("Test_domain_1_1_sin_mesh_001");

            this->_dirichlet_lower_boundary_condition = [](const T x, const T y){
                return 10.;
            };

            this->_dirichlet_upper_boundary_condition = [](const T x, const T y){
                return 0.;
            };

            this->Set_domain(*this, "../domains/domain_1_1/mesh001");
            return *this;
        }


        Test_domain_1_1_sin<T> Set_mesh_0005_calfem(){

            this->_name = std::string ("Test_domain_1_1_sin_mesh_0005_calfem");

            this->_dirichlet_lower_boundary_condition = [](const T x, const T y){
                return -7.;
            };
            this->_dirichlet_upper_boundary_condition = [](const T x, const T y){
                return 12.;
            };

            this->Set_domain(*this, "../domains/domain_1_1/mesh0005_calfem");
            return *this;
        }


        Test_domain_1_1_sin<T> Set_mesh_0005_3_in_row_calfem(){

            this->_name = std::string ("Test_domain_1_1_sin_mesh_0005_3_in_row_calfem");

            this->_dirichlet_lower_boundary_condition = [](const T x, const T y){
                return -7.;
            };
            this->_dirichlet_upper_boundary_condition = [](const T x, const T y){
                return 12.;
            };

            this->Set_domain(*this, "../domains/domain_1_1/mesh0005_3_in_row_calfem");
            return *this;
        }


        Test_domain_1_1_sin<T> Set_mesh_0001_calfem(){

            this->_name = std::string ("Test_domain_1_1_sin_mesh_0001_calfem");

            this->_dirichlet_lower_boundary_condition = [](const T x, const T y){
                return -7.;
            };
            this->_dirichlet_upper_boundary_condition = [](const T x, const T y){
                return 12.;
            };

            this->Set_domain(*this, "../domains/domain_1_1/mesh0001_calfem");
            return *this;
        }


        Test_domain_1_1_sin<T> Set_mesh_0001_3_in_row_calfem(){

            this->_name = std::string ("Test_domain_1_1_sin_mesh_0001_3_in_row_calfem");

            this->_dirichlet_lower_boundary_condition = [](const T x, const T y){
                return -7.;
            };
            this->_dirichlet_upper_boundary_condition = [](const T x, const T y){
                return 12.;
            };

            this->Set_domain(*this, "../domains/domain_1_1/mesh0001_3_in_row_calfem");
            return *this;
        }

};