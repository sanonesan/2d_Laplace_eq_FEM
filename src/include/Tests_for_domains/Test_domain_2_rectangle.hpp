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
class Test_domain_2_rectangle: virtual public Class_2d_Laplace_equation<T>{

    public:

        Test_domain_2_rectangle(){            
            //test name
            this->_name = std::string (__func__);

            this->_dirichlet_lower_boundary_condition = [](const T x, const T y){
                return 10.;
            };

            this->_dirichlet_upper_boundary_condition = [](const T x, const T y){
                return 0.;
            };

            this->Set_domain(*this, "../domains/domain_2/mesh001");
        }


        Test_domain_2_rectangle<T> Set_mesh_001(){

            this->_name = std::string ("Test_domain_1_sin_mesh_001");
            this->Set_domain(*this, "../domains/domain_1/mesh001");
            return *this;
        }


};