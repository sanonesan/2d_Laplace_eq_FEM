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
class Test_domain_3_sin_sin: virtual public Class_2d_Laplace_equation<T>{

    public:

        Test_domain_3_sin_sin(){            
            //test name
            //this->_name = std::string (__func__);
            this->_name = std::string ("Test_domain_3_sin_sin_mesh_001");


            this->_dirichlet_lower_boundary_condition = [](const T x, const T y){
                return 10.;
            };

            this->_dirichlet_upper_boundary_condition = [](const T x, const T y){
                return 0.;
            };

            this->Set_domain(*this, "../domains/domain_3/mesh001");
        }


};