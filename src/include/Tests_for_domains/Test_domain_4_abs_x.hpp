#pragma once

#include <cmath>
#include <functional>
#include <vector>

#include "../Class_2d_Laplace_equation.hpp"

/**
 * Test in domain:
 * {(x, y) | w(x) <= y <= 2, x \in [0, 2]}
 * w(x) = abs(x - 1) + 1
 * V(x, 2) = 0.
 * V(x, w(x)) = 10.
 * V(0., y) = V(2, y)
 *
 */
template <class T>
class Test_domain_4_abs_x : virtual public Class_2d_Laplace_equation<T> {
    public:
        Test_domain_4_abs_x() {
            
            this->_name = std::string("Test_domain_4_mesh001_3_in_row_calfem");
            
            this->_dirichlet_lower_boundary_condition = [](const T x, const T y) {
                return 10.;
            };

            this->_dirichlet_upper_boundary_condition = [](const T x, const T y) {
                return 0.;
            };

            this->Set_domain(*this, "../domains/domain_4/mesh001_3_in_row_calfem");
        }


        Test_domain_4_abs_x<T> Set_mesh001_calfem() {

            this->_name = std::string("Test_domain_4_mesh001_calfem");
            this->Set_domain(*this, "../domains/domain_4/mesh001_calfem");
            return *this;
        }


        Test_domain_4_abs_x<T> Set_mesh001_3_in_row_calfem() {

            this->_name = std::string("Test_domain_4_mesh001_3_in_row_calfem");
            this->Set_domain(*this, "../domains/domain_4/mesh001_3_in_row_calfem");
            return *this;
        }
};