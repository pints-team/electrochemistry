#ifndef E_IMPLICIT_EXPONENTIAL_MESH_HPP
#define E_IMPLICIT_EXPONENTIAL_MESH_HPP

#include "utilities.hpp"

#include <map>
#include <string>
#include <vector>

namespace electrochemistry {
void e_implicit_exponential_mesh(py::dict params, py::array_t<double> Itot,
                                 py::array_t<double> t);
}

#endif
