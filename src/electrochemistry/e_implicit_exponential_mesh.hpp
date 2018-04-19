#ifndef E_IMPLICIT_EXPONENTIAL_MESH_HPP
#define E_IMPLICIT_EXPONENTIAL_MESH_HPP

#include <map>
#include <string>
#include <vector>

#include "utilities.hpp"

namespace pints {
void e_implicit_exponential_mesh(py::dict params, py::array Itot, py::array t);
}

#endif
