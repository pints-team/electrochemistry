#ifndef SEQ_ELECTRON_TRANSFER3_EXPLICIT_HPP
#define SEQ_ELECTRON_TRANSFER3_EXPLICIT_HPP

#include <map>
#include <string>
#include <vector>

#include "utilities.hpp"

namespace electrochemistry {
void seq_electron_transfer3_explicit(py::dict params,
                                     py::array_t<double> Itot_numpy,
                                     py::array_t<double> t_numpy);
}

#endif
