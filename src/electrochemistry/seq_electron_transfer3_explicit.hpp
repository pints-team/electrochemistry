#ifndef SEQ_ELECTRON_TRANSFER3_EXPLICIT_HPP
#define SEQ_ELECTRON_TRANSFER3_EXPLICIT_HPP

#include "utilities.hpp"

#include <map>
#include <string>
#include <vector>

namespace electrochemistry {
void seq_electron_transfer3_explicit(py::dict params,
                                     py::array_t<double> Itot_numpy,
                                     py::array_t<double> t_numpy);
}

#endif
