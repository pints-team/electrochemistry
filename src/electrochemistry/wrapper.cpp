#include <pybind11/pybind11.h>
#include "e_implicit_exponential_mesh.hpp"
#include "seq_electron_transfer3_explicit.hpp"

namespace py = pybind11;
using namespace electrochemistry;

PYBIND11_MODULE(electrochemistry, m) {
    m.def("e_implicit_exponential_mesh", &e_implicit_exponential_mesh);
    m.def("seq_electron_transfer1_explicit",
          &seq_electron_transfer_explicit<1>);
    m.def("seq_electron_transfer2_explicit",
          &seq_electron_transfer_explicit<2>);
    m.def("seq_electron_transfer3_explicit",
          &seq_electron_transfer_explicit<3>);
}
