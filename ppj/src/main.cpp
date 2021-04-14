#include <pybind11/pybind11.h>

double axpb(int a, float x, double b) {
    return a * x + b;
}

namespace py = pybind11;

PYBIND11_MODULE(my_module, m) {
    m.doc() = "docstring";

    m.def("axpb", &axpb, "calculate a*x + b");

    m.attr("__version__") = "dev";
}
