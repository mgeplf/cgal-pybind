#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_point(py::module&);
void bind_triangle_mesh(py::module&);

PYBIND11_MODULE(cgal_pybind, m)
{
    bind_point(m);
    bind_triangle_mesh(m);
}
