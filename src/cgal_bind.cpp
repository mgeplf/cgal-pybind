#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_point(py::module&);
void bind_triangle_mesh(py::module&);
void bind_facet_with_id_pmap(py::module&);
void bind_polyhedron_inc_builder(py::module&);

PYBIND11_MODULE(cgal_pybind, m)
{
    m.doc() = "python binding for CGAL (https://doc.cgal.org/latest/Manual/index.html)";
    bind_point(m);
    bind_triangle_mesh(m);
}
