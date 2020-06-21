#include <pybind11/pybind11.h>

#include <boost/lexical_cast.hpp> // boost::lexical_cast

#include <CGAL/Simple_cartesian.h>
using Point_3 = CGAL::Simple_cartesian<double>::Point_3;

namespace py = pybind11;

void bind_point(py::module& m)
{
    py::class_<Point_3>(m, "Point_3")
        .def(py::init<int, int, int>(), py::arg("x"), py::arg("y"),
             py::arg("z"))
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"),
             py::arg("z"))
        .def_property_readonly("x", &Point_3::x)
        .def_property_readonly("y", &Point_3::y)
        .def_property_readonly("z", &Point_3::z)
        .def("__repr__", [](const Point_3& self) {
            std::string r("Point_3(");
            r += boost::lexical_cast<std::string>(self.x());
            r += ", ";
            r += boost::lexical_cast<std::string>(self.y());
            r += ", ";
            r += boost::lexical_cast<std::string>(self.z());
            r += ")";
            return r;
        });
}
