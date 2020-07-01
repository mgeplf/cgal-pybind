#include <pybind11/pybind11.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>
#include <iostream>


namespace py = pybind11;

py::array return_numpy_array(){
    std::vector<double> cpp_vec;
    cpp_vec.push_back(1.);
    cpp_vec.push_back(2.);
    for (auto val:cpp_vec) {
        std::cout << val << std::endl;
    }
    ssize_t  ndim    = 2;


    std::vector<ssize_t> shape   = { 2, 1 };
    std::vector<ssize_t> strides = { sizeof(double)*3 , sizeof(double) };

    // return 2-D NumPy array
    arr1 = py::array(py::buffer_info(
    cpp_vec.data(),                           /* data as contiguous array  */
    sizeof(double),                          /* size of one scalar        */
    py::format_descriptor<double>::format(), /* data type                 */
    ndim,                                    /* number of dimensions      */
    shape,                                   /* shape of the matrix       */
    strides                                  /* strides for each axis     */
  ));
  arr2 = py::array(py::buffer_info(
    cpp_vec.data(),                           /* data as contiguous array  */
    sizeof(double),                          /* size of one scalar        */
    py::format_descriptor<double>::format(), /* data type                 */
    ndim,                                    /* number of dimensions      */
    shape,                                   /* shape of the matrix       */
    strides                                  /* strides for each axis     */
  ));
  rerturn arr1

}

void bind_test(py::module& m)
{
    py::class_<Test>(m, "Test")
        .def(py::init())
        })
     .def("return_numpy_array", &return_numpy_array)
     ;

}
