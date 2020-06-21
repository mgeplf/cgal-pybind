#include <string>
#include <vector>
#include <utility> // std::pair

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/extract_mean_curvature_flow_skeleton.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;

void bind_triangle_mesh(py::module& m)
{
    using Kernel = CGAL::Simple_cartesian<double>;
    using Point_2 = Kernel::Point_2;
    using Point_3 = Kernel::Point_3;

    using SurfaceMesh = CGAL::Surface_mesh<Point_3>;
    using Skeletonization = CGAL::Mean_curvature_flow_skeletonization<SurfaceMesh>;

    using IntIndices = std::vector<std::tuple<int, int, int>>;

    py::class_<SurfaceMesh>(m, "SurfaceMesh")
        .def(py::init())
        .def("add_vertices",
             [](SurfaceMesh& self, const std::vector<Point_3>& points) {
                 for(const auto& p : points) {
                     self.add_vertex(p);
                 }
             })
        .def("add_faces",
             [](SurfaceMesh& self, const IntIndices& indices) {
                for(auto [v0, v1, v2] : indices) {
                     self.add_face(static_cast<SurfaceMesh::Vertex_index>(v0),
                                   static_cast<SurfaceMesh::Vertex_index>(v1),
                                   static_cast<SurfaceMesh::Vertex_index>(v2));
                 }
             })
        .def("number_of_vertices", &SurfaceMesh::number_of_vertices)
        .def("number_of_faces", &SurfaceMesh::number_of_faces)
        .def("area", [](SurfaceMesh& self) {
             return CGAL::Polygon_mesh_processing::area(self);
        })
        .def("contract", [](SurfaceMesh& self) {
            using Skeletonization = CGAL::Mean_curvature_flow_skeletonization<SurfaceMesh>;

            Skeletonization mean_curve_skeletonizer(self);
            mean_curve_skeletonizer.contract_until_convergence();

            Skeletonization::Skeleton skeleton;
            mean_curve_skeletonizer.convert_to_skeleton(skeleton);

            std::vector<std::tuple<double, double, double>> return_vertices;
            for (const auto& v : CGAL::make_range(vertices(skeleton)))
            {
                return_vertices.emplace_back(skeleton[v].point[0],
                                             skeleton[v].point[1],
                                             skeleton[v].point[2]
                                             );
            }

            std::vector<std::pair<int, int>> return_edges;
            for (const auto& e : CGAL::make_range(edges(skeleton)))
            {
                return_edges.emplace_back(source(e, skeleton), target(e, skeleton));
            }

            return py::make_tuple(return_vertices, return_edges);
        })
        ;
}
