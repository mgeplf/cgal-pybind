#include <exception> // std::runtime_error
#include <utility> // std::pair
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

// { for contract
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
// }

// { for authalic
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
// }

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;

using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;

using SurfaceMesh = CGAL::Surface_mesh<Point_3>;

namespace {
std::vector<std::pair<double, double>>
run_discrete_authalic(SurfaceMesh &sm)
{
    namespace SMP = CGAL::Surface_mesh_parameterization;

    using BorderParameterizer = SMP::Square_border_arc_length_parameterizer_3<SurfaceMesh>;
    using Parameterizer = SMP::Discrete_authalic_parameterizer_3<SurfaceMesh, BorderParameterizer>;
    using halfedge_descriptor = boost::graph_traits<SurfaceMesh>::halfedge_descriptor;

    auto [bhd, length] = CGAL::Polygon_mesh_processing::longest_border(sm);

    // XXX: proper length test here? check bhd's count?
    if(length == 0.) {
        throw std::runtime_error("Longest border is 0.");
    }

    /*
    std::vector<int> return_vertices;
    for(auto vertex_index : CGAL::vertices_around_face(bhd, sm)) {
        return_vertices.push_back(vertex_index);
    }
    */

    //return return_vertices;

    using vertex_descriptor = boost::graph_traits<SurfaceMesh>::vertex_descriptor;
    using UV_pmap = SurfaceMesh::Property_map<vertex_descriptor, Point_2>;
    UV_pmap uv_map = sm.add_property_map<vertex_descriptor, Point_2>("v:uv").first;

    SMP::Error_code err = SMP::parameterize(sm, Parameterizer(), bhd, uv_map);
    if(err != SMP::OK) {
        std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
        throw std::runtime_error("Error parameterize-izing");
    }

    // https://github.com/CGAL/cgal/issues/2994
    using Vertex_index_map = boost::unordered_map<vertex_descriptor, std::size_t>;
    Vertex_index_map vium;
    boost::associative_property_map<Vertex_index_map> vimap(vium);

    std::vector<std::pair<double, double>> return_vertices;
    size_t vertices_counter = 0;
    for (const auto& v : CGAL::make_range(vertices(sm)))
    {
        auto& p = get(uv_map, v);
        return_vertices.emplace_back(p.x(), p.y());
        put(vimap, v, vertices_counter++);
    }
    return return_vertices;

    /*
    BOOST_FOREACH(face_descriptor fd, faces(sm)){
        halfedge_descriptor hd = halfedge(fd, sm);
        out << "3";
        BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd, sm)){
            out << " " << get(vimap, vd);
        }
        out << '\n';
        faces_counter++;
    }
    */
}
} // anonymous namespace


void bind_surface_mesh(py::module& m)
{
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
        .def("authalic", [](SurfaceMesh& self) {
            return py::make_tuple(run_discrete_authalic(self));

            /*
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
            */
        })
        ;
}
