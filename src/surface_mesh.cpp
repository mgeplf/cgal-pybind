#include <string>
#include <vector>
#include <utility> // std::pair

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/mesh_segmentation.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


// Property map associating a facet with an integer as id to an
// element in a vector stored internally
template<class ValueType>
struct Facet_with_id_pmap
        : public boost::put_get_helper<ValueType&, Facet_with_id_pmap<ValueType> >
{
    using key_type =  CGAL::Polyhedron_3<
        CGAL::Exact_predicates_inexact_constructions_kernel,
        CGAL::Polyhedron_items_with_id_3
    >::Facet_const_handle;
    using category = boost::lvalue_property_map_tag;

    Facet_with_id_pmap(std::vector<ValueType>& internal_vector) : internal_vector(internal_vector) { }

    ValueType& operator[](key_type key) const
    { return internal_vector[key->id()]; }

private:
    std::vector<ValueType>& internal_vector;
};


namespace py = pybind11;

void bind_triangle_mesh(py::module& m)
{
/*
    using Kernel = CGAL::Simple_cartesian<double>;

    using Point_2 = Kernel::Point_2;
    using Point_3 = Kernel::Point_3;
    using Polyhedron = CGAL::Polyhedron_3<Kernel>;

    using SurfaceMesh = CGAL::Surface_mesh<Point_3>;
    using Skeletonization = CGAL::Mean_curvature_flow_skeletonization<SurfaceMesh>;
*/
   using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
   using Point_3 = Kernel::Point_3;
   using Polyhedron = CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>;
   using IntIndices = std::vector<std::tuple<int, int, int>>;

    py::class_<Polyhedron>(m, "Polyhedron")
        .def(py::init())
        .def("load_from_off",[](Polyhedron& self, const std::string& full_path_name) {
            std::ifstream input(full_path_name);
            input >> self;
        })
        .def("contract", [](Polyhedron& self) {
            using Skeletonization = CGAL::Mean_curvature_flow_skeletonization<Polyhedron> ;
            using Skeleton  = Skeletonization::Skeleton;
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
        .def("segmentation", [](Polyhedron& self){

            // assign a unique id to each facet in the mesh
            std::size_t facet_id = 0;
            for(auto facet = self.facets_begin(); facet != self.facets_end(); ++facet) {
                facet->id() = facet_id;
                ++facet_id;
            }

            // create a property-map for signed distance function values
            std::vector<double> sdf_values(self.size_of_facets());
            Facet_with_id_pmap<double> sdf_property_map(sdf_values);

            std::pair<double, double> min_max_sdf = CGAL::sdf_values(self, sdf_property_map);
            // mapping between the skeleton points and the surface points
            /*
            std::ofstream output(out_filename);
            // access SDF values (with constant-complexity)
            for(auto facet = mesh.facets_begin(); facet != mesh.facets_end(); ++facet){
                output << sdf_property_map[facet] << " ";
            }
            output << std::endl;
            */
            // create a property-map for segment-ids
            std::vector<std::size_t> segment_ids(self.size_of_facets());
            Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);

            std::size_t num_segments = CGAL::segmentation_from_sdf_values(
                self, sdf_property_map, segment_property_map);
            /*
            // access segment-ids (with constant-complexity)
            for(auto facet = mesh.facets_begin(); facet != mesh.facets_end(); ++facet){
                output << segment_property_map[facet] << " ";
            }
            output << std::endl;
            */
            return py::make_tuple(sdf_property_map, segment_property_map);
        })
        ;
}
