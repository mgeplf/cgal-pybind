#include <string>
#include <vector>
#include <utility> // std::pair

//#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


// Property map associating a facet with an integer as id to an
// element in a vector stored internally
using key_type =  CGAL::Polyhedron_3<
    CGAL::Exact_predicates_inexact_constructions_kernel,
    CGAL::Polyhedron_items_with_id_3
>::Facet_const_handle;
template<class ValueType>
struct Facet_with_id_pmap
        : public boost::put_get_helper<ValueType&, Facet_with_id_pmap<ValueType> >
{
    using category = boost::lvalue_property_map_tag;
    Facet_with_id_pmap(std::vector<ValueType>& internal_vector) : internal_vector(internal_vector) { }
    ValueType& operator[](key_type key) const
    { return internal_vector[key->id()]; }

private:
    std::vector<ValueType>& internal_vector;
};



// A modifier creating a Polyhedron from vectors of vertex and faces(vertex Ids).
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS> {
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = CGAL::Polyhedron_3<Kernel> ;
using HalfedgeDS = Polyhedron::HalfedgeDS;
public:
    std::vector<double> &coords;
    std::vector<int>    &tris;
    polyhedron_builder( std::vector<double> &_coords, std::vector<int> &_tris ) : coords(_coords), tris(_tris) {}
    void operator()( HDS& hds) {
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;

        // create a cgal incremental builder
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( coords.size()/3, tris.size()/3 );

        // add the polyhedron vertices
        for( int i=0; i<(int)coords.size(); i+=3 ){
        B.add_vertex( Point( coords[i+0], coords[i+1], coords[i+2] ) );
        }
        // add the polyhedron triangles
        for( int i=0; i<(int)tris.size(); i+=3 ){
            B.begin_facet();
            B.add_vertex_to_facet( tris[i+0] );
            B.add_vertex_to_facet( tris[i+1] );
            B.add_vertex_to_facet( tris[i+2] );
            B.end_facet();
        }
        // finish up the surface
        B.end_surface();
    }
};
/*---------------------------------------------------*/

namespace py = pybind11;
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Polyhedron = CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>;


py::tuple contract(Polyhedron& polyhendron) {
    /**
     * Convert an interactively contracted skeleton to a skeleton curve
     *
     * @param Polyhendron Polyhedron reference.
     * @return py::tuple containing to vector (view as Python list):
           return_vertices: vector containing skeleton curve vertices describe as 3D point
           return_edges: vector containing skeleton curve edges
     */
    using Skeletonization = CGAL::Mean_curvature_flow_skeletonization<Polyhedron> ;
    using Skeleton  = Skeletonization::Skeleton;
    Skeletonization mean_curve_skeletonizer(polyhendron);
    // Iteratively calls CGAL::contract() until the change of surface area of the meso-skeleton after one iteration
    // is smaller than area_variation_factor()*original_area where original_area is the area of the input triangle mesh,
    // or if the maximum number of iterations has been reached.
    mean_curve_skeletonizer.contract_until_convergence();

    Skeletonization::Skeleton skeleton;
    // Converts the contracted surface mesh to a skeleton curve.
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
}


py::tuple segmentation(Polyhedron& polyhendron) {
    /**
     * assign a unique id to each facet in the mesh
     *
     * @param Polyhendron Polyhedron reference.
     * @return py::tuple containing to vector (view as Python list):
           return_sdf_property_map: vector containing Shape Diameter Function property map
           return_segment_property_map: vector containing segment property map
     */
    std::size_t facet_id = 0;
    for(auto facet = polyhendron.facets_begin(); facet != polyhendron.facets_end(); ++facet) {
        facet->id() = facet_id;
        ++facet_id;
    }
    // create a property-map for signed distance function values
    std::vector<double> sdf_values_vec(polyhendron.size_of_facets());
    Facet_with_id_pmap<double> sdf_property_map(sdf_values_vec);

    // computing the Shape Diameter Function over a surface mesh into sdf_property_map .
    CGAL::sdf_values(polyhendron, sdf_property_map);

    // mapping between the skeleton points and the surface points
    // Fill return_sdf_property_map
    using vector_double = std::vector<std::tuple<double>>;
    vector_double return_sdf_property_map(polyhendron.size_of_facets());
    for(auto facet = polyhendron.facets_begin(); facet != polyhendron.facets_end(); ++facet){
        return_sdf_property_map.emplace_back(sdf_property_map[facet]);
     }

    // create a property-map for segment-ids
    std::vector<std::size_t> segment_ids(polyhendron.size_of_facets());
    Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);

    // computes the segmentation of a surface mesh given an SDF value per facet.
    CGAL::segmentation_from_sdf_values(polyhendron, sdf_property_map, segment_property_map);
    // Fill return_segment_property_map
    using vector_size_t = std::vector<std::tuple<size_t>>;
    vector_size_t return_segment_property_map(polyhendron.size_of_facets());
    for(auto facet = polyhendron.facets_begin(); facet != polyhendron.facets_end(); ++facet){
            return_segment_property_map.emplace_back(segment_property_map[facet]);
    }

    return py::make_tuple(return_sdf_property_map, return_segment_property_map);
}


void bind_triangle_mesh(py::module& m)
    {
    using HalfedgeDS=Polyhedron::HalfedgeDS;
    /**
     * Create Polyhedron class and its methods
     *
     * @param Polyhendron Polyhedron.
     */
    py::class_<Polyhedron>(m, "Polyhedron",
    "CGAL CGAL::Polyhedron_3<Exact_predicates_inexact_constructions_kernel, Polyhedron_items_with_id_3> binding")
        .def(py::init())
        .def("size_of_vertices", &Polyhedron::size_of_vertices, "number of vertices")
        .def("size_of_facets", &Polyhedron::size_of_facets, "number of facet")
        .def("area", [](Polyhedron& self) {
             return CGAL::Polygon_mesh_processing::area(self);
        }, "polyhedron area")
        .def("load_from_off",[](Polyhedron& self, const std::string& full_path_name) {
            std::ifstream input(full_path_name);
            input >> self;
        })
        /*
         * Default policy for method: return_value_policy::automatic. This policy falls back to the policy
         * return_value_policy::take_ownership when the return value is a pointer.
         * Otherwise, it uses return_value_policy::move or return_value_policy::copy
         * for rvalue and lvalue references, respectively. See above for a description of what
         * all of these different policies do.
         */
        .def("contract", &contract, "Convert an interactively contracted skeleton to a skeleton curve")
        .def("segmentation",&segmentation, "Assign a unique id to each facet in the mesh")
        .def("build",[](Polyhedron& self, std::vector<double>& coords, std::vector<int>& tris) {
            polyhedron_builder<HalfedgeDS> builder( coords, tris );
            self.delegate( builder );
        })
    ;
}
