#include <string>
#include <vector>
#include <utility> // std::pair

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

using Point_3 = CGAL::Point_3<CGAL::Exact_predicates_inexact_constructions_kernel>;

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

using IntIndices = std::vector<std::tuple<int, int, int>>;
// A modifier creating a Polyhedron from vectors of vertex and faces(vertex Ids).
template<class HDS>
class polyhedron_builder : public CGAL::Modifier_base<HDS> {
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = CGAL::Polyhedron_3<Kernel> ;
using HalfedgeDS = Polyhedron::HalfedgeDS;

public:
    const std::vector<Point_3>&vertexes;
    const IntIndices &tris;
    polyhedron_builder( const std::vector<Point_3> &_vertexes, const IntIndices &_tris ) : vertexes(_vertexes), tris(_tris) {}
    void operator()( HDS& hds) {
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        // create a cgal incremental builder
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( vertexes.size(), tris.size() );
        // add the polyhedron vertices
        for(const auto& p : vertexes) {
            B.add_vertex( Point(p));
        }
        // add the polyhedron triangles
        for(auto [v0, v1, v2] : tris) {
            B.begin_facet();
            B.add_vertex_to_facet(v0);
            B.add_vertex_to_facet(v1);
            B.add_vertex_to_facet(v2);
            B.end_facet();
        }
        // finish up the surface
        B.end_surface();
    }

};

/*-----------------------------------------------------*/


/*---------------------------------------------------*/
namespace py = pybind11;
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
//using Point_3 = Kernel::Point_3;
using Polyhedron = CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>;


py::tuple contract(Polyhedron& polyhendron) {
    /**
     * Convert an interactively contracted skeleton to a skeleton curve
     *
     * @param Polyhendron Polyhedron reference.
     * @return py::tuple containing to vector (view as Python list):
           return_vertices: vector containing skeleton curve vertices describe as 3D point
           return_edges: vector containing skeleton curve edges
           return_correspondence: map containing skeleton mapping between vertices id and vector of surface id.
     */
    using Skeletonization = CGAL::Mean_curvature_flow_skeletonization<Polyhedron> ;
    using Skeleton  = Skeletonization::Skeleton;
    Skeletonization mean_curve_skeletonizer(polyhendron);
    // Iteratively calls CGAL::contract() until the change of surface area of the meso-skeleton after one iteration
    // is smaller than area_variation_factor()*original_area where original_area is the area of the input triangle mesh,
    // or if the maximum number of iterations has been reached.
    mean_curve_skeletonizer.contract_until_convergence();

    Skeletonization::Skeleton skeleton;
    // Converts the contracted Polyhedron to a skeleton curve.
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

    // associate indices to the vertices using the "id()" field of the vertex.
    using vertex_descriptor = boost::graph_traits<Polyhedron>::vertex_descriptor ;
    typedef boost::graph_traits<Polyhedron>::vertex_iterator   vertex_iterator;
    typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
    vertex_iterator vb, ve;
    int index = 0;
    // boost::tie assigns the first and second element of the std::pair
    // returned by boost::vertices to the variables vb and ve
     for(boost::tie(vb,ve)=vertices(polyhendron); vb!=ve; ++vb ){
        vertex_descriptor  vd = *vb;
        vd->id() = index++;
    }

    // Output skeleton points and the corresponding surface points
    std::map<size_t, std::vector<size_t>> return_correspondence;
    for (const auto& v : CGAL::make_range(vertices(skeleton)))
    {
        std::vector<size_t> vector_surface_id;
        for(const auto& vd : skeleton[v].vertices ){
            vector_surface_id.push_back(get(boost::vertex_index, polyhendron, vd));
        }
        return_correspondence.insert(std::pair<size_t,std::vector<size_t>>(v,vector_surface_id));
    }
    return py::make_tuple(return_vertices, return_edges, return_correspondence);
}


py::tuple segmentation(Polyhedron& polyhendron) {
    /**
     * assign a unique id to each facet in the mesh.
     * computes property-map for SDF values.
     * computes property-map for segment-ids.
     * @param Polyhendron Polyhedron reference.
     * @return py::tuple containing to vector (view as Python list):
           return_sdf_property_map: vector containing Shape Diameter Function property map
           return_segment_property_map: vector containing segment property map
     */
    // assign id field for each facet
    std::size_t facet_id = 0;
    for(auto facet = polyhendron.facets_begin(); facet != polyhendron.facets_end(); ++facet) {
        facet->id() = facet_id;
        ++facet_id;
    }
    /* The Shape Diameter Function provides a connection between the surface mesh and the
    * volume of the subtended 3D bounded object. More specifically, the SDF is a scalar-valued function defined on
    * facets of the mesh that measures the corresponding local object diameter.
    */
    // create a property-map for signed distance function values
    std::vector<double> sdf_values_vec(polyhendron.size_of_facets());
    Facet_with_id_pmap<double> sdf_property_map(sdf_values_vec);
    // computing the Shape Diameter Function over a surface mesh into sdf_property_map .
    CGAL::sdf_values(polyhendron, sdf_property_map);

    // Fill return_sdf_property_map
    std::vector<std::tuple<double>> return_sdf_property_map(polyhendron.size_of_facets());
    for(auto facet = polyhendron.facets_begin(); facet != polyhendron.facets_end(); ++facet){
        return_sdf_property_map.emplace_back(sdf_property_map[facet]);
     }

    // create a property-map for segment-ids
    std::vector<std::size_t> segment_ids(polyhendron.size_of_facets());
    Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);

    // computes the segmentation of a surface mesh given an SDF value per facet.
    /* CGAL::segmentation_from_sdf_values function computes the segmentation of a surface mesh given an SDF value per facet.
    * This function fills a property map which associates a segment-id (in [0, number of segments -1])
    * or a cluster-id (in [0, number_of_clusters -1]) to each facet. A segment is a set of connected facets
    * which are placed under the same cluster (see Figure 66.5).
    */
    CGAL::segmentation_from_sdf_values(polyhendron, sdf_property_map, segment_property_map);

    // Fill return_segment_property_map
    std::vector<std::tuple<size_t>> return_segment_property_map(polyhendron.size_of_facets());
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
        .def("build",[](Polyhedron& self, const std::vector<Point_3>& vertexes,const IntIndices& indices) {

            polyhedron_builder<HalfedgeDS> builder( vertexes, indices );
            self.delegate( builder );
        })
    ;
}
