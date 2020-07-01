from cgal_pybind import Polyhedron, Point_3
import trimesh
import numpy as np

def test_polyhedron_contract():
    mesh = trimesh.primitives.Capsule(transform=np.array([[1, 0, 0, 2],
                                                          [0, 1, 0, 2],
                                                          [0, 0, 1, 2],
                                                          [0, 0, 0, 1], ]))
    # TODO: add_vertices/add_faces should be more efficient
    mesh_vertices = [Point_3(v[0], v[1], v[2]) for v in mesh.vertices]
    mesh_face_indices = [tuple(f) for f in mesh.faces]

    polyhedron = Polyhedron()
    polyhedron.build(mesh_vertices, mesh_face_indices)
    vertices, edges, correspondence, = polyhedron.contract()
    vertices = np.asarray(vertices).reshape((-1, 3))
    edges = np.asarray(edges)
    assert vertices.shape[0] > edges.shape[0]
    assert len(correspondence.keys()) == vertices.shape[0]


def test_polyhedron_segmentation():
    mesh = trimesh.primitives.Capsule(transform=np.array([[1, 0, 0, 2],
                                                          [0, 1, 0, 2],
                                                          [0, 0, 1, 2],
                                                          [0, 0, 0, 1], ]))
    # TODO: add_vertices/add_faces should be more efficient
    mesh_vertices = [Point_3(v[0], v[1], v[2]) for v in mesh.vertices]
    mesh_face_indices = [tuple(f) for f in mesh.faces]

    polyhedron = Polyhedron()
    polyhedron.build(mesh_vertices, mesh_face_indices)
    segment_result = np.array(polyhedron.segmentation()).reshape(2, -1)

    assert segment_result.shape == np.ndarray((2, 8064)).shape

