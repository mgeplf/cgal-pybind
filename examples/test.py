# XXX: THIS IS A TERRIBLE HACK
import sys
sys.path.insert(0, '../build/')

import cgal_pybind

import numpy as np
import trimesh

mesh = trimesh.primitives.Capsule(transform=np.array([[1, 0, 0, 2],
                                                      [0, 1, 0, 2],
                                                      [0, 0, 1, 2],
                                                      [0, 0, 0, 1],]))
surface_mesh = cgal_pybind.SurfaceMesh()

#TODO: add_vertices/add_faces should be more efficient
surface_mesh.add_vertices([cgal_pybind.Point_3(v[0], v[1], v[2]) for v in mesh.vertices])
surface_mesh.add_faces([tuple(f) for f in mesh.faces])

print('number_of_vertices: ', surface_mesh.number_of_vertices())
print('number_of_faces: ', surface_mesh.number_of_faces())
print('area: ', surface_mesh.area())

def contract():
    vertices, edges = surface_mesh.contract()
    print('vertices\n', vertices)
    print('edges\n', edges)

    #XXX: need to shape this right on exit
    vertices = np.array(vertices).reshape((-1, 3))

    print('number_of_vertices: ', len(vertices))
    print('number_of_faces: ', len(edges))

    # display the skeleton shifted over, so one can see it
    path_visual = trimesh.load_path(vertices[edges, :] + [2, 0, 0])
    scene = trimesh.Scene([path_visual, mesh ])
    scene.show()


def authalic():
    edges = surface_mesh.authalic()
    edges = list(zip(edges[:-1], edges[1:]))

    path_visual = trimesh.load_path(mesh.vertices[edges, :] + [2, 0, 0])
    scene = trimesh.Scene([path_visual, mesh ])
    scene.show()


#contract()
authalic()
