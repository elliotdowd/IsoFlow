import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

from collections import defaultdict

import meshpy

def make_nacamesh():

    def round_trip_connect(seq):
        result = []
        for i in range(len(seq)):
            result.append((i, (i+1)%len(seq)))
        return result
 
    pt_back = np.array([1,0])
 
    #def max_area(pt):
        #max_area_front = 1e-2*la.norm(pt)**2 + 1e-5
        #max_area_back = 1e-2*la.norm(pt-pt_back)**2 + 1e-4
        #return min(max_area_front, max_area_back)
 
    def max_area(pt):
        x = pt[0]
 
        if x < 0:
            return 2e-2*la.norm(pt)**2 + 1e-5
        elif x > 1:
            return 2e-2*la.norm(pt-pt_back)**2 + 1e-5
        else:
            return 2e-2*pt[1]**2 + 1e-5
 
    def needs_refinement(vertices, area):
        barycenter =  sum(np.array(v) for v in vertices)/3
        return bool(area > max_area(barycenter))
 
    from meshpy.naca import get_naca_points
    points = get_naca_points(naca_digits="24012", number_of_points=120)
 
    from meshpy.geometry import GeometryBuilder, Marker
    from meshpy.triangle import write_gnuplot_mesh
 
    profile_marker = Marker.FIRST_USER_MARKER
    builder = GeometryBuilder()
    builder.add_geometry(points=points,
            facets=round_trip_connect(points),
            facet_markers=profile_marker)

    # add circle geometry
    builder.add_geometry(*meshpy.geometry.make_circle(0.05, center=(1.5, 0), subdivisions=60, marker=1001))

    builder.wrap_in_box((12, 9), (16, 12))
 
    from meshpy.triangle import MeshInfo, build
    mi = MeshInfo()
    builder.set(mi)
    mi.set_holes(np.vstack([builder.center(), [1.5, 0]]) )
    # mi.set_holes( [builder.center()] )

 
    mesh = build(mi, refinement_func=needs_refinement,
            #allow_boundary_steiner=False,
            generate_faces=True)

    # allocate memory for element volumes list
    mesh.element_volumes.setup()
 
    # write_gnuplot_mesh("naca.dat", mesh)
 
    print( "%d vertices" % len(mesh.points))
    print( "%d elements" % len(mesh.elements) )
 
    fvi2fm = mesh.face_vertex_indices_to_face_marker
 
    face_marker_to_tag = {
            profile_marker: "noslip",
            Marker.MINUS_X: "inflow",
            Marker.PLUS_X: "outflow",
            Marker.MINUS_Y: "inflow",
            Marker.PLUS_Y: "inflow"
            }

    # create list of face tags 
    mesh.face_tags = []
    for face in enumerate(fvi2fm.values()):
        # tag for specified face
        tag = face[1]

        for key in enumerate(face_marker_to_tag.keys()):
            if tag == key:
                tag = face_marker_to_tag[key]

        mesh.face_tags.append( tag )


    return mesh

# call mesh creation function
mesh = make_nacamesh()


# calculate element areas and cell face areas
def calc_cell_metrics( mesh ):

    i = 0
    for elem in mesh.elements:
        x1 = mesh.points[elem[0]][0]
        y1 = mesh.points[elem[0]][1]
        x2 = mesh.points[elem[1]][0]
        y2 = mesh.points[elem[1]][1]
        x3 = mesh.points[elem[2]][0]
        y3 = mesh.points[elem[2]][1]

        mesh.element_volumes[i] = (1/2) * ( (x1-x2)*(y1+y2) + (x2-x3)*(y2+y3) + (x3-x1)*(y3+y1) )
        i = i+1

    return mesh

mesh = calc_cell_metrics(mesh)


# determine neighbors of each cell face
def find_facepairs(mesh):

    mesh.face_pairs = np.zeros( [len(mesh.faces), 2] )
    mesh.face_pts = np.zeros( [len(mesh.faces), 2] )

    # loop through faces to convert face points to numpy array
    for i, face in enumerate( mesh.faces ):

        mesh.face_pts[i,:] = face


    for i, neighbor in enumerate( mesh.neighbors ):

        for j, elem in enumerate( neighbor ):

            if elem == -1:
                pass
            else:
                elem_pts = mesh.elements[i]
                neigh_pts = mesh.elements[elem]

                # find points which belong to face between neighbors
                combine = np.hstack( [elem_pts, neigh_pts] )
                unq, unq_id, unq_ct = np.unique( combine, return_counts=True, return_inverse=True )
                mask = unq_ct > 1
                face_pts = unq[mask]

                # find where face points == data structure face definition
                face_id_mask = np.logical_or( face_pts == mesh.face_pts, np.flipud(face_pts) ==  mesh.face_pts )
                face_id = np.where( np.logical_and( face_id_mask[:,0] == True, face_id_mask[:,1] == True ) )

            mesh.face_pairs[face_id,:] = (i, elem)

    return mesh

mesh = find_facepairs(mesh)


# plot unstructured mesh
for face in enumerate( mesh.faces ):
    i = face[0]
    pts = []
    for pt in enumerate( face[1] ):
        pts.append( mesh.points[pt[1]] )

    tag = mesh.face_tags[i]

    if tag == 0:
        plt.plot( (pts[0][0], pts[1][0]), (pts[0][1], pts[1][1]), 'g-', linewidth=0.15 )
    else:
        plt.plot( (pts[0][0], pts[1][0]), (pts[0][1], pts[1][1]), 'k-', linewidth=0.5 )

# plot mesh lines
plt.axis('equal')
plt.show()

mesh
