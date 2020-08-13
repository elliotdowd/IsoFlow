import numpy as np
import numpy.linalg as la

def gen_nacamesh(domain, airfoil):

    def round_trip_connect(seq):
        result = []
        for i in range(len(seq)):
            result.append((i, (i+1)%len(seq)))
        return result
 
    pt_back = np.array([1,0])
 
    def max_area(pt):
        x = pt[0]
 
        if x < -airfoil.L/2:
            return 2e-2*la.norm(pt)**2 + 1e-5
        elif x > airfoil.L/2:
            return 2e-2*la.norm(pt-pt_back)**2 + 1e-5
        else:
            return 2e-2*pt[1]**2 + 1e-5
 
    def needs_refinement(vertices, area):
        barycenter =  sum(np.array(v) for v in vertices)/3
        return bool(area > max_area(barycenter))
 
    # from meshpy.naca import get_naca_points
    # points = get_naca_points(naca_digits=airfoil.naca, number_of_points=airfoil.N)
    points = airfoil.vertices
 
    from meshpy.geometry import GeometryBuilder, Marker
    from meshpy.triangle import write_gnuplot_mesh
 
    profile_marker = Marker.FIRST_USER_MARKER
    builder = GeometryBuilder()
    builder.add_geometry(points=points,
            facets=round_trip_connect(points),
            facet_markers=profile_marker)

    builder.wrap_in_box((domain.L, domain.h), (domain.M, domain.N))
 
    from meshpy.triangle import MeshInfo, build
    mi = MeshInfo()
    builder.set(mi)
    mi.set_holes( [builder.center()] )

 
    mesh = build(mi, refinement_func=needs_refinement,
            #allow_boundary_steiner=False,
            generate_faces=True)

    # allocate memory for element volumes list
    mesh.element_volumes.setup()
    mesh.normals.setup()
 
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
