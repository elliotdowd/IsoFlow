import numpy as np
import numpy.linalg as la
import meshpy

def gen_nacamesh(domain, airfoil):

    def round_trip_connect(seq):
        result = []
        for i in range(len(seq)):
            result.append((i, (i+1)%len(seq)))
        return result
 
    pt_back = np.array([1,0])
 
    def max_area(pt):
        x = pt[0]

        pt_front = (-airfoil.L/2, 0)
        pt_back = (airfoil.L/2, 0)
 
        if x < -airfoil.L/2:
            return 5e-2*la.norm(pt-pt_front)**2 + 1e-2 # 2e-2, 4e-3
        elif x > airfoil.L/2:
            return 5e-2*la.norm(pt-pt_back)**2 + 1e-2 # 2e-2, 4e-3
        else:
            return 5e-2*pt[1]**2 + 1e-2 # 2e-2, 4e-3
 
    def needs_refinement(vertices, area):
        barycenter =  sum(np.array(v) for v in vertices)/3
        return bool(area > max_area(barycenter))
 
    from meshpy.naca import get_naca_points
    points = get_naca_points(naca_digits=airfoil.naca, number_of_points=airfoil.M)

    for k, pt in enumerate( points ):
        points[k][0] = pt[0]*airfoil.L - airfoil.L/2
        points[k][1] = pt[1]*airfoil.L
 
    from meshpy.geometry import GeometryBuilder, Marker
    from meshpy.triangle import write_gnuplot_mesh
 
    builder = GeometryBuilder()
    builder.add_geometry(points=points,
            facets=round_trip_connect(points),
            facet_markers=domain.profile_marker)

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

    # create list of face tags 
    mesh.face_tags = []
    for face in enumerate(fvi2fm.values()):
        # tag for specified face
        tag = face[1]

        for key in enumerate(domain.face_marker_to_tag.keys()):
            if tag == key:
                tag = face_marker_to_tag[key]

        mesh.face_tags.append( tag )


    return mesh


def gen_biconvexmesh(domain, airfoil):

    def round_trip_connect(seq):
        result = []
        for i in range(len(seq)):
            result.append((i, (i+1)%len(seq)))
        return result
 
    pt_back = np.array([1,0])
 
    def max_area(pt):
        x = pt[0]

        pt_front = (-airfoil.L/2, 0)
        pt_back = (airfoil.L/2, 0)
 
        if x < -airfoil.L/2:
            return 5e-2*la.norm(pt-pt_front)**2 + 1e-2 # 2e-2, 4e-3
        elif x > airfoil.L/2:
            return 5e-2*la.norm(pt-pt_back)**2 + 1e-2 # 2e-2, 4e-3
        else:
            return 5e-2*pt[1]**2 + 1e-2 # 2e-2, 4e-3
 
    def needs_refinement(vertices, area):
        barycenter =  sum(np.array(v) for v in vertices)/3
        return bool(area > max_area(barycenter))
 
    points = airfoil.vertices

    for k, pt in enumerate( points ):
        points[k][0] = pt[0]*airfoil.L - airfoil.L/2
        points[k][1] = pt[1]*airfoil.L
 
    from meshpy.geometry import GeometryBuilder, Marker
    from meshpy.triangle import write_gnuplot_mesh
 
    builder = GeometryBuilder()
    builder.add_geometry(points=points,
            facets=round_trip_connect(points),
            facet_markers=domain.profile_marker)

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

    # create list of face tags 
    mesh.face_tags = []
    for face in enumerate(fvi2fm.values()):
        # tag for specified face
        tag = face[1]

        for key in enumerate(domain.face_marker_to_tag.keys()):
            if tag == key:
                tag = face_marker_to_tag[key]

        mesh.face_tags.append( tag )


    return mesh



def gen_circlemesh(domain, circle):

    def round_trip_connect(seq):
        result = []
        for i in range(len(seq)):
            result.append((i, (i+1)%len(seq)))
        return result
 
    pt_back = np.array([1,0])
 
    def max_area(pt):
        x = pt[0]

        pt_front = (-circle.L/2, 0)
        pt_back = (circle.L/2, 0)
 
        if x < -circle.L/2:
            return 2e-2*la.norm(pt-pt_front)**2 + 4e-3
        elif x > circle.L/2:
            return 2e-2*la.norm(pt-pt_back)**2 + 4e-3
        else:
            return 2e-2*pt[1]**2 + 4e-3
 
    def needs_refinement(vertices, area):
        barycenter =  sum(np.array(v) for v in vertices)/3
        return bool(area > max_area(barycenter))
 
 
    from meshpy.geometry import GeometryBuilder, Marker, make_circle
    from meshpy.triangle import write_gnuplot_mesh
 
    builder = GeometryBuilder()
    meshpy.geometry.make_circle( circle.L/2, (0,0) )
    builder.wrap_in_box( (domain.L, domain.h), (domain.M, domain.N) )
 
    from meshpy.triangle import MeshInfo, build
    mesh_mod = builder.mesher_module()
    mi = mesh_mod.MeshInfo()
    # builder.set(mi)
    mi.set_holes( [0,0] )

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

    # create list of face tags 
    mesh.face_tags = []
    for face in enumerate(fvi2fm.values()):
        # tag for specified face
        tag = face[1]

        for key in enumerate(domain.face_marker_to_tag.keys()):
            if tag == key:
                tag = face_marker_to_tag[key]

        mesh.face_tags.append( tag )


    return mesh
