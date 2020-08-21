import numpy as np

# inverse metrics, projected cell face areas, cell volumes, centroids
def cellmetrics(xx, yy, domain):

    import numpy as np
    import python.mesh.metrics.cells as cells

    # format grid variable for input into fortran subroutines
    grid = np.array((xx, yy), dtype='float', order='F')
    grid = np.array(grid).reshape(domain.M+3,domain.N+3,2)

    # calculate projected cell face areas and cell face areas 
    # S_proj = [S_zeta.x, S_zeta.y, S_eta.x, S_eta.y S_zeta, S_eta]

    s_proj = np.zeros((domain.M+2, domain.N+2, 6), dtype='float', order='F')
    s_proj[:,:,0] = yy[1:domain.M+3, 1:domain.N+3] - yy[1:domain.M+3, 0:domain.N+2]
    s_proj[:,:,1] = -( xx[1:domain.M+3, 1:domain.N+3] - xx[1:domain.M+3, 0:domain.N+2] )
    s_proj[:,:,2] = -( yy[1:domain.M+3, 1:domain.N+3] - yy[0:domain.M+2, 1:domain.N+3] )
    s_proj[:,:,3] = xx[1:domain.M+3, 1:domain.N+3] - xx[0:domain.M+2, 1:domain.N+3]
    s_proj[:,:,4] = np.sqrt( s_proj[:,:,0]**2 + s_proj[:,:,1]**2 )
    s_proj[:,:,5] = np.sqrt( s_proj[:,:,2]**2 + s_proj[:,:,3]**2 )

    # calculate cell areas
    area = np.zeros((domain.M+2, domain.N+2), dtype='float', order='F')
    cells.calc_cellarea(4, xx, yy, area, domain.M+2, domain.N+2)

    # calculate cell centroids
    xxc = np.zeros((domain.M+2, domain.N+2), dtype='float', order='F')
    yyc = np.zeros((domain.M+2, domain.N+2), dtype='float', order='F')

    cells.calc_cellcentroids(xx, yy, xxc, yyc, area, domain.M+2, domain.N+2)

    # create mesh class
    class mesh:
        pass

    mesh.xx = xx
    mesh.yy = yy
    mesh.xxc = xxc
    mesh.yyc = yyc
    mesh.dV = area
    mesh.s_proj = s_proj

    return mesh


# face normals and cell volumes for unstructured mesh
def unstruct_cellmetrics( mesh ):

    mesh.face_areas = []
    mesh.dV = np.zeros( len(mesh.elements) )
    mesh.Sx = np.zeros( len(mesh.faces) )
    mesh.Sy = np.zeros( len(mesh.faces) )
    mesh.dS = np.zeros( len(mesh.faces) )
    mesh.n = np.zeros( [len(mesh.faces), 2] )

    for i, elem in enumerate( mesh.elements ):

        # triangular elements
        if len( elem ) == 3:
            x1 = mesh.points[elem[0]][0]
            y1 = mesh.points[elem[0]][1]
            x2 = mesh.points[elem[1]][0]
            y2 = mesh.points[elem[1]][1]
            x3 = mesh.points[elem[2]][0]
            y3 = mesh.points[elem[2]][1]

            mesh.element_volumes[i] = (1/2) * ( (x1-x2)*(y1+y2) + (x2-x3)*(y2+y3) + (x3-x1)*(y3+y1) )
            # also save to numpy array
            mesh.dV[i] = mesh.element_volumes[i]

        # quadrilateral elements
        elif len( elem ) == 4:
            x1 = mesh.points[elem[0]][0]
            y1 = mesh.points[elem[0]][1]
            x2 = mesh.points[elem[1]][0]
            y2 = mesh.points[elem[1]][1]
            x3 = mesh.points[elem[2]][0]
            y3 = mesh.points[elem[2]][1]
            x4 = mesh.points[elem[3]][0]
            y4 = mesh.points[elem[3]][1]

            mesh.element_volumes[i] = (1/2) * ( (x1-x3)*(y2-y4) + (x4-x2)*(y1-y3) )
            mesh.dV[i] = mesh.element_volumes[i]

    # face area magnitude and normal vector calculations
    for j, pts in enumerate( mesh.faces ):

        # outward pointing face vector components
        Sx = mesh.points[pts[1]][1] - mesh.points[pts[0]][1]
        Sy = mesh.points[pts[0]][0] - mesh.points[pts[1]][0]

        dS = np.sqrt( Sx**2 + Sy**2 )

        mesh.face_areas.append( dS )
        mesh.Sx[j] = Sx
        mesh.Sy[j] = Sy
        mesh.dS[j] = dS

        mesh.normals[j] = ( Sx/dS, Sy/dS )
        mesh.n[j,:] = mesh.normals[j]

    return mesh


# determine neighbors of each cell face, pointers from elements to faces (left and right cell states)
def find_facepairs( mesh ):

    mesh.face_pairs = np.zeros( [len(mesh.faces), 2] )
    mesh.face_pts = np.zeros( [len(mesh.faces), 2] )
    mesh.elem_to_face = np.zeros( [len(mesh.elements), 3], dtype=int ) - 1
    elem_face_allocated = np.zeros( len(mesh.elements), dtype=int )

    # loop through faces to convert face points to numpy array
    for i, face in enumerate( mesh.faces ):

        mesh.face_pts[i,:] = face

    # mask for boundary faces
    mesh.bdry_ind = np.nonzero(mesh.face_tags)

    # loop through element neighbors to find matching face for each
    for i, neighbor in enumerate( mesh.neighbors ):

        for j, elem in enumerate( neighbor ):

            if elem == -1:

                for k in range( 0, len(mesh.elements[i]) ):
                    # find points which belong to face on boundary
                    face_pts = np.array( [mesh.elements[i][k], mesh.elements[i][np.mod(k+1, len(mesh.elements[i]))] ] )

                    # find where face points == data structure face definition
                    face_id_mask = np.logical_or( face_pts == mesh.face_pts, np.flipud(face_pts) ==  mesh.face_pts )
                    face_id = np.where( np.logical_and( face_id_mask[:,0] == True, face_id_mask[:,1] == True ) )
                    
                    if np.any( np.isin( np.nonzero(mesh.face_tags), face_id ) ):
                        mesh.face_pairs[face_id,:] = ( i, -mesh.face_tags[int(face_id[0])] )

                        if elem_face_allocated[i] < 3:
                            # point from elements to faces
                            mesh.elem_to_face[i,elem_face_allocated[i]] = face_id[0]
                            # tracking of how many element -> face pointers have been allocated for a given element
                            elem_face_allocated[i] = int(elem_face_allocated[i] + 1)
                
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

                # inside(left) element, then right(outside) element in mesh.face_pairs
                mesh.face_pairs[face_id,:] = (elem, i)

                # point from elements to faces
                mesh.elem_to_face[i,elem_face_allocated[i]] = face_id[0]
                # tracking of how many element->face pointers have been allocated for a given element
                elem_face_allocated[i] = int(elem_face_allocated[i] + 1)

    return mesh

