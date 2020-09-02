import numpy as np
from python.finite_volume.helper import thermo

def unstruct_boundary_cond( domain, mesh, state, parameters, gas ):
    
    for i, bdry_id in enumerate( mesh.bdry_ind[0] ):

        # numeric boundary type id
        type_id = mesh.face_markers[bdry_id]
        # pointers from faces to left (interior) and right (exterior) elements
        L = int( mesh.face_pairs[bdry_id,0] )
        R = int( mesh.face_pairs[bdry_id,1] )

        bound_type = domain.face_marker_to_tag[type_id]

        if bound_type == 'slip':

            # update wall values
            Qwall, pwall, Twall = slip_wall( state.Q[L,:], state.p[L], state.T[L], mesh.n[bdry_id,0], mesh.n[bdry_id,1], gas )

            state.Qbound[i,:] = np.transpose( Qwall )
            state.pbound[i] = pwall
            state.Tbound[i] = Twall

        elif bound_type == 'noslip':

            # update wall values
            Qwall, pwall, Twall = noslip_wall( state.Q[L,:], state.p[L], state.T[L], gas )

            state.Qbound[i,:] = np.transpose( Qwall )
            state.pbound[i] = pwall
            state.Tbound[i] = Twall
            
        elif bound_type == 'inflow':

            state.Qbound[i,0] = parameters.p_in / (gas.R * parameters.T_in)
            state.Qbound[i,1] = state.Qbound[i,0] * parameters.M_in * np.sqrt( gas.gamma * parameters.p_in / state.Qbound[i,0] )
            state.Qbound[i,2] = 0
            state.Qbound[i,3] = parameters.p_in/(gas.gamma-1) + \
                            (1/2)*state.Qbound[i,0]*( (state.Qbound[i,1]/state.Qbound[i,0])**2 + (state.Qbound[i,2]/state.Qbound[i,0])**2 )

            state.pbound[i] = parameters.p_in
            state.Tbound[i] = parameters.T_in

        elif bound_type == 'outflow':

            state.Qbound[i,:] = state.Q[L,:]

            state.pbound[i] = state.p[L]
            state.Tbound[i] = state.T[L]

        state.ubound = state.Qbound[:,1] / state.Qbound[:,0]
        state.vbound = state.Qbound[:,2] / state.Qbound[:,0]
        

    return state


# inviscid wall function
def noslip_wall( Q, p, T, gas ):

    # adiabatic wall
    Twall = T

    # no pressure gradient 
    pwall = p

    # u0 = -u1 and v0 = -v1
    Qwall = np.zeros( 4 )
    
    Qwall[0] = pwall / (gas.R * Twall)
    Qwall[1] = -Qwall[0] * Q[1]/Q[0]
    Qwall[2] = -Qwall[0] * Q[2]/Q[0]
    Qwall[3] = pwall/(gas.gamma-1) + (1/2)*Qwall[0]*( (Qwall[1]/Qwall[0])**2 + (Qwall[2]/Qwall[0])**2 )

    return Qwall, pwall, Twall


def slip_wall( Q, p, T, nx, ny, gas ):

    # adiabatic wall
    Twall = T

    # no pressure gradient 
    pwall = p

    # u0 = -u1 and v0 = -v1
    Qwall = np.zeros( 4 )

    # inviscid velocity projection
    transform = np.array( [ 1-2*nx**2,   -2*nx*ny, \
                           -2*nx*ny,      1-2*ny**2 ] ).reshape( [2, 2] )

    uv = np.matmul( transform, np.array( [Q[1]/Q[0], Q[2]/Q[0] ] ).reshape([2,1]) )
    
    Qwall[0] = pwall / (gas.R * Twall)
    Qwall[1] = uv[0] * Qwall[0]
    Qwall[2] = uv[1] * Qwall[0]
    # Qwall[3] = pwall/(gas.gamma-1) + (1/2)*Qwall[0]*( (Qwall[1]/Qwall[0])**2 + (Qwall[2]/Qwall[0])**2 )
    Qwall[3] = thermo.calc_rho_et( pwall, Qwall[0], Qwall[1]/Qwall[0], Qwall[2]/Qwall[0], gas.gamma )

    return Qwall, pwall, Twall


# convective velocity at each cell face
def unstruct_covariant( mesh, state ):

    # pointers from faces to left (interior) and right (exterior) elements
    L = mesh.face_pairs[:,0]
    R = mesh.face_pairs[:,1]

    # calculate average velocities at each cell face
    for i, face in enumerate( mesh.face_pairs ):

        # pointers from faces to left (interior) and right (exterior) elements
        L = int( mesh.face_pairs[i,0] )
        R = int( mesh.face_pairs[i,1] )

        # account for boundary cells
        if R < 0:
            bound_i = int(np.where(np.isin(mesh.bdry_ind,i))[1])
            u_L = np.array( [ state.u[int(mesh.face_pairs[i,0])], state.v[int(mesh.face_pairs[i,0])] ] ).reshape([1,2])
            u_R = np.array( [ state.ubound[bound_i], state.vbound[bound_i] ] ).reshape([1,2])

            # calculate covariant velocity at each cell face
            state.U[i,0] = np.dot( u_L, mesh.n[i,:] )
            state.Ubound[bound_i,1] = np.dot( u_R, mesh.n[i,:] )

        else:
            u_L = np.array( [ state.u[int(mesh.face_pairs[i,0])], state.v[int(mesh.face_pairs[i,0])] ] ).reshape([1,2])
            u_R = np.array( [ state.u[int(mesh.face_pairs[i,1])], state.v[int(mesh.face_pairs[i,1])] ] ).reshape([1,2])

            # calculate left and right covariant velocity at each cell face
            state.U[i,0] = np.dot( u_L, mesh.n[i,:] )
            state.U[i,1] = np.dot( u_R, mesh.n[i,:] )


    return state