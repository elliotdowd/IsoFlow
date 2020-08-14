import numpy as np
from python.finite_volume.helper import thermo
from python.boundary.unstruct_boundary_cond import unstruct_boundary_cond, unstruct_covariant

# initialize Q vector given simulation parameters
def init_state(domain, mesh, boundary, parameters, gas):

    # specific heat ratio calculation
    gas.Cp = np.zeros( (domain.M+2, domain.N+2 ) ) + gas.Cp_fn( gas.gamma_p, gas.Cp_p, gas.theta, parameters.T_out )
    gas.Cv = np.zeros( (domain.M+2, domain.N+2 ) ) + gas.Cv_fn( gas.gamma_p, gas.Cv_p, gas.theta, parameters.T_out )

    # initialize Q state vector at each point
    Q = np.zeros((domain.M+2, domain.N+2, 4), dtype='float', order='F')
    #soln_vars.init_q(Q, parameters.p_in, parameters.T_in, parameters.M_in, gas.R, gas.gamma, domain.M+2, domain.M+2)

    Q[:,:,0] = parameters.p_out / (gas.R_p * parameters.T_out)
    Q[:,:,1] = Q[:,:,0] * parameters.M_in * np.cos(domain.alpha) * np.sqrt(gas.gamma_fn(gas.Cp, gas.Cv)*parameters.p_out/Q[:,:,0])
    Q[:,:,2] = Q[:,:,0] * parameters.M_in * np.sin(domain.alpha) * np.sqrt(gas.gamma_fn(gas.Cp, gas.Cv)*parameters.p_out/Q[:,:,0])
    Q[:,:,3] = thermo.calc_rho_et(parameters.p_out, Q[:,:,0], Q[:,:,1]/Q[:,:,0], Q[:,:,2]/Q[:,:,0], gas.gamma_fn(gas.Cp, gas.Cv))

    class state:
        pass

    state.u = Q[:,:,1]/Q[:,:,0]
    state.v = Q[:,:,2]/Q[:,:,0]
    state.Q = Q
    state.Qn = Q
    state.p = thermo.calc_p( Q[:,:,0], Q[:,:,3], Q[:,:,1]/Q[:,:,0], Q[:,:,2]/Q[:,:,0], gas.gamma_fn(gas.Cp, gas.Cv) )
    state.T = state.p / (gas.R_p * Q[:,:,0])
    state.Mach = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 ) / \
                           thermo.calc_c( state.p, state.Q[:,:,0], gas.gamma_fn(gas.Cp, gas.Cv) )
    state.vel = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 )
    state.p0 = (1+((gas.gamma_fn(gas.Cp, gas.Cv)-1)/2)*state.Mach**2)**(gas.gamma_fn(gas.Cp, gas.Cv)/(gas.gamma_fn(gas.Cp, gas.Cv)-1)) * state.p
    state.res = np.array([1, 1, 1, 1])
    state.T0 = (1+((gas.gamma_fn(gas.Cp, gas.Cv)-1)/2)*state.Mach**2) * state.T
    state.n = 0

    # find temperature at centroids
    c = thermo.calc_c( state.p, state.Q[:,:,0], gas.gamma_fn(gas.Cp, gas.Cv) )
    state.c = c

    # boundary conditions
    state = enforce_bc(domain, mesh, boundary, parameters, state, gas)

    state = covariant(mesh, state)

    return state


# initialization for unstructured meshes
def init_unstruct_state( mesh, parameters, gas ):

    # initialize state vector at each element
    Q = np.zeros( [ len(mesh.elements), 4 ] )
    c = np.zeros( len(mesh.elements) ) + np.sqrt( gas.gamma*parameters.p_in / (parameters.p_in / (gas.R * parameters.T_in)) )
    u = np.zeros( len(mesh.elements) ) + parameters.M_in * c
    v = np.zeros( len(mesh.elements) ) + 0
    p = np.zeros( len(mesh.elements) ) + parameters.p_in
    T = np.zeros( len(mesh.elements) ) + parameters.T_in

    Q[:, 0] = parameters.p_in / (gas.R * parameters.T_in)
    Q[:, 1] = Q[:, 0] * u
    Q[:, 2] = Q[:, 0] * v
    Q[:, 3] = parameters.p_in/(gas.gamma-1) + (1/2)*Q[:, 0]*( u**2 + v**2 )

    class state:
        pass

    state.Q = Q
    state.u = u
    state.v = v
    state.p = p
    state.T = T
    state.c = c
    state.n = 0

    # create boundary "halo" cell state vector and variables with indices pointing to related cell face, Q
    state.Qbound = np.zeros( [len( mesh.bdry_ind[0] ), 4] )
    state.pbound = np.zeros( len( mesh.bdry_ind[0] ) )
    state.Tbound = np.zeros( len( mesh.bdry_ind[0] ) )
    state.ubound = np.zeros( len( mesh.bdry_ind[0] ) )
    state.vbound = np.zeros( len( mesh.bdry_ind[0] ) )
    state.Ubound = np.zeros( [len( mesh.bdry_ind[0] ), 2] )

    # enforce boundary conditions
    state = unstruct_boundary_cond( mesh, state, parameters, gas )

    # calculate covariant velocities 
    state.U = np.zeros( [len(mesh.faces), 2] )
    state = unstruct_covariant( mesh, state )

    return mesh, state