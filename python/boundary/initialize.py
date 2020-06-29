

# initialize Q vector given simulation parameters
def init_state(domain, mesh, parameters, gas):

    import numpy as np
    from python.finite_volume.helper import thermo

    # initialize Q state vector at each point
    Q = np.zeros((domain.M+2, domain.N+2, 4), dtype='float', order='F')
    #soln_vars.init_q(Q, parameters.p_in, parameters.T_in, parameters.M_in, gas.R, gas.gamma, domain.M+2, domain.M+2)

    Q[:,:,0] = parameters.p_in / (gas.R * parameters.T_in)
    Q[:,:,1] = Q[:,:,0] * parameters.M_in * np.sqrt(gas.gamma*parameters.p_in/Q[:,:,0])
    Q[:,:,2] = Q[:,:,0] * 0
    Q[:,:,3] = thermo.calc_rho_et(parameters.p_in, Q[:,:,0], Q[:,:,1]/Q[:,:,0], Q[:,:,2]/Q[:,:,0], gas.gamma)

    class state:
        pass

    state.u = Q[:,:,1]/Q[:,:,0]
    state.v = Q[:,:,2]/Q[:,:,0]
    state.Q = Q
    state.Qn = Q
    state.p = thermo.calc_p( Q[:,:,0], Q[:,:,3], Q[:,:,1]/Q[:,:,0], Q[:,:,2]/Q[:,:,0], gas.gamma )
    state.T = state.p / (gas.R * Q[:,:,0])

    # boundary conditions

    from python.boundary.boundary_cond import enforce_bc, covariant
    state = enforce_bc(domain, mesh, parameters, state, gas)

    state = covariant(mesh, state)

    return state