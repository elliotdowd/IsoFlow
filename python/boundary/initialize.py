

# initialize Q vector given simulation parameters
def init_state(domain, mesh, parameters, gas):

    import numpy as np
    from python.finite_volume.helper import thermo

    # initialize Q state vector at each point
    Q = np.zeros((domain.M+2, domain.N+2, 4), dtype='float', order='F')
    #soln_vars.init_q(Q, parameters.p_in, parameters.T_in, parameters.M_in, gas.R, gas.gamma, domain.M+2, domain.M+2)

    Q[:,:,0] = parameters.p_in / (gas.R_p * parameters.T_in)
    Q[:,:,1] = Q[:,:,0] * parameters.M_in * np.cos(domain.alpha) * np.sqrt(gas.gamma_p*parameters.p_in/Q[:,:,0])
    Q[:,:,2] = Q[:,:,0] * parameters.M_in * np.sin(domain.alpha) * np.sqrt(gas.gamma_p*parameters.p_in/Q[:,:,0])
    Q[:,:,3] = thermo.calc_rho_et(parameters.p_in, Q[:,:,0], Q[:,:,1]/Q[:,:,0], Q[:,:,2]/Q[:,:,0], gas.gamma_p)

    class state:
        pass

    state.u = Q[:,:,1]/Q[:,:,0]
    state.v = Q[:,:,2]/Q[:,:,0]
    state.Q = Q
    state.Qn = Q
    state.p = thermo.calc_p( Q[:,:,0], Q[:,:,3], Q[:,:,1]/Q[:,:,0], Q[:,:,2]/Q[:,:,0], gas.gamma_p )
    state.T = state.p / (gas.R_p * Q[:,:,0])
    state.Mach = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 ) / thermo.calc_c( state.p, state.Q[:,:,0], gas.gamma_p )
    state.vel = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 )
    state.p0 = (1+((gas.gamma_p-1)/2)*state.Mach**2)**(gas.gamma_p/(gas.gamma_p-1)) * state.p
    state.res = np.array([1, 1, 1, 1])
    state.T0 = (1+((gas.gamma_p-1)/2)*state.Mach**2) * state.T
    state.n = 1

    # specific heat ratio calculation
    gas.Cp = gas.Cp_fn( gas.gamma_p, gas.Cp_p, gas.theta, state.T )
    gas.Cv = gas.Cv_fn( gas.gamma_p, gas.Cv_p, gas.theta, state.T )

    # boundary conditions
    from python.boundary.boundary_cond import enforce_bc, covariant
    state = enforce_bc(domain, mesh, parameters, state, gas)

    state = covariant(mesh, state)

    return state