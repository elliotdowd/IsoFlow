

def init_state(domain, mesh, parameters, gas):

    import numpy as np

    # define anonymous functions
    rho_et = lambda p, rho, u, v, gam: (p/(gam-1)) + 0.5*rho*(u**2 + v**2)

    # initialize Q state vector at each point
    Q = np.zeros((domain.M+2, domain.N+2, 4), dtype='float', order='F')
    #soln_vars.init_q(Q, parameters.p_in, parameters.T_in, parameters.M_in, gas.R, gas.gamma, domain.M+2, domain.M+2)

    Q[:,:,0] = parameters.p_in / (gas.R * parameters.T_in)
    Q[:,:,1] = Q[:,:,0] * parameters.M_in * np.sqrt(gas.gamma*parameters.p_in/Q[:,:,0])
    Q[:,:,2] = Q[:,:,0] * 0
    Q[:,:,3] = rho_et(parameters.p_in, Q[:,:,0], Q[:,:,1]/Q[:,:,0], Q[:,:,2]/Q[:,:,0], gas.gamma)


    class state:
        pass

    state.Q = Q
    #state.U = U
    #state.V = V

    return state