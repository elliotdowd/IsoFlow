def TicTocGenerator():
    # Generator that returns time differences
    import time
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "init_state time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)


# initialize Q vector given simulation parameters
def init_state(domain, mesh, parameters, gas):

    import numpy as np
    from helper import thermo

    tic()

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

    from boundary_cond import enforce_bc, covariant
    state = enforce_bc(domain, mesh, state, gas)

    state = covariant(mesh, state)

    toc()

    return state