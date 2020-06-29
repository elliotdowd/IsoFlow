

# set boundary conditions (inviscid walls at bottom, outlet at top and right)

def enforce_bc(domain, mesh, parameters, state, gas):

    import numpy as np
    from finite_volume.helper import thermo

    # enforce inlet condition
    state.Q[0,:,0] = parameters.p_in / (gas.R * parameters.T_in)
    state.Q[0,:,1] = state.Q[0,:,0] * parameters.M_in * np.sqrt(gas.gamma*parameters.p_in/state.Q[0,:,0])
    state.Q[0,:,2] = state.Q[0,:,0] * 0
    state.Q[0,:,3] = thermo.calc_rho_et(parameters.p_in, state.Q[0,:,0], state.Q[0,:,1]/state.Q[0,:,0], state.Q[0,:,2]/state.Q[0,:,0], gas.gamma)

    wedge_i = abs(mesh.yy[:,1]) > 0.0000001
    wedge_i = wedge_i[0]

    state.p[:, 0] = state.p[:, 1]
    state.T[0:wedge_i-1,0] = state.T[0:wedge_i-1,1]
    state.T[wedge_i:, 0] = 300
    state.Q[:, 0:2, :] = invisc_wall(state.Q[:, 0:2, :], state.p[:, 0], state.T[:, 0], mesh.s_proj[:, 0:2, :], domain.M+2, gas)


    # topwall = invisc_wall(np.array([state.Q[:, domain.N+1, :], state.Q[:, domain.N, :]]).reshape([domain.M+2, 2, 4]), state.p[:, -1], state.T[:, -1], \
    #                       np.array([mesh.s_proj[:, domain.N+1, :], mesh.s_proj[:, domain.N, :]]).reshape([domain.M+2, 2, 6]), domain.M+2, gas)
    # state.Q[:, domain.N:, :] = np.array([topwall[:,1,:], topwall[:,0,:]]).reshape([domain.M+2, 2, 4])
    
    state.Q[-1, :, :] = state.Qn[-1, :, :]
    # state.p[domain.M+1, :] = state.p[domain.M, :]
    # state.T[domain.M+1, :] = state.T[domain.M, :]

    state.Q[:, -1, :] = state.Qn[:, -1, :]
    #state.p[:, domain.N+1] = state.p[:, domain.N]
    #state.T[:, domain.N+1] = state.T[:, domain.N]

    return state


# compute velocities at inviscid slip wall, input Qwall[M+2, 2, 4]
def invisc_wall(Qwall, pwall, Twall, s_proj, M, gas):

    import numpy as np

    from finite_volume.helper import thermo
    import boundary.boundary as boundary

    u1 = Qwall[:, 1, 1] / Qwall[:, 1, 0]
    v1 = Qwall[:, 1, 2] / Qwall[:, 1, 0]

    # slip wall boundary condition

    u0 = Qwall[:, 0, 1] / Qwall[:, 0, 0]
    v0 = Qwall[:, 0, 2] / Qwall[:, 0, 0]

    # run Fortran 90 subroutine to determine wall velocities
    boundary.slip(u0, v0, u1, v1, s_proj, M)

    Qwall[:, 0, 0] = pwall / (gas.R * Twall)
    Qwall[:, 0, 1] = u0 * Qwall[:, 0, 0]
    Qwall[:, 0, 2] = v0 * Qwall[:, 0, 0]
    Qwall[:, 0, 3] = thermo.calc_rho_et( pwall, Qwall[:, 0, 0], Qwall[:, 0, 1]/Qwall[:, 0, 0], Qwall[:, 0, 2]/Qwall[:, 0, 0], gas.gamma )
    
    return Qwall


# compute velocities at viscous no-slip wall, input Qwall[M+2, 2, 4]
def visc_wall(Qwall, pwall, Twall, s_proj, M, gas):

    import numpy as np

    from helper import thermo
    import boundary

    u1 = Qwall[:, 1, 1] / Qwall[:, 1, 0]
    v1 = Qwall[:, 1, 2] / Qwall[:, 1, 0]

    # no-slip wall condition

    u0 = -u1
    v0 = -v1

    Qwall[:, 0, 0] = pwall / (gas.R * Twall)
    Qwall[:, 0, 1] = u0 * Qwall[:, 0, 0]
    Qwall[:, 0, 2] = v0 * Qwall[:, 0, 0]
    Qwall[:, 0, 3] = thermo.calc_rho_et( pwall, Qwall[:, 0, 0], Qwall[:, 0, 1]/Qwall[:, 0, 0], Qwall[:, 0, 2]/Qwall[:, 0, 0], gas.gamma )
    
    return Qwall


# covariant velocities 

def covariant(mesh, state):

    state.U = (1/mesh.s_proj[:,:,4]) * \
              (state.u*mesh.s_proj[:,:,0] + state.v*mesh.s_proj[:,:,1])
       
    state.V = (1/mesh.s_proj[:,:,5]) * \
              (state.u*mesh.s_proj[:,:,2] + state.v*mesh.s_proj[:,:,3])

    return state