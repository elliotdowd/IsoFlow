

# set boundary conditions (inviscid walls at bottom, outlet at top and right)

def enforce_bc(domain, mesh, parameters, state, gas):

    import numpy as np
    from python.finite_volume.helper import thermo

    gas.Cp = gas.Cp_fn( gas.gamma_p, gas.Cp_p, gas.theta, state.T )
    gas.Cv = gas.Cv_fn( gas.gamma_p, gas.Cv_p, gas.theta, state.T )

    if domain.name == 'Wedge' or domain.name == 'Corner':

        if domain.name == 'Wedge':
            obj_i = mesh.xx[:,1] > domain.obj_start
            obj_i = np.where(obj_i)
            obj_i = obj_i[0][0]
            obj_f = domain.M+2
        elif domain.name == 'Corner':
            obj_i = mesh.xx[:,1] > domain.obj_start
            obj_i = np.where(obj_i) 
            obj_i = obj_i[0][0]
            obj_f = mesh.xx[:,1] > domain.obj_end
            if np.any(obj_f):
                obj_f = np.where(obj_f)
                obj_f = obj_f[0][0]
            else:
                obj_f = domain.M+2

        # update state variables at bottom wall
        state.p[:, 0] = state.p[:, 1]
        state.T[:, 0] = state.T[:,0]

        if parameters.botwall_thermal == 'Adiabatic':
            state.T[:,0] = state.T[:,1]
        elif parameters.botwall_thermal == 'Isothermal':
            state.T[obj_i:obj_f,0] = 2 * parameters.botwall_temp - state.T[obj_i:obj_f,1]
        elif parameters.botwall_thermal == 'Fixed Temperature':
            state.T[obj_i:obj_f,0] = parameters.botwall_temp

        if parameters.botwall == 'Inviscid Wall':
            state.Q[:, 0:2, :] = invisc_wall(state.Q[:, 0:2, :], state.p[:, 0], state.T[:, 0], mesh.s_proj[:, 0:2, :], \
                                            gas, 0, domain.M+2, domain.N+1, False)
        elif parameters.botwall == 'Viscous Wall':
            state.Q[obj_i:obj_f, 0:2, :] = visc_wall(state.Q[obj_i:obj_f, 0:2, :], \
                    state.p[obj_i:obj_f, 0], state.T[obj_i:obj_f, 0], mesh.s_proj[obj_i:obj_f, 0:2, :], \
                                                    gas, obj_i-1, obj_f, 0, False)

        if parameters.topwall_thermal == 'Adiabatic':
            state.T[:,domain.N+1] = state.T[:,domain.N]
        elif parameters.topwall_thermal == 'Isothermal':
            state.T[:,domain.N+1] = 2 * parameters.topwall_temp - state.T[:,domain.N]
        elif parameters.topwall_thermal == 'Fixed Temperature':
            state.T[:,domain.N+1] = parameters.topwall_temp

        if parameters.topwall == 'Outflow':
            state.Q[:, domain.N+1, :] = state.Qn[:, domain.N+1, :]
        elif parameters.topwall == 'Inviscid Wall':
            # update state variables at top wall
            state.p[:, domain.N+1] = state.p[:, domain.N]
            state.T[:,domain.N+1] = state.T[:,domain.N]
            state.T[:, domain.N+1] = 300

            state.Q[:, domain.N:, :] = invisc_wall(state.Q[:, domain.N:, :], state.p[:, domain.N+1], state.T[:, domain.N+1], \
                                                mesh.s_proj[:, domain.N:, :], gas, 0, domain.M+2, domain.N+1, True)
        elif parameters.topwall == 'Viscous Wall':
            # update state variables at top wall
            state.p[:, domain.N+1] = state.p[:, domain.N]
            state.T[:,domain.N+1] = state.T[:,domain.N]
            state.T[:, domain.N+1] = 300

            state.Q[:, domain.N:, :] = visc_wall(state.Q[:, domain.N:, :], \
                                    state.p[:, domain.N+1], state.T[:, domain.N+1], mesh.s_proj[:, domain.N:, :], gas, -1, domain.M+2, domain.N+1, True)

        # enforce inlet condition
        state.Q[0,:,0] = parameters.p_in / (gas.R_fn(gas.Cp[0,:], gas.Cv[0,:]) * parameters.T_in)
        state.Q[0,:,1] = state.Q[0,:,0] * parameters.M_in * np.sqrt(gas.gamma_fn(gas.Cp[0,:], gas.Cv[0,:])*parameters.p_in/state.Q[0,:,0])
        state.Q[0,:,2] = state.Q[0,:,0] * 0
        state.Q[0,:,3] = thermo.calc_rho_et(parameters.p_in, state.Q[0,:,0], state.Q[0,:,1]/state.Q[0,:,0], state.Q[0,:,2]/state.Q[0,:,0], gas.gamma_fn(gas.Cp[0,:], gas.Cv[0,:]))

        # right side outflow
        state.Q[domain.M+1, :, :] = state.Qn[domain.M, :, :]

    elif domain.name == 'Cylinder':

        obj_i = 0
        obj_f = domain.M+2

        # update state variables at cylinder  wall
        state.p[:, 0] = state.p[:, 1]
        state.T[:, 0] = state.T[:,0]

        if parameters.botwall_thermal == 'Adiabatic':
            state.T[:,0] = state.T[:,1]
        elif parameters.botwall_thermal == 'Isothermal':
            state.T[:,0] = 2 * parameters.botwall_temp - state.T[:,1]
        elif parameters.botwall_thermal == 'Fixed Temperature':
            state.T[:,0] = parameters.botwall_temp

        if parameters.botwall == 'Inviscid Wall':
            state.Q[:, 0:2, :] = invisc_wall(state.Q[:, 0:2, :], state.p[:, 0], state.T[:, 0], mesh.s_proj[:, 0:2, :], \
                                            gas, 0, domain.M+2, 0, False)
        elif parameters.botwall == 'Viscous Wall':
            state.Q[:, 0:2, :] = visc_wall(state.Q[:, 0:2, :], \
                    state.p[:, 0], state.T[:, 0], mesh.s_proj[:, 0:2, :], \
                                                    gas, -1, domain.M+2, 0, False)

        # enforce inlet condition at front half of outer circle
        half = int(domain.M/2)
        state.Q[0:half,-1,0] = parameters.p_in / (gas.R_fn(gas.Cp[0:half,-1], gas.Cv[0:half,-1]) * parameters.T_in)
        state.Q[0:half,-1,1] = state.Q[0:half,-1,0] * parameters.M_in * np.sqrt(gas.gamma_fn(gas.Cp[0:half,-1], gas.Cv[0:half,-1])*parameters.p_in/state.Q[0:half,-1,0])
        state.Q[0:half,-1,2] = state.Q[0:half,-1,0] * 0
        state.Q[0:half,-1,3] = thermo.calc_rho_et(parameters.p_in, state.Q[0:half,-1,0], state.Q[0:half,-1,1]/state.Q[0:half,-1,0], \
                                                  state.Q[0:half,-1,2]/state.Q[0:half,-1,0], gas.gamma_fn(gas.Cp[0:half,-1], gas.Cv[0:half,-1]))

        state.Q[half+1:-1,-1,:] = state.Qn[half+1:-1,-2,:]

        # symmetry conditions
        state.p[:,0] = state.p[:,-1]
        state.T[:,0] = state.T[:,-1]
        state.Q[:,0,:] = state.Q[:,-1,:]

    elif domain.name == 'NACA XXXX Airfoil':

        obj_i = domain.obj_i
        obj_f = domain.obj_f
        wallL = domain.wallL
        wallU = domain.wallU-1

        state.p[obj_i:obj_f,wallU] = state.p[obj_i:obj_f,wallU+1]
        state.p[obj_i:obj_f,wallL] = state.p[obj_i:obj_f,wallL-1]

        if parameters.botwall_thermal == 'Adiabatic':
            state.T[obj_i:obj_f,wallU] = state.T[obj_i:obj_f,wallU+1]
            state.T[obj_i:obj_f,wallL] = state.T[obj_i:obj_f,wallL-1]
        elif parameters.botwall_thermal == 'Isothermal':
            state.T[obj_i:obj_f,wallL] = 2 * parameters.botwall_temp - state.T[obj_i:obj_f,wallL-1]
            state.T[obj_i:obj_f,wallU] = 2 * parameters.botwall_temp - state.T[obj_i:obj_f,wallU+1]
        elif parameters.botwall_thermal == 'Fixed Temperature':
            state.T[obj_i:obj_f,wallL] = parameters.botwall_temp
            state.T[obj_i:obj_f,wallU] = parameters.botwall_temp

        if parameters.botwall == 'Inviscid Wall':
            state.Q[obj_i:obj_f, wallU:wallU+2, :] = invisc_wall(state.Q[obj_i:obj_f, wallU:wallU+2, :], state.p[obj_i:obj_f, wallU], \
                                                                 state.T[obj_i:obj_f, wallU], mesh.s_proj[obj_i:obj_f, wallU:wallU+2, :], \
                                                                 gas, obj_i, obj_f, wallL, False)
            state.Q[obj_i:obj_f, wallL-1:wallL+1, :] = invisc_wall(state.Q[obj_i:obj_f, wallL-1:wallL+1, :], state.p[obj_i:obj_f, wallL], \
                                                                 state.T[obj_i:obj_f, wallL], mesh.s_proj[obj_i:obj_f, wallL-1:wallL+1, :], \
                                                                 gas, obj_i, obj_f, wallL, True)
        elif parameters.botwall == 'Viscous Wall':
            state.Q[obj_i:obj_f, wallU:wallU+2, :] = visc_wall(state.Q[obj_i:obj_f, wallU:wallU+2, :], \
                    state.p[obj_i:obj_f, wallU], state.T[obj_i:obj_f, wallU], mesh.s_proj[obj_i:obj_f, wallU:wallU+2, :], \
                                                    gas, obj_i-1, obj_f, wallU, False)
            state.Q[obj_i:obj_f, wallL-1:wallL+1, :] = visc_wall(state.Q[obj_i:obj_f, wallL-1:wallL+1, :], \
                    state.p[obj_i:obj_f, wallL], state.T[obj_i:obj_f, wallL], mesh.s_proj[obj_i:obj_f, wallL-1:wallL+1, :], \
                                                    gas, obj_i-1, obj_f, wallL, True)

        # enforce inlet condition
        state.Q[0,:,0] = parameters.p_in / (gas.R_fn(gas.Cp[0,:], gas.Cv[0,:]) * parameters.T_in)
        state.Q[0,:,1] = state.Q[0,:,0] * parameters.M_in * np.sqrt(gas.gamma_fn(gas.Cp[0,:], gas.Cv[0,:])*parameters.p_in/state.Q[0,:,0])
        state.Q[0,:,2] = state.Q[0,:,0] * 0
        state.Q[0,:,3] = thermo.calc_rho_et(parameters.p_in, state.Q[0,:,0], state.Q[0,:,1]/state.Q[0,:,0], state.Q[0,:,2]/state.Q[0,:,0], gas.gamma_fn(gas.Cp[0,:], gas.Cv[0,:]))                       

    return state


# compute velocities at inviscid slip wall, input Qwall[M+2, 2, 4]
def invisc_wall(Qwall, pwall, Twall, s_proj, gas, obj_i, obj_f, N, flip):

    import numpy as np

    from python.finite_volume.helper import thermo
    import python.boundary.boundary as boundary

    if flip:
        i = 1
        j = 0
    else:
        i = 0
        j = 1

    u1 = Qwall[:, j, 1] / Qwall[:, j, 0]
    v1 = Qwall[:, j, 2] / Qwall[:, j, 0]

    # slip wall boundary condition

    u0 = Qwall[:, i, 1] / Qwall[:, i, 0]
    v0 = Qwall[:, i, 2] / Qwall[:, i, 0]

    # run Fortran 90 subroutine to determine wall velocities
    boundary.slip(u0, v0, u1, v1, s_proj, obj_f-obj_i)

    Qwall[:, i, 0] = pwall / (gas.R_fn(gas.Cp[obj_i:obj_f,N], gas.Cv[obj_i:obj_f, N]) * Twall)
    Qwall[:, i, 1] = u0 * Qwall[:, i, 0]
    Qwall[:, i, 2] = v0 * Qwall[:, i, 0]
    Qwall[:, i, 3] = thermo.calc_rho_et( pwall, Qwall[:, i, 0], Qwall[:, i, 1]/Qwall[:, i, 0], \
                                         Qwall[:, i, 2]/Qwall[:, i, 0], gas.gamma_fn(gas.Cp[obj_i:obj_f,N], gas.Cv[obj_i:obj_f,N]) )
    
    return Qwall


# compute velocities at viscous no-slip wall, input Qwall[M+2, 2, 4]
def visc_wall(Qwall, pwall, Twall, s_proj, gas, obj_i, obj_f, N, flip):

    import numpy as np
    from python.finite_volume.helper import thermo

    if flip:
        i = 1
        j = 0
    else:
        i = 0
        j = 1

    u1 = Qwall[:, j, 1] / Qwall[:, j, 0]
    v1 = Qwall[:, j, 2] / Qwall[:, j, 0]

    # no-slip wall condition

    u0 = -u1
    v0 = -v1

    Qwall[:, i, 0] = pwall / (gas.R_fn(gas.Cp[obj_i+1:obj_f,N], gas.Cv[obj_i+1:obj_f,N]) * Twall)
    Qwall[:, i, 1] = u0 * Qwall[:, i, 0]
    Qwall[:, i, 2] = v0 * Qwall[:, i, 0]
    Qwall[:, i, 3] = thermo.calc_rho_et( pwall, Qwall[:, i, 0], Qwall[:, i, 1]/Qwall[:, i, 0], Qwall[:, i, 2]/Qwall[:, i, 0], \
                                         gas.gamma_fn(gas.Cp[obj_i+1:obj_f,N], gas.Cv[obj_i+1:obj_f,N]) )
    
    return Qwall


# covariant velocities 

def covariant(mesh, state):

    state.U = (1/mesh.s_proj[:,:,4]) * \
              (state.u*mesh.s_proj[:,:,0] + state.v*mesh.s_proj[:,:,1])
       
    state.V = (1/mesh.s_proj[:,:,5]) * \
              (state.u*mesh.s_proj[:,:,2] + state.v*mesh.s_proj[:,:,3])

    return state