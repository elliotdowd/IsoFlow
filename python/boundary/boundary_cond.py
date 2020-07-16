import numpy as np
from python.finite_volume.helper import thermo

# set boundary conditions (inviscid walls at bottom, outlet at top and right)

def enforce_bc_old(domain, mesh, parameters, state, gas):

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

    elif domain.name == 'NACA XXXX Airfoil' or domain.name == 'Biconvex Airfoil':

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
                                                                 gas, obj_i, obj_f, wallU, False)
            state.Q[obj_i:obj_f, wallL-1:wallL+1, :] = invisc_wall(state.Q[obj_i:obj_f, wallL-1:wallL+1, :], state.p[obj_i:obj_f, wallL], \
                                                                 state.T[obj_i:obj_f, wallL], mesh.s_proj[obj_i:obj_f, wallL-1:wallL+1, :], \
                                                                 gas, obj_i, obj_f, wallL, True)
        elif parameters.botwall == 'Viscous Wall':
            state.Q[obj_i:obj_f, wallU:wallU+2, :]  =  visc_wall(state.Q[obj_i:obj_f, wallU:wallU+2, :], state.p[obj_i:obj_f, wallU], \
                                                                 state.T[obj_i:obj_f, wallU], mesh.s_proj[obj_i:obj_f, wallU:wallU+2, :], \
                                                                 gas, obj_i-1, obj_f, wallU, False)
            state.Q[obj_i:obj_f, wallL-1:wallL+1, :] = visc_wall(state.Q[obj_i:obj_f, wallL-1:wallL+1, :], state.p[obj_i:obj_f, wallL], \
                                                                 state.T[obj_i:obj_f, wallL], mesh.s_proj[obj_i:obj_f, wallL-1:wallL+1, :], \
                                                                 gas, obj_i-1, obj_f, wallL, True)
            #state.Q[obj_i:obj_f, wallL-1:wallL+1, :] = np.array( (state.Q[obj_i:obj_f, wallL, :], state.Q[obj_i:obj_f, wallL-1, :]) ).reshape([obj_f-obj_i,2,4])

        # enforce inlet condition
        state.Q[0,:,0] = parameters.p_in / (gas.R_fn(gas.Cp[0,:], gas.Cv[0,:]) * parameters.T_in)
        state.Q[0,:,1] = state.Q[0,:,0] * parameters.M_in * np.cos(domain.alpha) * np.sqrt(gas.gamma_fn(gas.Cp[0,:], gas.Cv[0,:])*parameters.p_in/state.Q[0,:,0])
        state.Q[0,:,2] = state.Q[0,:,0] * parameters.M_in * np.sin(domain.alpha) * np.sqrt(gas.gamma_fn(gas.Cp[0,:], gas.Cv[0,:])*parameters.p_in/state.Q[0,:,0])
        state.Q[0,:,3] = thermo.calc_rho_et(parameters.p_in, state.Q[0,:,0], state.Q[0,:,1]/state.Q[0,:,0], state.Q[0,:,2]/state.Q[0,:,0], gas.gamma_fn(gas.Cp[0,:], gas.Cv[0,:]))                       

    return state


blend = lambda r: np.minimum( np.sqrt((1-r**2)**2)/(1+r**2), 1 )


def enforce_bc(domain, mesh, boundary, parameters, state, gas):

    for obj in boundary:
        n = obj.wall_n
        x0 = obj.wall_x
        x1 = obj.wall_x + n[0]
        y0 = obj.wall_y
        y1 = obj.wall_y + n[1]

        if obj.type == 'Outflow':

            # dynamic outflow boundary condition, leave alone if subsonic, propagate if supersonic
            state.Q[x0,y0,0] = parameters.p_in / (gas.R_fn(gas.Cp[x0,y0], gas.Cv[x0,y0]) * parameters.T_in)
            state.Q[x0,y0,1] = state.Q[x0,y0,0] * parameters.M_in * np.cos(domain.alpha) * np.sqrt(gas.gamma_fn(gas.Cp[x0,y0], gas.Cv[x0,y0])*parameters.p_in/state.Q[x0,y0,0])
            state.Q[x0,y0,2] = state.Q[x0,y0,0] * parameters.M_in * np.sin(domain.alpha) * np.sqrt(gas.gamma_fn(gas.Cp[x0,y0], gas.Cv[x0,y0])*parameters.p_in/state.Q[x0,y0,0])
            state.Q[x0,y0,3] = thermo.calc_rho_et(parameters.p_in, state.Q[x0,y0,0], state.u[x0,y0], state.v[x0,y0], gas.gamma_fn(gas.Cp[x0,y0], gas.Cv[x0,y0]))

            M = np.sqrt( (state.Q[x1,y1,1]/state.Q[x1,y1,0])**2 + (state.Q[x1,y1,2]/state.Q[x1,y1,0])**2 ) / \
                          thermo.calc_c( state.p[x1,y1], state.Q[x1,y1,0], gas.gamma_fn(gas.Cp[x1,y1], gas.Cv[x1,y1]) )
            # mask for state vector positions where M > 1
            sonic_mask = np.where( M > 1 )

            # mask for left/right vs. top or bottom outlets
            if n[0] == 0:
                state.Q[sonic_mask,y0,:] = state.Qn[sonic_mask,y1,:]
                state.p[sonic_mask,y0] = state.p[sonic_mask,y1]
                state.T[sonic_mask,y0] = state.T[sonic_mask,y1]
            elif n[0] != 1:
                pass
                # state.Q[x0,sonic_mask,:] = state.Qn[x0,sonic_mask,:]
                # state.p[x0,sonic_mask] = state.p[x1,sonic_mask]
                # state.T[x0,sonic_mask] = state.T[x1,sonic_mask]

        else:
            if y0 > y1:
                y = np.array( (y1, y0) )
            else: 
                y = np.array( (y0, y1) )

            # zero pressure gradient at wall
            state.p[x0,y0] = state.p[x1,y1]

            # thermal condition
            if obj.thermal == 'Adiabatic':
                state.T[x0,y0] = state.T[x1,y1]
            elif obj.thermal == 'Isothermal':
                state.T[x0,y0] = 2 * obj.Tw - state.T[x1,y1]
            elif obj.thermal == 'Fixed Temperature':
                state.T[x0,y0] = obj.Tw

            flip = n[1] == -1
            if obj.type == 'Inviscid Wall':
                state.Q[x0, y[0]:y[1]+1, :] = invisc_wall(state.Q[x0, y[0]:y[1]+1, :], state.p[x0, y0], state.T[x0, y0], mesh.s_proj[x0, y[0]:y[1]+1, :], \
                                                          gas, x0, y0, flip)
            elif obj.type == 'Viscous Wall':
                state.Q[x0, y[0]:y[1]+1, :] = visc_wall(state.Q[x0, y[0]:y[1]+1, :], state.p[x0, y0], state.T[x0, y0], mesh.s_proj[x0, y[0]:y[1]+1, :], \
                                                        gas, x0, y0, flip)

    return state


# compute velocities at inviscid slip wall, input Qwall[M+2, 2, 4]
def invisc_wall(Qwall, pwall, Twall, s_proj, gas, M, N, flip):

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
    boundary.slip(u0, v0, u1, v1, s_proj, len(M))

    Qwall[:, i, 0] = pwall / (gas.R_fn(gas.Cp[M,N], gas.Cv[M, N]) * Twall)
    Qwall[:, i, 1] = u0 * Qwall[:, i, 0]
    Qwall[:, i, 2] = v0 * Qwall[:, i, 0]
    Qwall[:, i, 3] = thermo.calc_rho_et( pwall, Qwall[:, i, 0], Qwall[:, i, 1]/Qwall[:, i, 0], \
                                         Qwall[:, i, 2]/Qwall[:, i, 0], gas.gamma_fn(gas.Cp[M,N], gas.Cv[M,N]) )
    
    return Qwall


# compute velocities at viscous no-slip wall, input Qwall[M+2, 2, 4]
def visc_wall(Qwall, pwall, Twall, s_proj, gas, M, N, flip):

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

    Qwall[:, i, 0] = pwall / (gas.R_fn(gas.Cp[M,N], gas.Cv[M,N]) * Twall)
    Qwall[:, i, 1] = u0 * Qwall[:, i, 0]
    Qwall[:, i, 2] = v0 * Qwall[:, i, 0]
    Qwall[:, i, 3] = thermo.calc_rho_et( pwall, Qwall[:, i, 0], Qwall[:, i, 1]/Qwall[:, i, 0], \
                                         Qwall[:, i, 2]/Qwall[:, i, 0], gas.gamma_fn(gas.Cp[M,N], gas.Cv[M,N]) )
    
    return Qwall


# covariant velocities 

def covariant(mesh, state):

    state.U = (1/mesh.s_proj[:,:,4]) * \
              (state.u*mesh.s_proj[:,:,0] + state.v*mesh.s_proj[:,:,1])
       
    state.V = (1/mesh.s_proj[:,:,5]) * \
              (state.u*mesh.s_proj[:,:,2] + state.v*mesh.s_proj[:,:,3])

    return state