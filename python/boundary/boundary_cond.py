import numpy as np
from python.finite_volume.helper import thermo

# set boundary conditions (inviscid walls at bottom, outlet at top and right)

blend = lambda r: np.minimum( np.sqrt((1-r**2)**2)/(1+r**2), 1 )

def enforce_bc(domain, mesh, boundary, parameters, state, gas):

    for obj in boundary:
        n = obj.wall_n
        x0 = obj.wall_x
        x1 = obj.wall_x + n[0]
        y0 = obj.wall_y
        y1 = obj.wall_y + n[1]

        if obj.type == 'Outflow':

            if parameters.outlet == 'supersonic':
                if n[0] == 0:
                    state.Q[x0,y0,:] = state.Qn[x1,y1,:]
                    state.p[x0,y0] = state.p[x1,y1]
                    state.T[x0,y0] = state.T[x1,y1]
                else:
                    state.Q[x0,y0,:] = state.Qn[x1,y1,:]
                    state.p[x0,y0] = state.p[x1,y1]
                    state.T[x0,y0] = state.T[x1,y1]

            elif parameters.outlet == 'hybrid':
                # dynamic outflow boundary condition, leave alone if subsonic, propagate if supersonic
                state.Q[x0,y0,0] = parameters.p_out / (gas.R_fn(gas.Cp[x0,y0], gas.Cv[x0,y0]) * parameters.T_out)
                state.Q[x0,y0,1] = state.Q[x0,y0,0] * parameters.M_out * np.cos(domain.alpha) * np.sqrt(gas.gamma_fn(gas.Cp[x0,y0], gas.Cv[x0,y0])*parameters.p_out/state.Q[x0,y0,0])
                state.Q[x0,y0,2] = state.Q[x0,y0,0] * parameters.M_out * np.sin(domain.alpha) * np.sqrt(gas.gamma_fn(gas.Cp[x0,y0], gas.Cv[x0,y0])*parameters.p_out/state.Q[x0,y0,0])
                state.Q[x0,y0,3] = thermo.calc_rho_et(parameters.p_out, state.Q[x0,y0,0], state.u[x0,y0], state.v[x0,y0], gas.gamma_fn(gas.Cp[x0,y0], gas.Cv[x0,y0]))

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
                    state.Q[x0,sonic_mask,:] = state.Qn[x0,sonic_mask,:]
                    state.p[x0,sonic_mask] = state.p[x1,sonic_mask]
                    state.T[x0,sonic_mask] = state.T[x1,sonic_mask]

            elif parameters.outlet == 'subsonic':
                state.Q[x0,y0,0] = parameters.p_out / (gas.R_fn(gas.Cp[x0,y0], gas.Cv[x0,y0]) * parameters.T_out)
                state.Q[x0,y0,1] = state.Q[x0,y0,0] * parameters.M_out * np.cos(domain.alpha) * np.sqrt(gas.gamma_fn(gas.Cp[x0,y0], gas.Cv[x0,y0])*parameters.p_out/state.Q[x0,y0,0])
                state.Q[x0,y0,2] = state.Q[x0,y0,0] * parameters.M_out * np.sin(domain.alpha) * np.sqrt(gas.gamma_fn(gas.Cp[x0,y0], gas.Cv[x0,y0])*parameters.p_out/state.Q[x0,y0,0])
                state.Q[x0,y0,3] = thermo.calc_rho_et(parameters.p_out, state.Q[x0,y0,0], state.u[x0,y0], state.v[x0,y0], gas.gamma_fn(gas.Cp[x0,y0], gas.Cv[x0,y0]))

        elif obj.type == 'Inlet':
            state.Q[x0,y0,0] = parameters.p_in / (gas.R_fn(gas.Cp[x0,y0], gas.Cv[x0,y0]) * parameters.T_in)
            state.Q[x0,y0,1] = state.Q[x0,y0,0] * parameters.M_in * np.cos(domain.alpha) * np.sqrt(gas.gamma_fn(gas.Cp[x0,y0], gas.Cv[x0,y0])*parameters.p_in/state.Q[x0,y0,0])
            state.Q[x0,y0,2] = state.Q[x0,y0,0] * parameters.M_in * np.sin(domain.alpha) * np.sqrt(gas.gamma_fn(gas.Cp[x0,y0], gas.Cv[x0,y0])*parameters.p_in/state.Q[x0,y0,0])
            state.Q[x0,y0,3] = thermo.calc_rho_et(parameters.p_in, state.Q[x0,y0,0], state.u[x0,y0], state.v[x0,y0], gas.gamma_fn(gas.Cp[x0,y0], gas.Cv[x0,y0]))

        # elif obj.type == 'Symmetry':
        #     if len(y0) > 1:
                

        #     else:
        #         x = 
        #         if y0 > y1:
        #             y = np.array( (y1, y0) )
        #         else: 
        #             y = np.array( (y0, y1) )

        #     state.T[x0,y0] = state.T[x1,y1]
        #     state.p[x0,y0] = state.p[x1,y1]

        #     if obj.wall_n[0] == 0:
        #         pass
            # state.Q[x0, y[0]:y[1]+1, :] = invisc_wall(state.Q[x0, y[0]:y[1]+1, :], state.p[x0, y0], state.T[x0, y0], mesh.s_proj[x0, y[0]:y[1]+1, :], \
            #                                             gas, x0, y0, flip)
            # state.Q[x0,y0] = 

            # state.Q[x0, y[0]:y[1]+1, :] = invisc_wall(state.Q[x0, y[0]:y[1]+1, :], state.p[x0, y0], state.T[x0, y0], mesh.s_proj[x0, y[0]:y[1]+1, :], \
            #                                             gas, x0, y0, flip)

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