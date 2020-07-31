# module importation
import numpy as np
from python.finite_volume.helper import thermo, split
from python.boundary.boundary_cond import enforce_bc, covariant
from python.finite_volume.timestepping import local_timestep
import python.finite_volume.soln_vars as soln_vars
import python.finite_volume.flux as flux
import python.finite_volume.muscl as muscl
from pytictoc import TicToc

# AUSM flux vector splitting scheme 
def AUSM( domain, mesh, boundary, parameters, state, gas ):

    print('AUSM Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # initialize speed of sound arrays
    c_half_zeta = np.zeros( (domain.M+1, domain.N+2), dtype='float', order='F' )
    c_half_eta = np.zeros( (domain.M+2, domain.N+1), dtype='float', order='F' )

    # initialize phi and p state vectors
    Phi = np.zeros( (domain.M+2, domain.N+2, 4) )
    Phi[:,:,0] = np.ones( (domain.M+2, domain.N+2) )

    P_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    P_zeta[:,:,0] = np.zeros( (domain.M+1, domain.N+2) )
    P_zeta[:,:,3] = np.zeros( (domain.M+1, domain.N+2) )

    P_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    P_eta[:,:,0] = np.zeros( (domain.M+2, domain.N+1) )
    P_eta[:,:,3] = np.zeros( (domain.M+2, domain.N+1) )

    E_hat_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    Ev_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    E_hat_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    Ev_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    Fv_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    Fv_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( [domain.M, domain.N, 4], dtype='float', order='F' )
    if state.n == 0:
        state.res = np.ones( [parameters.iterations + 1, 4] )
    else:
        state.res = np.vstack( (state.res, np.ones( [parameters.iterations + 1, 4] )) )

    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])

    t = TicToc()

    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        t.tic()

        n = n+1
        state.n = state.n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # local timestepping
        state = local_timestep( domain, mesh, state, parameters, gas )

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        # state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / \
        #                               state.Q[:,:,0] + state.p/state.Q[:,:,0]
        state.ht = gas.ht_fn( gas, state )
        # state.Q[:,:,3] = state.Q[:,:,0]*(state.ht - state.p/state.Q[:,:,0])

        # speed of sound at cell interfaces
        # from Liou 2006 (JCP 214)
        flux.c_entropy_fixed( c_half_zeta, c_half_eta, state.ht, gas.gamma_fn(gas.Cp, gas.Cv), \
                              state.U, state.V, domain.M, domain.N )

        # cell face Mach numbers
        M_L = state.U[0:-1,:] / c_half_zeta
        M_R = state.U[1:,:] / c_half_zeta
        M_D = state.V[:,0:-1] / c_half_eta
        M_U = state.V[:,1:] / c_half_eta

        # split interface Mach numbers in the zeta and eta directions
        M_half_zeta = split.Mvp( M_L ) + split.Mvm( M_R )
        M_half_eta = split.Mvp( M_D ) + split.Mvm( M_U )

        # calculate mass flux at cell interfaces
        mdot_half_zeta = c_half_zeta * M_half_zeta * ( np.double(M_half_zeta>0) * state.Q[0:-1,:,0] + np.double(M_half_zeta<=0) * state.Q[1:,:,0] )
        mdot_half_eta =  c_half_eta  * M_half_eta  * ( np.double(M_half_eta>0)  * state.Q[:,0:-1,0] + np.double(M_half_eta<=0)  * state.Q[:,1:,0] )

        cr = 1e-60
        # calculate pressure flux at cell interfaces
        p_half_zeta = split.P1p( M_L+cr )*state.p[0:-1,:] + split.P1m( M_R+cr )*state.p[1:,:]
        p_half_eta =  split.P1p( M_D+cr )*state.p[:,0:-1] + split.P1m( M_U+cr )*state.p[:,1:]

        # initialize Phi vector components
        Phi[:,:,1] = state.u
        Phi[:,:,2] = state.v
        Phi[:,:,3] = state.ht
        
        # pressure flux vector
        P_zeta[:,:,1] = p_half_zeta * mesh.s_proj[0:-1,:,0] / mesh.s_proj[0:-1,:,4]
        P_zeta[:,:,2] = p_half_zeta * mesh.s_proj[0:-1,:,1] / mesh.s_proj[0:-1,:,4]

        P_eta[:,:,1] = p_half_eta * mesh.s_proj[:,0:-1,2] / mesh.s_proj[:,0:-1,5]
        P_eta[:,:,2] = p_half_eta * mesh.s_proj[:,0:-1,3] / mesh.s_proj[:,0:-1,5]

        flux.face_flux( mdot_half_zeta, mdot_half_eta, Phi, P_zeta, P_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, domain.M, domain.N)

        if parameters.viscModel == 'viscous laminar':
            # calculate viscous stress terms
            Txx_m, Tyy_m, Txy_m = calc_stress( mesh, state, gas, '-' )
            Txx_p, Tyy_p, Txy_p = calc_stress( mesh, state, gas, '+' )

            Ev_left[:,:,1] = Txx_m
            Ev_left[:,:,2] = Txy_m
            Ev_left[:,:,3] = state.u[0:-2,1:-1]*Txx_m + state.v[0:-2,1:-1]*Txy_m
            Ev_right[:,:,1] = Txx_p
            Ev_right[:,:,2] = Txy_p
            Ev_right[:,:,3] = state.u[1:-1,1:-1]*Txx_p + state.v[1:-1,1:-1]*Txy_p

            Fv_bot[:,:,1] = Txy_m
            Fv_bot[:,:,2] = Tyy_m
            Fv_bot[:,:,3] = state.u[1:-1,0:-2]*Txy_m + state.v[1:-1,0:-2]*Tyy_m
            Fv_top[:,:,1] = Txy_p
            Fv_top[:,:,2] = Tyy_p
            Fv_top[:,:,3] = state.u[1:-1,1:-1]*Txy_p + state.v[1:-1,1:-1]*Tyy_p

            E_hat_left = E_hat_left - Ev_left
            E_hat_right = E_hat_right - Ev_right
            F_hat_bot = F_hat_bot - Fv_bot
            F_hat_top = F_hat_top - Fv_top

        elif parameters.viscModel == 'inviscid':
            pass

        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        # state.residual = -np.dstack( (state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1]) ) * \
        #                           ( ( E_hat_right * np.dstack([mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4]]) \
        #                             - E_hat_left * np.dstack([mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4]]) ) + \
        #                             ( F_hat_top * np.dstack([mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5]])\
        #                             - F_hat_bot * np.dstack([mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5]]) ) )
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj, domain.M, domain.N ) 
        state.Q[1:-1,1:-1,:] = state.Qn[1:-1,1:-1,:] + state.residual / mesh.dV4

        # L_inf-norm residual
        state.res[state.n-1,0] = np.log10( np.max(state.residual[:,:,0] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,1] = np.log10( np.max(state.residual[:,:,1] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,2] = np.log10( np.max(state.residual[:,:,2] * mesh.dV[1:-1,1:-1]) )
        state.res[state.n-1,3] = np.log10( np.max(state.residual[:,:,3] * mesh.dV[1:-1,1:-1]) ) 

        #state.res[n-1] = np.log10( np.max(state.residual * mesh.dV4) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], \
                                 state.Q[:,:,1]/state.Q[:,:,0], state.Q[:,:,2]/state.Q[:,:,0], gas.gamma_fn(gas.Cp, gas.Cv) )
        state.T = state.p / (gas.R_fn(gas.Cp, gas.Cv) * state.Q[:,:,0])

        # update covariant velocities
        #state = covariant(mesh, state)
        soln_vars.calc_covariant(mesh.s_proj, state.Q[:,:,1]/state.Q[:,:,0], state.Q[:,:,2]/state.Q[:,:,0], state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, boundary, parameters, state, gas)

        # print iteration output
        if n % 10 == 0:
            print('Iteration: ' + str(state.n) + '    ' + str(round(state.res[n-1,0],3)) + '    ' + str(round(state.res[n-1,1],3)) + \
                                                 '    ' + str(round(state.res[n-1,2],3)) + '    ' + str(round(state.res[n-1,3],3)) )
            t.toc('Iteration time:')

    print('________________________________________________________________________________________________________________________________________')

    # post processing variables
    state = calc_postvars(state, gas)

    return state


# AUSM+up flux vector splitting scheme
def AUSMplusup( domain, mesh, boundary, parameters, state, gas ):

    print('AUSM+up Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # initialize speed of sound arrays
    c_half_zeta = np.zeros( (domain.M+1, domain.N+2), dtype='float', order='F' )
    c_half_eta = np.zeros( (domain.M+2, domain.N+1), dtype='float', order='F' )

    # initialize phi and p state vectors
    Phi = np.zeros( (domain.M+2, domain.N+2, 4) )
    Phi[:,:,0] = np.ones( (domain.M+2, domain.N+2) )

    P_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    P_zeta[:,:,0] = np.zeros( (domain.M+1, domain.N+2) )
    P_zeta[:,:,3] = np.zeros( (domain.M+1, domain.N+2) )

    P_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    P_eta[:,:,0] = np.zeros( (domain.M+2, domain.N+1) )
    P_eta[:,:,3] = np.zeros( (domain.M+2, domain.N+1) )

    E_hat_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    E_hat_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( [domain.M, domain.N, 4], dtype='float', order='F' )

    if state.n == 0:
        state.res = np.ones( [parameters.iterations + 1, 4] )
    else:
        state.res = np.vstack( (state.res, np.ones( [parameters.iterations + 1, 4] )) )

    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])

    t = TicToc()

    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        t.tic()

        n = n+1
        state.n = state.n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # local timestepping
        state = local_timestep( domain, mesh, state, parameters, gas )

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        # state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / state.Q[:,:,0] + state.p/state.Q[:,:,0]
        state.ht = gas.ht_fn( gas, state )
        # state.Q[:,:,3] = state.Q[:,:,0]*(state.ht - state.p/state.Q[:,:,0])

        # density at cell interfaces, upwinded
        rho_half_zeta = ( state.Q[0:-1,:,0] + state.Q[1:,:,0] ) / 2
        rho_half_eta =  ( state.Q[:,0:-1,0] + state.Q[:,1:,0] ) / 2

        # speed of sound at cell interfaces
        # from Liou 2006 (JCP 214)

        # c_st = thermo.calc_c_star( state.ht, gas.gamma_fn(gas.Cp, gas.Cv) )

        # c_L = c_st[0:-1,:] / np.maximum( np.sqrt(c_st[0:-1,:]), state.U[0:-1,:] )
        # c_R = c_st[1:,:] / np.maximum( np.sqrt(c_st[1:,:]), -state.U[1:,:] ) 
        # c_D = c_st[:,0:-1] / np.maximum( np.sqrt(c_st[:,0:-1]), state.V[:,0:-1] )
        # c_U = c_st[:,1:] / np.maximum( np.sqrt(c_st[:,1:]), -state.V[:,1:] )

        # c_half_zeta = np.minimum( c_L, c_R )
        # c_half_eta =  np.minimum( c_D, c_U )

        flux.c_entropy_fixed( c_half_zeta, c_half_eta, state.ht, gas.gamma_fn(gas.Cp, gas.Cv), \
                              state.U, state.V, domain.M, domain.N )

        # cell face Mach numbers
        M_L = state.U[0:-1,:] / c_half_zeta
        M_R = state.U[1:,:] / c_half_zeta
        M_D = state.V[:,0:-1] / c_half_eta
        M_U = state.V[:,1:] / c_half_eta

        # local mean mach number in each computational direction (Liou 2005 eq. 13)
    
        M_bar_zeta = np.sqrt((1/2)*(state.U[0:-1,:]**2 + state.U[1:,:]**2)/c_half_zeta)
        M_bar_eta = np.sqrt((1/2)*(state.V[:,0:-1]**2 + state.V[:,1:]**2)/c_half_eta)

        # global cutoff Mach numbers
        M_co_zeta = parameters.M_in
        M_co_eta = 0.25

        # split interface Mach numbers in the zeta and eta directions
        M_half_zeta = split.M4p( M_L ) + split.M4m( M_R ) + split.Mp(M_bar_zeta, split.M0(M_bar_zeta, M_co_zeta), \
                                                                     state.p[0:-1,:], state.p[1:,:], rho_half_zeta, c_half_zeta)
        M_half_eta = split.M4p( M_D ) + split.M4m( M_U ) + split.Mp(M_bar_eta, split.M0(M_bar_eta, M_co_eta), \
                                                                     state.p[:,0:-1], state.p[:,1:], rho_half_eta, c_half_eta)

        # calculate mass flux at cell interfaces
        mdot_half_zeta = c_half_zeta * M_half_zeta * ( np.double(M_half_zeta>0) * state.Q[0:-1,:,0] + np.double(M_half_zeta<=0) * state.Q[1:,:,0] )
        mdot_half_eta =  c_half_eta  * M_half_eta  * ( np.double(M_half_eta>0)  * state.Q[:,0:-1,0] + np.double(M_half_eta<=0)  * state.Q[:,1:,0] )

        cr = 1e-60
        # calculate pressure flux at cell interfaces
        p_half_zeta = split.P5p( M_L+cr, 3/16 )*state.p[0:-1,:] + split.P5m( M_R+cr, 3/16 )*state.p[1:,:] + \
                      split.Pu( M_L+cr, M_R+cr, state.U[0:-1,:], state.U[1:,:], state.Q[0:-1,:,0], state.Q[1:,:,0], c_half_zeta, 3/16 )
        p_half_eta =  split.P5p( M_D+cr, 3/16 )*state.p[:,0:-1] + split.P5m( M_U+cr, 3/16 )*state.p[:,1:] + \
                      split.Pu( M_D+cr, M_U+cr, state.V[:,0:-1], state.V[:,1:], state.Q[:,0:-1,0], state.Q[:,1:,0], c_half_eta, 3/16 )

        # initialize Phi vector components
        Phi[:,:,1] = state.u
        Phi[:,:,2] = state.v
        Phi[:,:,3] = state.ht
        
        # pressure flux vector
        P_zeta[:,:,1] = p_half_zeta * mesh.s_proj[0:-1,:,0] / mesh.s_proj[0:-1,:,4]
        P_zeta[:,:,2] = p_half_zeta * mesh.s_proj[0:-1,:,1] / mesh.s_proj[0:-1,:,4]

        P_eta[:,:,1] = p_half_eta * mesh.s_proj[:,0:-1,2] / mesh.s_proj[:,0:-1,5]
        P_eta[:,:,2] = p_half_eta * mesh.s_proj[:,0:-1,3] / mesh.s_proj[:,0:-1,5]

        # flux vector reconstruction
        flux.face_flux( mdot_half_zeta, mdot_half_eta, Phi, P_zeta, P_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, domain.M, domain.N)

        # state.residual = -np.dstack( (state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1]) ) * \
        #                           ( ( E_hat_right * np.dstack([mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4]]) \
        #                             - E_hat_left * np.dstack([mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4]]) ) + \
        #                             ( F_hat_top * np.dstack([mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5]])\
        #                             - F_hat_bot * np.dstack([mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5]]) ) )
        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj, domain.M, domain.N ) 
        state.Q[1:-1,1:-1,:] = state.Qn[1:-1,1:-1,:] + state.residual / mesh.dV4

        # L_inf-norm residual
        state.res[state.n-1,0] = np.log10( np.max(state.residual[:,:,0] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,1] = np.log10( np.max(state.residual[:,:,1] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,2] = np.log10( np.max(state.residual[:,:,2] * mesh.dV[1:-1,1:-1]) )
        state.res[state.n-1,3] = np.log10( np.max(state.residual[:,:,3] * mesh.dV[1:-1,1:-1]) ) 

        #state.res[n-1] = np.log10( np.max(state.residual * mesh.dV4) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv) )
        state.T = state.p / (gas.R_fn(gas.Cp, gas.Cv) * state.Q[:,:,0])

        # update covariant velocities
        #state = covariant(mesh, state)
        soln_vars.calc_covariant(mesh.s_proj, state.u, state.v, state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, boundary, parameters, state, gas)

        # print iteration output
        if n % 10 == 0:
            print('Iteration: ' + str(state.n) + '    ' + str(round(state.res[n-1,0],3)) + '    ' + str(round(state.res[n-1,1],3)) + \
                                                 '    ' + str(round(state.res[n-1,2],3)) + '    ' + str(round(state.res[n-1,3],3)) )
            t.toc('Iteration time:')

    print('________________________________________________________________________________________________________________________________________')

    # post processing variables
    state = calc_postvars(state, gas)

    return state


# AUSMDV flux vector splitting scheme 
def AUSMDV( domain, mesh, boundary, parameters, state, gas ):

    print('AUSMDV Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # initialize speed of sound arrays
    c_half_zeta = np.zeros( (domain.M+1, domain.N+2), dtype='float', order='F' )
    c_half_eta = np.zeros( (domain.M+2, domain.N+1), dtype='float', order='F' )

    # initialize phi and p state vectors
    Phi = np.zeros( (domain.M+2, domain.N+2, 4) )
    Phi[:,:,0] = np.ones( (domain.M+2, domain.N+2) )

    P_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    P_zeta[:,:,0] = np.zeros( (domain.M+1, domain.N+2) )
    P_zeta[:,:,3] = np.zeros( (domain.M+1, domain.N+2) )

    P_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    P_eta[:,:,0] = np.zeros( (domain.M+2, domain.N+1) )
    P_eta[:,:,3] = np.zeros( (domain.M+2, domain.N+1) )

    E_hat_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    E_hat_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( [domain.M, domain.N, 4], dtype='float', order='F' )
    if state.n == 0:
        state.res = np.ones( [parameters.iterations + 1, 4] )
    else:
        state.res = np.vstack( (state.res, np.ones( [parameters.iterations + 1, 4] )) )
        
    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])

    t = TicToc()


    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        t.tic()

        n = n+1
        state.n = state.n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # local timestepping
        state = local_timestep( domain, mesh, state, parameters, gas )

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        # state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / state.Q[:,:,0] + state.p/state.Q[:,:,0]
        state.ht = gas.ht_fn( gas, state )
        # state.Q[:,:,3] = state.Q[:,:,0]*(state.ht - state.p/state.Q[:,:,0])

        # speed of sound at cell interfaces
        # from Liou 2006 (JCP 214)

        flux.c_entropy_fixed( c_half_zeta, c_half_eta, state.ht, gas.gamma_fn(gas.Cp, gas.Cv), \
                              state.U, state.V, domain.M, domain.N )

        # cell face Mach numbers
        M_L = state.U[0:-1,:] / c_half_zeta
        M_R = state.U[1:,:] / c_half_zeta
        M_D = state.V[:,0:-1] / c_half_eta
        M_U = state.V[:,1:] / c_half_eta

        # blending functions
        omega_zeta_plus = 2*(state.p[0:-1,:]/state.Q[0:-1,:,0]) / ((state.p[0:-1,:]/state.Q[0:-1,:,0]) + (state.p[1:,:]/state.Q[1:,:,0]))
        omega_zeta_minus = 2*(state.p[1:,:]/state.Q[1:,:,0]) / ((state.p[0:-1,:]/state.Q[0:-1,:,0]) + (state.p[1:,:]/state.Q[1:,:,0]))
        omega_eta_plus = 2*(state.p[:,0:-1]/state.Q[:,0:-1,0]) / ((state.p[:,0:-1]/state.Q[:,0:-1,0]) + (state.p[:,1:]/state.Q[:,1:,0]))
        omega_eta_minus = 2*(state.p[:,1:]/state.Q[:,1:,0]) / ((state.p[:,0:-1]/state.Q[:,0:-1,0]) + (state.p[:,1:]/state.Q[:,1:,0]))

        # split interface Mach numbers in the zeta and eta directions
        M_half_zeta_plus = omega_zeta_plus*split.M4p( M_L ) + (1-omega_zeta_plus)*split.M1p( M_L )
        M_half_zeta_minus = omega_zeta_minus*split.M4m( M_R ) + (1-omega_zeta_minus)*split.M1m( M_R )
        M_half_eta_plus = omega_eta_plus*split.M4p( M_D ) + (1-omega_eta_plus)*split.M1p( M_D )
        M_half_eta_minus = omega_eta_minus*split.M4m( M_U ) + (1-omega_eta_minus)*split.M1m( M_U )

        # calculate mass flux at cell interfaces
        mdot_half_zeta = c_half_zeta * ( state.Q[0:-1,:,0]*M_half_zeta_plus + state.Q[1:,:,0]*M_half_zeta_minus )
        mdot_half_eta = c_half_eta * ( state.Q[:,0:-1,0]*M_half_eta_plus + state.Q[:,1:,0]*M_half_eta_minus )

        cr = 1e-60
        # calculate pressure flux at cell interfaces
        p_half_zeta = split.P1p( M_L+cr )*state.p[0:-1,:] + split.P1m( M_R+cr )*state.p[1:,:]
        p_half_eta =  split.P1p( M_D+cr )*state.p[:,0:-1] + split.P1m( M_U+cr )*state.p[:,1:]

        # initialize Phi vector components
        Phi[:,:,1] = state.u
        Phi[:,:,2] = state.v
        Phi[:,:,3] = state.ht
        
        # pressure flux vector
        P_zeta[:,:,1] = p_half_zeta * mesh.s_proj[0:-1,:,0] / mesh.s_proj[0:-1,:,4]
        P_zeta[:,:,2] = p_half_zeta * mesh.s_proj[0:-1,:,1] / mesh.s_proj[0:-1,:,4]

        P_eta[:,:,1] = p_half_eta * mesh.s_proj[:,0:-1,2] / mesh.s_proj[:,0:-1,5]
        P_eta[:,:,2] = p_half_eta * mesh.s_proj[:,0:-1,3] / mesh.s_proj[:,0:-1,5]

        # flux vector reconstruction
        flux.face_flux( mdot_half_zeta, mdot_half_eta, Phi, P_zeta, P_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, domain.M, domain.N)

        # state.residual = -np.dstack( (state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1]) ) * \
        #                           ( ( E_hat_right * np.dstack([mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4]]) \
        #                             - E_hat_left * np.dstack([mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4]]) ) + \
        #                             ( F_hat_top * np.dstack([mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5]])\
        #                             - F_hat_bot * np.dstack([mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5]]) ) )

        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj, domain.M, domain.N ) 
        state.Q[1:-1,1:-1,:] = state.Qn[1:-1,1:-1,:] + state.residual / mesh.dV4

        # L_inf-norm residual
        state.res[state.n-1,0] = np.log10( np.max(state.residual[:,:,0] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,1] = np.log10( np.max(state.residual[:,:,1] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,2] = np.log10( np.max(state.residual[:,:,2] * mesh.dV[1:-1,1:-1]) )
        state.res[state.n-1,3] = np.log10( np.max(state.residual[:,:,3] * mesh.dV[1:-1,1:-1]) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv) )
        state.T = state.p / (gas.R_fn(gas.Cp, gas.Cv) * state.Q[:,:,0])

        # update covariant velocities
        #state = covariant(mesh, state)
        soln_vars.calc_covariant(mesh.s_proj, state.u, state.v, state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, boundary, parameters, state, gas)

        # print iteration output
        if n % 10 == 0:
            print('Iteration: ' + str(state.n) + '    ' + str(round(state.res[n-1,0],3)) + '    ' + str(round(state.res[n-1,1],3)) + \
                                                 '    ' + str(round(state.res[n-1,2],3)) + '    ' + str(round(state.res[n-1,3],3)) )
            t.toc('Iteration time:')

    print('________________________________________________________________________________________________________________________________________')


    # post processing variables
    state = calc_postvars(state, gas)
    
    return state


# Simple Low-dissipation AUSM (SLAU) scheme {Shima and Kitamura 2009}
def SLAU( domain, mesh, boundary, parameters, state, gas ):

    print('SLAU Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # initialize speed of sound arrays
    c_half_zeta = np.zeros( (domain.M+1, domain.N+2), dtype='float', order='F' )
    c_half_eta = np.zeros( (domain.M+2, domain.N+1), dtype='float', order='F' )

    # initialize phi and p state vectors
    Phi = np.zeros( (domain.M+2, domain.N+2, 4) )
    Phi[:,:,0] = np.ones( (domain.M+2, domain.N+2) )

    P_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    P_zeta[:,:,0] = np.zeros( (domain.M+1, domain.N+2) )
    P_zeta[:,:,3] = np.zeros( (domain.M+1, domain.N+2) )

    P_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    P_eta[:,:,0] = np.zeros( (domain.M+2, domain.N+1) )
    P_eta[:,:,3] = np.zeros( (domain.M+2, domain.N+1) )

    E_hat_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    E_hat_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( [domain.M, domain.N, 4], dtype='float', order='F' )
    if state.n == 0:
        state.res = np.ones( [parameters.iterations + 1, 4] )
    else:
        state.res = np.vstack( (state.res, np.ones( [parameters.iterations + 1, 4] )) )

    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])

    t = TicToc()

    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        t.tic()

        n = n+1
        state.n = state.n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # local timestepping
        state = local_timestep( domain, mesh, state, parameters, gas )

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        # state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / state.Q[:,:,0] + state.p/state.Q[:,:,0]
        state.ht = gas.ht_fn( gas, state )
        # state.Q[:,:,3] = state.Q[:,:,0]*(state.ht - state.p/state.Q[:,:,0])

        # density at cell interfaces, upwinded
        rho_half_zeta = ( state.Q[0:-1,:,0] + state.Q[1:,:,0] ) / 2
        rho_half_eta =  ( state.Q[:,0:-1,0] + state.Q[:,1:,0] ) / 2

        # density at cell centroids
        rhoL = state.Q[0:-1,:,0]
        rhoR = state.Q[1:,:,0]
        rhoD = state.Q[:,0:-1,0]
        rhoU = state.Q[:,1:,0]

        # speed of sound at cell interfaces
        # from Liou 2006 (JCP 214)
        flux.c_entropy_fixed( c_half_zeta, c_half_eta, state.ht, gas.gamma_fn(gas.Cp, gas.Cv), \
                              state.U, state.V, domain.M, domain.N )

        # cell face Mach numbers
        M_L = state.U[0:-1,:] / c_half_zeta
        M_R = state.U[1:,:] / c_half_zeta
        M_D = state.V[:,0:-1] / c_half_eta
        M_U = state.V[:,1:] / c_half_eta

        # corrective factor g in each computational direction
        g_zeta = -np.max( np.min(split.M4p( M_L ),0), -1 ) * np.min( np.max( split.M4m( M_R ), 1) )
        g_eta = -np.max( np.min(split.M4p( M_D ),0), -1 ) * np.min( np.max( split.M4m( M_U ), 1) )

        # abs(Vbar_n), density weighted covariant velocity, A8 in Shima and Kitamura 2009
        Ubar = ( rhoR*np.abs(state.U[1:,:]) + rhoL*np.abs(state.U[0:-1,:]) ) / ( rhoL + rhoR )
        Vbar = ( rhoU*np.abs(state.V[:,1:]) + rhoD*np.abs(state.V[:,0:-1]) ) / ( rhoD + rhoU )

        # corrected covariant terms
        U_bar_plus = (1-g_zeta)*Ubar + g_zeta*state.U[1:,:]
        U_bar_minus = (1-g_zeta)*Ubar + g_zeta*state.U[0:-1,:]
        V_bar_plus = (1-g_eta)*Vbar + g_eta*state.V[:,1:]
        V_bar_minus = (1-g_eta)*Vbar + g_eta*state.V[:,0:-1]

        # non-dimensional function, O(M)
        M_bar_zeta = np.minimum( 1.0, (1/c_half_zeta) * np.sqrt( ( state.u[1:,:]**2 + state.v[1:,:]**2 + state.u[0:-1,:]**2 + state.v[0:-1,:]**2 ) / 2 ) )
        M_bar_eta = np.minimum( 1.0, (1/c_half_eta) * np.sqrt( ( state.u[:,1:]**2 + state.v[:,1:]**2 + state.u[:,0:-1]**2 + state.v[:,0:-1]**2 ) / 2 ) )

        chi_zeta = ( 1 - M_bar_zeta ) ** 2
        chi_eta = ( 1 - M_bar_eta ) ** 2

        # calculate mass flux at cell interfaces
        mdot_half_zeta = (1/2) * ( rhoR*( state.U[1:,:]+U_bar_plus ) + rhoL*(state.U[0:-1,:]+U_bar_minus) - (chi_zeta/c_half_zeta)*(state.p[1:,:]-state.p[0:-1,:]) )
        mdot_half_eta = (1/2) * ( rhoU*( state.V[:,1:]+V_bar_plus ) + rhoD*(state.V[:,0:-1]+V_bar_minus) - (chi_eta/c_half_eta)*(state.p[:,1:]-state.p[:,0:-1]) )

        cr = 1e-60
        # calculate pressure flux at cell interfaces
        alpha = 0
        p_half_zeta = (state.p[0:-1,:]+state.p[1:,:])/2 + (((split.P5p( M_L+cr, alpha )+split.P5m( M_R+cr, alpha )))/2)*(state.p[0:-1,:]+state.p[1:,:]) + \
                      (1 - chi_zeta) * (split.P5p( M_L+cr, alpha )+split.P5m( M_R+cr, alpha )-1 )*(state.p[0:-1,:]+state.p[1:,:])/2
        p_half_eta = (state.p[:,0:-1]+state.p[:,1:])/2 + (((split.P5p( M_D+cr, alpha )+split.P5m( M_U+cr, alpha )))/2)*(state.p[:,0:-1]+state.p[:,1:]) + \
                      (1 - chi_eta) * (split.P5p( M_D+cr, alpha )+split.P5m( M_U+cr, alpha )-1 )*(state.p[:,0:-1]+state.p[:,1:])/2


        # initialize Phi vector components
        Phi[:,:,1] = state.u
        Phi[:,:,2] = state.v
        Phi[:,:,3] = state.ht
        
        # pressure flux vector
        P_zeta[:,:,1] = p_half_zeta * mesh.s_proj[0:-1,:,0] / mesh.s_proj[0:-1,:,4]
        P_zeta[:,:,2] = p_half_zeta * mesh.s_proj[0:-1,:,1] / mesh.s_proj[0:-1,:,4]

        P_eta[:,:,1] = p_half_eta * mesh.s_proj[:,0:-1,2] / mesh.s_proj[:,0:-1,5]
        P_eta[:,:,2] = p_half_eta * mesh.s_proj[:,0:-1,3] / mesh.s_proj[:,0:-1,5]

        # flux vector reconstruction
        flux.face_flux( mdot_half_zeta, mdot_half_eta, Phi, P_zeta, P_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, domain.M, domain.N)

        # state.residual = -np.dstack( (state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1]) ) * \
        #                           ( ( E_hat_right * np.dstack([mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4]]) \
        #                             - E_hat_left * np.dstack([mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4]]) ) + \
        #                             ( F_hat_top * np.dstack([mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5]])\
        #                             - F_hat_bot * np.dstack([mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5]]) ) )

        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj, domain.M, domain.N ) 
        state.Q[1:-1,1:-1,:] = state.Qn[1:-1,1:-1,:] + state.residual / mesh.dV4

        # L_inf-norm residual
        state.res[state.n-1,0] = np.log10( np.max(state.residual[:,:,0] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,1] = np.log10( np.max(state.residual[:,:,1] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,2] = np.log10( np.max(state.residual[:,:,2] * mesh.dV[1:-1,1:-1]) )
        state.res[state.n-1,3] = np.log10( np.max(state.residual[:,:,3] * mesh.dV[1:-1,1:-1]) ) 

        #state.res[n-1] = np.log10( np.max(state.residual * mesh.dV4) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv) )
        state.T = state.p / (gas.R_fn(gas.Cp, gas.Cv) * state.Q[:,:,0])

        # update covariant velocities
        #state = covariant(mesh, state)
        soln_vars.calc_covariant(mesh.s_proj, state.u, state.v, state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, boundary, parameters, state, gas)

        # print iteration output
        if n % 10 == 0:
            print('Iteration: ' + str(state.n) + '    ' + str(round(state.res[n-1,0],3)) + '    ' + str(round(state.res[n-1,1],3)) + \
                                                 '    ' + str(round(state.res[n-1,2],3)) + '    ' + str(round(state.res[n-1,3],3)) )
            t.toc('Iteration time:')

    print('________________________________________________________________________________________________________________________________________')

    # post processing variables
    state = calc_postvars(state, gas)

    return state


# AUSM flux vector splitting scheme with MUSCL interpolation
def AUSMmuscl( domain, mesh, boundary, parameters, state, gas ):

    print('AUSM Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # initialize speed of sound arrays
    c_half_zeta = np.zeros( (domain.M+1, domain.N+2), dtype='float', order='F' )
    c_half_eta = np.zeros( (domain.M+2, domain.N+1), dtype='float', order='F' )

    # initialize phi and p state vectors
    Phi = np.zeros( (domain.M+2, domain.N+2, 4) )
    Phi[:,:,0] = np.ones( (domain.M+2, domain.N+2) )

    PhiL = np.zeros( (domain.M+1, domain.N+2, 4) )
    PhiL[:,:,0] = np.ones( (domain.M+1, domain.N+2) )
    PhiR = np.zeros( (domain.M+1, domain.N+2, 4) )
    PhiR[:,:,0] = np.ones( (domain.M+1, domain.N+2) )
    PhiB = np.zeros( (domain.M+2, domain.N+1, 4) )
    PhiB[:,:,0] = np.ones( (domain.M+2, domain.N+1) )
    PhiT = np.zeros( (domain.M+2, domain.N+1, 4) )
    PhiT[:,:,0] = np.ones( (domain.M+2, domain.N+1) )

    P_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    P_zeta[:,:,0] = np.zeros( (domain.M+1, domain.N+2) )
    P_zeta[:,:,3] = np.zeros( (domain.M+1, domain.N+2) )

    P_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    P_eta[:,:,0] = np.zeros( (domain.M+2, domain.N+1) )
    P_eta[:,:,3] = np.zeros( (domain.M+2, domain.N+1) )

    E_hat_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    E_hat_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( [domain.M, domain.N, 4], dtype='float', order='F' )
    if state.n == 0:
        state.res = np.ones( [parameters.iterations + 1, 4] )
    else:
        state.res = np.vstack( (state.res, np.ones( [parameters.iterations + 1, 4] )) )

    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])

    t = TicToc()

    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        t.tic()

        n = n+1
        state.n = state.n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # local timestepping
        state = local_timestep( domain, mesh, state, parameters, gas )

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        # state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / \
        #                               state.Q[:,:,0] + state.p/state.Q[:,:,0]
        state.ht = gas.ht_fn( gas, state )
        # state.Q[:,:,3] = state.Q[:,:,0]*(state.ht - state.p/state.Q[:,:,0])

        # MUSCL interpolation

        QL, QR, QB, QT = muscl.MUSCL( state.Q, parameters.epsilon, parameters.kappa, parameters.limiter )

        uL = QL[:,:,1] / QL[:,:,0]
        uR = QR[:,:,1] / QR[:,:,0]
        uB = QB[:,:,1] / QB[:,:,0]
        uT = QT[:,:,1] / QT[:,:,0]

        vL = QL[:,:,2] / QL[:,:,0]
        vR = QR[:,:,2] / QR[:,:,0]
        vB = QB[:,:,2] / QB[:,:,0]
        vT = QT[:,:,2] / QT[:,:,0]

        # soln_vars.calc_covariant(mesh.s_proj[0:domain.M+1,:,:], uL, vL, UL, VL, domain.M+1, domain.N+2)
        # soln_vars.calc_covariant(mesh.s_proj[1:domain.M+2,:,:], uR, vR, UR, VR, domain.M+1, domain.N+2)
        # soln_vars.calc_covariant(mesh.s_proj[:,0:domain.N+1,:], uB, vB, UB, VB, domain.M+2, domain.N+1)
        # soln_vars.calc_covariant(mesh.s_proj[:,1:domain.N+2,:], uT, vT, UT, VT, domain.M+2, domain.N+1)

        UL = (1 / mesh.s_proj[0:domain.M+1,:, 4]) * ( (QL[:,:,1]/QL[:,:,0])*mesh.s_proj[0:domain.M+1,:,0] + \
                                                      (QL[:,:,2]/QL[:,:,0])*mesh.s_proj[0:domain.M+1,:,1] )
        UR = (1 / mesh.s_proj[1:domain.M+2,:,4]) * ( (QR[:,:,1]/QR[:,:,0])*mesh.s_proj[1:domain.M+2,:,0] + \
                                                     (QR[:,:,2]/QR[:,:,0])*mesh.s_proj[1:domain.M+2,:,1] ) 
        VB = (1 / mesh.s_proj[:, 0:domain.N+1, 5]) * ( (QB[:,:,1]/QB[:,:,0])*mesh.s_proj[:,0:domain.N+1,2] + \
                                                       (QB[:,:,2]/QB[:,:,0])*mesh.s_proj[:,0:domain.N+1,3] )
        VT = (1 / mesh.s_proj[:, 1:domain.N+2, 5]) * ( (QT[:,:,1]/QT[:,:,0])*mesh.s_proj[:,1:domain.N+2,2] + \
                                                       (QT[:,:,2]/QT[:,:,0])*mesh.s_proj[:,1:domain.N+2,3] )   

        pL = thermo.calc_p( QL[:,:,0], QL[:,:,3], uL, vL, \
                            gas.gamma_fn(gas.Cp[0:domain.M+1,:], gas.Cv[0:domain.M+1,:]) )
        pR = thermo.calc_p( QR[:,:,0], QR[:,:,3], uR, vR, \
                            gas.gamma_fn(gas.Cp[1:domain.M+2,:], gas.Cv[1:domain.M+2,:]) )
        pB = thermo.calc_p( QB[:,:,0], QB[:,:,3], uB, vB, \
                            gas.gamma_fn(gas.Cp[:,0:domain.N+1], gas.Cv[:,0:domain.N+1]) )
        pT = thermo.calc_p( QT[:,:,0], QT[:,:,3], uT, vT, \
                            gas.gamma_fn(gas.Cp[:,1:domain.N+2], gas.Cv[:,1:domain.N+2]) )     

        # speed of sound at cell interfaces
        # from Liou 2006 (JCP 214)

        # c_st = thermo.calc_c_star( state.ht, gas.gamma_fn(gas.Cp, gas.Cv) )

        c_L = thermo.calc_c_star(QL[:,:,3]/QL[:,:,0], gas.gamma_fn(gas.Cp[0:-1,:], gas.Cv[0:-1,:])) / \
                                 np.maximum( np.sqrt(thermo.calc_c_star(QL[:,:,3]/QL[:,:,0], gas.gamma_fn(gas.Cp[0:-1,:], gas.Cv[0:-1,:]))), UL )
        c_R = thermo.calc_c_star(QR[:,:,3]/QR[:,:,0], gas.gamma_fn(gas.Cp[1:,:], gas.Cv[1:,:])) / \
                                 np.maximum( np.sqrt(thermo.calc_c_star(QR[:,:,3]/QR[:,:,0], gas.gamma_fn(gas.Cp[1:,:], gas.Cv[1:,:]))), -UR )
        c_B = thermo.calc_c_star(QB[:,:,3]/QB[:,:,0], gas.gamma_fn(gas.Cp[:,0:-1], gas.Cv[:,0:-1])) / \
                                 np.maximum( np.sqrt(thermo.calc_c_star(QB[:,:,3]/QB[:,:,0], gas.gamma_fn(gas.Cp[:,0:-1], gas.Cv[:,0:-1]))), VB )
        c_T = thermo.calc_c_star(QT[:,:,3]/QT[:,:,0], gas.gamma_fn(gas.Cp[:,1:], gas.Cv[:,1:])) / \
                                 np.maximum( np.sqrt(thermo.calc_c_star(QT[:,:,3]/QT[:,:,0], gas.gamma_fn(gas.Cp[:,1:], gas.Cv[:,1:]))), -VT )

        c_half_zeta = np.minimum( c_L, c_R )
        c_half_eta =  np.minimum( c_B, c_T )

        # flux.c_entropy_fixed( c_half_zeta, c_half_eta, state.ht, gas.gamma_fn(gas.Cp, gas.Cv), \
        #                       state.U, state.V, domain.M, domain.N )

        # cell face Mach numbers
        M_L = state.U[0:-1,:] / c_half_zeta
        M_R = state.U[1:,:] / c_half_zeta
        M_B = state.V[:,0:-1] / c_half_eta
        M_U = state.V[:,1:] / c_half_eta

        # split interface Mach numbers in the zeta and eta directions
        M_half_zeta = split.Mvp( M_L ) + split.Mvm( M_R )
        M_half_eta = split.Mvp( M_B ) + split.Mvm( M_U )

        # calculate mass flux at cell interfaces
        mdot_half_zeta = c_half_zeta * M_half_zeta * \
                       ( np.double(M_half_zeta>0) * QL[:,:,0] + np.double(M_half_zeta<=0) * QR[:,:,0] )
        mdot_half_eta =  c_half_eta  * M_half_eta  * \
                       ( np.double(M_half_eta>0)  * QB[:,:,0] + np.double(M_half_eta<=0)  * QT[:,:,0] )

        cr = 1e-60
        # calculate pressure flux at cell interfaces
        p_half_zeta = split.P1p( M_L+cr )*pL + split.P1m( M_R+cr )*pR
        p_half_eta =  split.P1p( M_B+cr )*pB + split.P1m( M_U+cr )*pT

        # initialize Phi vector components

        PhiT[:,:,1] = QT[:,:,1]/QT[:,:,0]
        PhiT[:,:,2] = QT[:,:,2]/QT[:,:,0]
        PhiT[:,:,3] = QT[:,:,3]/QT[:,:,0] + pT/QT[:,:,0]
        
        PhiB[:,:,1] = QB[:,:,1]/QB[:,:,0]
        PhiB[:,:,2] = QB[:,:,2]/QB[:,:,0]
        PhiB[:,:,3] = QB[:,:,3]/QB[:,:,0] + pB/QB[:,:,0]
        
        PhiL[:,:,1] = QL[:,:,1]/QL[:,:,0]
        PhiL[:,:,2] = QL[:,:,2]/QL[:,:,0]
        PhiL[:,:,3] = QL[:,:,3]/QL[:,:,0] + pL/QL[:,:,0]
        
        PhiR[:,:,1] = QR[:,:,1]/QR[:,:,0]
        PhiR[:,:,2] = QR[:,:,2]/QR[:,:,0]
        PhiR[:,:,3] = QR[:,:,3]/QR[:,:,0] + pR/QR[:,:,0]
        
        # pressure flux vector
        P_zeta[:,:,1] = p_half_zeta * mesh.s_proj[0:-1,:,0] / mesh.s_proj[0:-1,:,4]
        P_zeta[:,:,2] = p_half_zeta * mesh.s_proj[0:-1,:,1] / mesh.s_proj[0:-1,:,4]

        P_eta[:,:,1] = p_half_eta * mesh.s_proj[:,0:-1,2] / mesh.s_proj[:,0:-1,5]
        P_eta[:,:,2] = p_half_eta * mesh.s_proj[:,0:-1,3] / mesh.s_proj[:,0:-1,5]

        # compute flux at each cell face
        flux.face_flux_muscl( mdot_half_zeta, mdot_half_eta, PhiL, PhiR, PhiB, PhiT, P_zeta, P_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, domain.M, domain.N)
        # E_hat_left = (1/2) * mdot_half_zeta[0:-1,1:-1] * ( Phi[0:-2,1:-1,:] + Phi[1:-1,1:-1,:] ) \
        #             -(1/2) * np.abs(mdot_half_zeta[0:-1,1:-1]) * ( Phi[1:-1,1:-1,:] - Phi[0:-2,1:-1,:] ) \
        #                    + P_zeta[0:-1,1:-1,:]
        # E_hat_right= (1/2) * mdot_half_zeta[1:,1:-1] * ( Phi[1:-1,1:-1,:] + Phi[2:,1:-1,:] ) \
        #             -(1/2) * np.abs(mdot_half_zeta[1:,1:-1]) * ( Phi[2:,1:-1,:] - Phi[1:-1,1:-1,:] ) \
        #                    + P_zeta[1:,1:-1,:]
        # F_hat_bot =  (1/2) * mdot_half_eta[1:-1,0:-1] * ( Phi[1:-1,0:-2,:] + Phi[1:-1,1:-1,:] ) \
        #             -(1/2) * np.abs(mdot_half_eta[1:-1,0:-1]) * ( Phi[1:-1,1:-1,:] - Phi[1:-1,0:-2,:] ) \
        #                    + P_eta[1:-1,0:-1,:]
        # F_hat_top =  (1/2) * mdot_half_eta[1:-1,1:] * ( Phi[1:-1,1:-1,:] + Phi[1:-1,2:,:] ) \
        #             -(1/2) * np.abs(mdot_half_eta[1:-1,1:]) * ( Phi[1:-1,2:,:] - Phi[1:-1,1:-1,:] ) \
        #                    + P_eta[1:-1,1:,:]

        # state.residual = -np.dstack( (state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1]) ) * \
        #                           ( ( E_hat_right * np.dstack([mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4]]) \
        #                             - E_hat_left * np.dstack([mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4]]) ) + \
        #                             ( F_hat_top * np.dstack([mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5]])\
        #                             - F_hat_bot * np.dstack([mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5]]) ) )

        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj, domain.M, domain.N ) 
        state.Q[1:-1,1:-1,:] = state.Qn[1:-1,1:-1,:] + state.residual / mesh.dV4

        # L_inf-norm residual
        state.res[state.n-1,0] = np.log10( np.max(state.residual[:,:,0] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,1] = np.log10( np.max(state.residual[:,:,1] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,2] = np.log10( np.max(state.residual[:,:,2] * mesh.dV[1:-1,1:-1]) )
        state.res[state.n-1,3] = np.log10( np.max(state.residual[:,:,3] * mesh.dV[1:-1,1:-1]) ) 

        #state.res[n-1] = np.log10( np.max(state.residual * mesh.dV4) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv) )
        state.T = state.p / (gas.R_fn(gas.Cp, gas.Cv) * state.Q[:,:,0])

        # update covariant velocities
        #state = covariant(mesh, state)
        soln_vars.calc_covariant(mesh.s_proj, state.u, state.v, state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, boundary, parameters, state, gas)

        # print iteration output
        if n % 10 == 0:
            print('Iteration: ' + str(state.n) + '    ' + str(round(state.res[n-1,0],3)) + '    ' + str(round(state.res[n-1,1],3)) + \
                                                 '    ' + str(round(state.res[n-1,2],3)) + '    ' + str(round(state.res[n-1,3],3)) )
            t.toc('Iteration time:')

    print('________________________________________________________________________________________________________________________________________')

    # post processing variables
    state = calc_postvars(state, gas)

    return state


# AUSM+up with MUSCL interpolation
def AUSMplusupmuscl( domain, mesh, boundary, parameters, state, gas ):

    print('AUSM+up Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # initialize speed of sound arrays
    c_half_zeta = np.zeros( (domain.M+1, domain.N+2), dtype='float', order='F' )
    c_half_eta = np.zeros( (domain.M+2, domain.N+1), dtype='float', order='F' )

    # initialize phi and p state vectors
    PhiL = np.zeros( (domain.M+1, domain.N+2, 4) )
    PhiL[:,:,0] = np.ones( (domain.M+1, domain.N+2) )
    PhiR = np.zeros( (domain.M+1, domain.N+2, 4) )
    PhiR[:,:,0] = np.ones( (domain.M+1, domain.N+2) )
    PhiB = np.zeros( (domain.M+2, domain.N+1, 4) )
    PhiB[:,:,0] = np.ones( (domain.M+2, domain.N+1) )
    PhiT = np.zeros( (domain.M+2, domain.N+1, 4) )
    PhiT[:,:,0] = np.ones( (domain.M+2, domain.N+1) )

    P_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    P_zeta[:,:,0] = np.zeros( (domain.M+1, domain.N+2) )
    P_zeta[:,:,3] = np.zeros( (domain.M+1, domain.N+2) )

    P_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    P_eta[:,:,0] = np.zeros( (domain.M+2, domain.N+1) )
    P_eta[:,:,3] = np.zeros( (domain.M+2, domain.N+1) )

    E_hat_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    E_hat_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( [domain.M, domain.N, 4], dtype='float', order='F' )

    if state.n == 0:
        state.res = np.ones( [parameters.iterations + 1, 4] )
    else:
        state.res = np.vstack( (state.res, np.ones( [parameters.iterations + 1, 4] )) )

    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])

    t = TicToc()

    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        t.tic()

        n = n+1
        state.n = state.n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # local timestepping
        state = local_timestep( domain, mesh, state, parameters, gas )

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        # state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / \
        #                               state.Q[:,:,0] + state.p/state.Q[:,:,0]
        state.ht = gas.ht_fn( gas, state )
        # state.Q[:,:,3] = state.Q[:,:,0]*(state.ht - state.p/state.Q[:,:,0])

        # MUSCL interpolation

        QL, QR, QB, QT = muscl.MUSCL( state.Q, parameters.epsilon, parameters.kappa, parameters.limiter )

        uL = QL[:,:,1] / QL[:,:,0]
        uR = QR[:,:,1] / QR[:,:,0]
        uB = QB[:,:,1] / QB[:,:,0]
        uT = QT[:,:,1] / QT[:,:,0]

        vL = QL[:,:,2] / QL[:,:,0]
        vR = QR[:,:,2] / QR[:,:,0]
        vB = QB[:,:,2] / QB[:,:,0]
        vT = QT[:,:,2] / QT[:,:,0]

        UL = (1 / mesh.s_proj[0:domain.M+1,:, 4]) * ( (QL[:,:,1]/QL[:,:,0])*mesh.s_proj[0:domain.M+1,:,0] + \
                                                      (QL[:,:,2]/QL[:,:,0])*mesh.s_proj[0:domain.M+1,:,1] )
        UR = (1 / mesh.s_proj[1:domain.M+2,:,4]) * ( (QR[:,:,1]/QR[:,:,0])*mesh.s_proj[1:domain.M+2,:,0] + \
                                                     (QR[:,:,2]/QR[:,:,0])*mesh.s_proj[1:domain.M+2,:,1] ) 
        VB = (1 / mesh.s_proj[:, 0:domain.N+1, 5]) * ( (QB[:,:,1]/QB[:,:,0])*mesh.s_proj[:,0:domain.N+1,2] + \
                                                       (QB[:,:,2]/QB[:,:,0])*mesh.s_proj[:,0:domain.N+1,3] )
        VT = (1 / mesh.s_proj[:, 1:domain.N+2, 5]) * ( (QT[:,:,1]/QT[:,:,0])*mesh.s_proj[:,1:domain.N+2,2] + \
                                                       (QT[:,:,2]/QT[:,:,0])*mesh.s_proj[:,1:domain.N+2,3] )   

        pL = thermo.calc_p( QL[:,:,0], QL[:,:,3], uL, vL, \
                            gas.gamma_fn(gas.Cp[0:domain.M+1,:], gas.Cv[0:domain.M+1,:]) )
        pR = thermo.calc_p( QR[:,:,0], QR[:,:,3], uR, vR, \
                            gas.gamma_fn(gas.Cp[1:domain.M+2,:], gas.Cv[1:domain.M+2,:]) )
        pB = thermo.calc_p( QB[:,:,0], QB[:,:,3], uB, vB, \
                            gas.gamma_fn(gas.Cp[:,0:domain.N+1], gas.Cv[:,0:domain.N+1]) )
        pT = thermo.calc_p( QT[:,:,0], QT[:,:,3], uT, vT, \
                            gas.gamma_fn(gas.Cp[:,1:domain.N+2], gas.Cv[:,1:domain.N+2]) )             

        # speed of sound at cell interfaces
        # from Liou 2006 (JCP 214)

        # c_st = thermo.calc_c_star( state.ht, gas.gamma_fn(gas.Cp, gas.Cv) )

        c_L = thermo.calc_c_star(QL[:,:,3]/QL[:,:,0], gas.gamma_fn(gas.Cp[0:-1,:], gas.Cv[0:-1,:])) / \
                                 np.maximum( np.sqrt(thermo.calc_c_star(QL[:,:,3]/QL[:,:,0], gas.gamma_fn(gas.Cp[0:-1,:], gas.Cv[0:-1,:]))), UL )
        c_R = thermo.calc_c_star(QR[:,:,3]/QR[:,:,0], gas.gamma_fn(gas.Cp[1:,:], gas.Cv[1:,:])) / \
                                 np.maximum( np.sqrt(thermo.calc_c_star(QR[:,:,3]/QR[:,:,0], gas.gamma_fn(gas.Cp[1:,:], gas.Cv[1:,:]))), -UR )
        c_B = thermo.calc_c_star(QB[:,:,3]/QB[:,:,0], gas.gamma_fn(gas.Cp[:,0:-1], gas.Cv[:,0:-1])) / \
                                 np.maximum( np.sqrt(thermo.calc_c_star(QB[:,:,3]/QB[:,:,0], gas.gamma_fn(gas.Cp[:,0:-1], gas.Cv[:,0:-1]))), VB )
        c_T = thermo.calc_c_star(QT[:,:,3]/QT[:,:,0], gas.gamma_fn(gas.Cp[:,1:], gas.Cv[:,1:])) / \
                                 np.maximum( np.sqrt(thermo.calc_c_star(QT[:,:,3]/QT[:,:,0], gas.gamma_fn(gas.Cp[:,1:], gas.Cv[:,1:]))), -VT )

        c_half_zeta = np.minimum( c_L, c_R )
        c_half_eta =  np.minimum( c_B, c_T )

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / state.Q[:,:,0] + state.p/state.Q[:,:,0]

        # density at cell interfaces, upwinded
        rho_half_zeta = ( QL[:,:,0] + QR[:,:,0] ) / 2
        rho_half_eta =  ( QB[:,:,0] + QT[:,:,0] ) / 2

        # cell face Mach numbers
        M_L = UL / c_half_zeta
        M_R = UR / c_half_zeta
        M_D = VB / c_half_eta
        M_U = VT / c_half_eta

        # local mean mach number in each computational direction (Liou 2005 eq. 13)
    
        M_bar_zeta = np.sqrt((1/2)*(UL**2 + UR**2)/c_half_zeta)
        M_bar_eta = np.sqrt((1/2)*(VB**2 + VT**2)/c_half_eta)

        # global cutoff Mach numbers
        M_co_zeta = parameters.M_in
        M_co_eta = 0.25

        # split interface Mach numbers in the zeta and eta directions
        M_half_zeta = split.M4p( M_L ) + split.M4m( M_R ) + split.Mp(M_bar_zeta, split.M0(M_bar_zeta, M_co_zeta), \
                                                                     pL, pR, rho_half_zeta, c_half_zeta)
        M_half_eta = split.M4p( M_D ) + split.M4m( M_U ) + split.Mp(M_bar_eta, split.M0(M_bar_eta, M_co_eta), \
                                                                     pB, pT, rho_half_eta, c_half_eta)

        # calculate mass flux at cell interfaces
        mdot_half_zeta = c_half_zeta * M_half_zeta * ( np.double(M_half_zeta>0) * QL[:,:,0] + np.double(M_half_zeta<=0) * QR[:,:,0] )
        mdot_half_eta =  c_half_eta  * M_half_eta  * ( np.double(M_half_eta>0)  * QB[:,:,0] + np.double(M_half_eta<=0)  * QT[:,:,0] )

        cr = 1e-60
        # calculate pressure flux at cell interfaces
        p_half_zeta = split.P5p( M_L+cr, 3/16 )*pL + split.P5m( M_R+cr, 3/16 )*pR + \
                      split.Pu( M_L+cr, M_R+cr, UL, UR, QL[:,:,0], QR[:,:,0], c_half_zeta, 3/16 )
        p_half_eta =  split.P5p( M_D+cr, 3/16 )*pB + split.P5m( M_U+cr, 3/16 )*pT + \
                      split.Pu( M_D+cr, M_U+cr, VB, VT, QB[:,:,0], QT[:,:,0], c_half_eta, 3/16 )

        # initialize Phi vector components

        PhiT[:,:,1] = QT[:,:,1]/QT[:,:,0]
        PhiT[:,:,2] = QT[:,:,2]/QT[:,:,0]
        PhiT[:,:,3] = QT[:,:,3]/QT[:,:,0] + pT/QT[:,:,0]
        
        PhiB[:,:,1] = QB[:,:,1]/QB[:,:,0]
        PhiB[:,:,2] = QB[:,:,2]/QB[:,:,0]
        PhiB[:,:,3] = QB[:,:,3]/QB[:,:,0] + pB/QB[:,:,0]
        
        PhiL[:,:,1] = QL[:,:,1]/QL[:,:,0]
        PhiL[:,:,2] = QL[:,:,2]/QL[:,:,0]
        PhiL[:,:,3] = QL[:,:,3]/QL[:,:,0] + pL/QL[:,:,0]
        
        PhiR[:,:,1] = QR[:,:,1]/QR[:,:,0]
        PhiR[:,:,2] = QR[:,:,2]/QR[:,:,0]
        PhiR[:,:,3] = QR[:,:,3]/QR[:,:,0] + pR/QR[:,:,0]
        
        # pressure flux vector
        P_zeta[:,:,1] = p_half_zeta * mesh.s_proj[0:-1,:,0] / mesh.s_proj[0:-1,:,4]
        P_zeta[:,:,2] = p_half_zeta * mesh.s_proj[0:-1,:,1] / mesh.s_proj[0:-1,:,4]

        P_eta[:,:,1] = p_half_eta * mesh.s_proj[:,0:-1,2] / mesh.s_proj[:,0:-1,5]
        P_eta[:,:,2] = p_half_eta * mesh.s_proj[:,0:-1,3] / mesh.s_proj[:,0:-1,5]

        # compute flux at each cell face
        flux.face_flux_muscl( mdot_half_zeta, mdot_half_eta, PhiL, PhiR, PhiB, PhiT, P_zeta, P_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, domain.M, domain.N)

        # state.residual = -np.dstack( (state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1]) ) * \
        #                           ( ( E_hat_right * np.dstack([mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4]]) \
        #                             - E_hat_left * np.dstack([mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4]]) ) + \
        #                             ( F_hat_top * np.dstack([mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5]])\
        #                             - F_hat_bot * np.dstack([mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5]]) ) )
        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj, domain.M, domain.N ) 
        state.Q[1:-1,1:-1,:] = state.Qn[1:-1,1:-1,:] + state.residual / mesh.dV4

        # L_inf-norm residual
        state.res[state.n-1,0] = np.log10( np.max(state.residual[:,:,0] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,1] = np.log10( np.max(state.residual[:,:,1] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,2] = np.log10( np.max(state.residual[:,:,2] * mesh.dV[1:-1,1:-1]) )
        state.res[state.n-1,3] = np.log10( np.max(state.residual[:,:,3] * mesh.dV[1:-1,1:-1]) ) 

        #state.res[n-1] = np.log10( np.max(state.residual * mesh.dV4) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv) )
        state.T = state.p / (gas.R_fn(gas.Cp, gas.Cv) * state.Q[:,:,0])

        # update covariant velocities
        #state = covariant(mesh, state)
        soln_vars.calc_covariant(mesh.s_proj, state.u, state.v, state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, boundary, parameters, state, gas)

        # print iteration output
        if n % 10 == 0:
            print('Iteration: ' + str(state.n) + '    ' + str(round(state.res[n-1,0],3)) + '    ' + str(round(state.res[n-1,1],3)) + \
                                                 '    ' + str(round(state.res[n-1,2],3)) + '    ' + str(round(state.res[n-1,3],3)) )
            t.toc('Iteration time:')

    print('________________________________________________________________________________________________________________________________________')

    # post processing variables
    state = calc_postvars(state, gas)

    return state


# AUSMDV with MUSCL interpolation
def AUSMDVmuscl( domain, mesh, boundary, parameters, state, gas ):

    print('AUSMDV Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # initialize speed of sound arrays
    c_half_zeta = np.zeros( (domain.M+1, domain.N+2), dtype='float', order='F' )
    c_half_eta = np.zeros( (domain.M+2, domain.N+1), dtype='float', order='F' )

    QL = np.zeros( (domain.M+1, domain.N+2, 4) )
    QR = np.zeros( (domain.M+1, domain.N+2, 4) )
    QB = np.zeros( (domain.M+1, domain.N+2, 4) )
    QT = np.zeros( (domain.M+1, domain.N+2, 4) )

    # initialize phi and p state vectors
    Phi = np.zeros( (domain.M+2, domain.N+2, 4) )
    Phi[:,:,0] = np.ones( (domain.M+2, domain.N+2) )

    PhiL = np.zeros( (domain.M+1, domain.N+2, 4) )
    PhiL[:,:,0] = np.ones( (domain.M+1, domain.N+2) )
    PhiR = np.zeros( (domain.M+1, domain.N+2, 4) )
    PhiR[:,:,0] = np.ones( (domain.M+1, domain.N+2) )
    PhiB = np.zeros( (domain.M+2, domain.N+1, 4) )
    PhiB[:,:,0] = np.ones( (domain.M+2, domain.N+1) )
    PhiT = np.zeros( (domain.M+2, domain.N+1, 4) )
    PhiT[:,:,0] = np.ones( (domain.M+2, domain.N+1) )

    P_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    P_zeta[:,:,0] = np.zeros( (domain.M+1, domain.N+2) )
    P_zeta[:,:,3] = np.zeros( (domain.M+1, domain.N+2) )

    P_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    P_eta[:,:,0] = np.zeros( (domain.M+2, domain.N+1) )
    P_eta[:,:,3] = np.zeros( (domain.M+2, domain.N+1) )

    E_hat_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    E_hat_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( [domain.M, domain.N, 4], dtype='float', order='F' )
    if state.n == 0:
        state.res = np.ones( [parameters.iterations + 1, 4] )
    else:
        state.res = np.vstack( (state.res, np.ones( [parameters.iterations + 1, 4] )) )

    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])

    t = TicToc()

    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        t.tic()

        n = n+1
        state.n = state.n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # local timestepping
        state = local_timestep( domain, mesh, state, parameters, gas )

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        # state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / \
        #                               state.Q[:,:,0] + state.p/state.Q[:,:,0]
        state.ht = gas.ht_fn( gas, state )
        # state.Q[:,:,3] = state.Q[:,:,0]*(state.ht - state.p/state.Q[:,:,0])

        # MUSCL interpolation

        QL, QR, QB, QT = muscl.MUSCL( state.Q, parameters.epsilon, parameters.kappa, parameters.limiter )

        UL = (1 / mesh.s_proj[0:domain.M+1,:, 4]) * ( (QL[:,:,1]/QL[:,:,0])*mesh.s_proj[0:domain.M+1,:,0] + \
                                                      (QL[:,:,2]/QL[:,:,0])*mesh.s_proj[0:domain.M+1,:,1] )
        UR = (1 / mesh.s_proj[1:domain.M+2,:,4]) * ( (QR[:,:,1]/QR[:,:,0])*mesh.s_proj[1:domain.M+2,:,0] + \
                                                     (QR[:,:,2]/QR[:,:,0])*mesh.s_proj[1:domain.M+2,:,1] ) 
        VB = (1 / mesh.s_proj[:, 0:domain.N+1, 5]) * ( (QB[:,:,1]/QB[:,:,0])*mesh.s_proj[:,0:domain.N+1,2] + \
                                                       (QB[:,:,2]/QB[:,:,0])*mesh.s_proj[:,0:domain.N+1,3] )
        VT = (1 / mesh.s_proj[:, 1:domain.N+2, 5]) * ( (QT[:,:,1]/QT[:,:,0])*mesh.s_proj[:,1:domain.N+2,2] + \
                                                       (QT[:,:,2]/QT[:,:,0])*mesh.s_proj[:,1:domain.N+2,3] )   

        pL = (gas.gamma_fn(gas.Cp[0:domain.M+1,:], gas.Cv[0:domain.M+1,:])-1) * \
             ( QL[:,:,3] - (1/2)*QL[:,:,0]*((QL[:,:,1]/QL[:,:,0])**2 + (QL[:,:,2]/QL[:,:,0])**2))
        pR = (gas.gamma_fn(gas.Cp[1:domain.M+2,:], gas.Cv[1:domain.M+2,:])-1) * \
             ( QR[:,:,3] - (1/2)*QR[:,:,0]*((QR[:,:,1]/QR[:,:,0])**2 + (QR[:,:,2]/QR[:,:,0])**2))  
        pB = (gas.gamma_fn(gas.Cp[:,0:domain.N+1], gas.Cv[:,0:domain.N+1])-1) * \
             ( QB[:,:,3] - (1/2)*QB[:,:,0]*((QB[:,:,1]/QB[:,:,0])**2 + (QB[:,:,2]/QB[:,:,0])**2))
        pT = (gas.gamma_fn(gas.Cp[:,1:domain.N+2], gas.Cv[:,1:domain.N+2])-1) * \
             ( QT[:,:,3] - (1/2)*QT[:,:,0]*((QT[:,:,1]/QT[:,:,0])**2 + (QT[:,:,2]/QT[:,:,0])**2))                

        # density at cell interfaces, upwinded
        rho_half_zeta = ( QL + QR ) / 2
        rho_half_eta =  ( QB + QT ) / 2

        # speed of sound at cell interfaces
        # from Liou 2006 (JCP 214)

        c_st = thermo.calc_c_star( state.ht, gas.gamma_fn(gas.Cp, gas.Cv) )

        c_L = thermo.calc_c_star(QL[:,:,3]/QL[:,:,0], gas.gamma_fn(gas.Cp[0:-1,:], gas.Cv[0:-1,:])) / \
                                 np.maximum( np.sqrt(thermo.calc_c_star(QL[:,:,3]/QL[:,:,0], gas.gamma_fn(gas.Cp[0:-1,:], gas.Cv[0:-1,:]))), UL )
        c_R = thermo.calc_c_star(QR[:,:,3]/QR[:,:,0], gas.gamma_fn(gas.Cp[1:,:], gas.Cv[1:,:])) / \
                                 np.maximum( np.sqrt(thermo.calc_c_star(QR[:,:,3]/QR[:,:,0], gas.gamma_fn(gas.Cp[1:,:], gas.Cv[1:,:]))), -UR )
        c_B = thermo.calc_c_star(QB[:,:,3]/QB[:,:,0], gas.gamma_fn(gas.Cp[:,0:-1], gas.Cv[:,0:-1])) / \
                                 np.maximum( np.sqrt(thermo.calc_c_star(QB[:,:,3]/QB[:,:,0], gas.gamma_fn(gas.Cp[:,0:-1], gas.Cv[:,0:-1]))), VB )
        c_T = thermo.calc_c_star(QT[:,:,3]/QT[:,:,0], gas.gamma_fn(gas.Cp[:,1:], gas.Cv[:,1:])) / \
                                 np.maximum( np.sqrt(thermo.calc_c_star(QT[:,:,3]/QT[:,:,0], gas.gamma_fn(gas.Cp[:,1:], gas.Cv[:,1:]))), -VT )

        c_half_zeta = np.dstack( ( np.minimum( c_L, c_R ), np.minimum( c_L, c_R ), np.minimum( c_L, c_R ), np.minimum( c_L, c_R ) ) )
        c_half_eta =  np.dstack( ( np.minimum( c_B, c_T ), np.minimum( c_B, c_T ), np.minimum( c_B, c_T ), np.minimum( c_B, c_T ) ) )

        # cell face Mach numbers
        M_L = UL / c_half_zeta[:,:,0]
        M_R = UR / c_half_zeta[:,:,0]
        M_B = VB / c_half_eta[:,:,0]
        M_U = VT / c_half_eta[:,:,0]

        # blending functions
        omega_zeta_plus = 2*(pL/QL[:,:,0]) / ((pL/QL[:,:,0]) + (pR/QR[:,:,0]))
        omega_zeta_minus = 2*(pR/QR[:,:,0]) / ((pL/QL[:,:,0]) + (pR/QR[:,:,0]))
        omega_eta_plus = 2*(pB/QB[:,:,0]) / ((pB/QB[:,:,0]) + (pT/QT[:,:,0]))
        omega_eta_minus = 2*(pT/QT[:,:,0]) / ((pB/QB[:,:,0]) + (pT/QT[:,:,0]))

        # split interface Mach numbers in the zeta and eta directions
        M_half_zeta_plus = np.dstack( ( omega_zeta_plus*split.M4p( M_L ) + (1-omega_zeta_plus)*split.M1p( M_L ), \
                                        omega_zeta_plus*split.M4p( M_L ) + (1-omega_zeta_plus)*split.M1p( M_L ), \
                                        omega_zeta_plus*split.M4p( M_L ) + (1-omega_zeta_plus)*split.M1p( M_L ), \
                                        omega_zeta_plus*split.M4p( M_L ) + (1-omega_zeta_plus)*split.M1p( M_L ) ) )
        M_half_zeta_minus = np.dstack( ( omega_zeta_minus*split.M4m( M_R ) + (1-omega_zeta_minus)*split.M1m( M_R ), \
                                         omega_zeta_minus*split.M4m( M_R ) + (1-omega_zeta_minus)*split.M1m( M_R ), \
                                         omega_zeta_minus*split.M4m( M_R ) + (1-omega_zeta_minus)*split.M1m( M_R ), \
                                         omega_zeta_minus*split.M4m( M_R ) + (1-omega_zeta_minus)*split.M1m( M_R ) ) )
        M_half_eta_plus = np.dstack( ( omega_eta_plus*split.M4p( M_B ) + (1-omega_eta_plus)*split.M1p( M_B ), \
                                       omega_eta_plus*split.M4p( M_B ) + (1-omega_eta_plus)*split.M1p( M_B ), \
                                       omega_eta_plus*split.M4p( M_B ) + (1-omega_eta_plus)*split.M1p( M_B ), \
                                       omega_eta_plus*split.M4p( M_B ) + (1-omega_eta_plus)*split.M1p( M_B ) ) )
        M_half_eta_minus = np.dstack( ( omega_eta_minus*split.M4m( M_U ) + (1-omega_eta_minus)*split.M1m( M_U ), \
                                        omega_eta_minus*split.M4m( M_U ) + (1-omega_eta_minus)*split.M1m( M_U ), \
                                        omega_eta_minus*split.M4m( M_U ) + (1-omega_eta_minus)*split.M1m( M_U ), \
                                        omega_eta_minus*split.M4m( M_U ) + (1-omega_eta_minus)*split.M1m( M_U ) ) )

        # calculate mass flux at cell interfaces
        mdot_half_zeta = c_half_zeta * ( M_half_zeta_plus * QL + \
                                         M_half_zeta_minus * QR )
        mdot_half_zeta = c_half_eta * ( M_half_eta_plus * QB + \
                                         M_half_eta_minus * QT )

        # split interface Mach numbers in the zeta and eta directions
        M_half_zeta = split.Mvp( M_L ) + split.Mvm( M_R )
        M_half_eta = split.Mvp( M_B ) + split.Mvm( M_U )

        # calculate mass flux at cell interfaces
        mdot_half_zeta = c_half_zeta[:,:,0] * M_half_zeta * ( np.double(M_half_zeta>0) * QL[:,:,0] + np.double(M_half_zeta<=0) * QR[:,:,0] )
        mdot_half_eta = c_half_eta[:,:,0]  * M_half_eta  * ( np.double(M_half_eta>0)  * QB[:,:,0] + np.double(M_half_eta<=0)  * QT[:,:,0] )

        cr = 1e-60
        # calculate pressure flux at cell interfaces
        p_half_zeta = split.P1p( M_L+cr )*pL + split.P1m( M_R+cr )*pR
        p_half_eta =  split.P1p( M_B+cr )*pB + split.P1m( M_U+cr )*pT

        # initialize Phi vector components

        PhiT[:,:,1] = QT[:,:,1]/QT[:,:,0]
        PhiT[:,:,2] = QT[:,:,2]/QT[:,:,0]
        PhiT[:,:,3] = QT[:,:,3]/QT[:,:,0] + pT/QT[:,:,0]
        
        PhiB[:,:,1] = QB[:,:,1]/QB[:,:,0]
        PhiB[:,:,2] = QB[:,:,2]/QB[:,:,0]
        PhiB[:,:,3] = QB[:,:,3]/QB[:,:,0] + pB/QB[:,:,0]
        
        PhiL[:,:,1] = QL[:,:,1]/QL[:,:,0]
        PhiL[:,:,2] = QL[:,:,2]/QL[:,:,0]
        PhiL[:,:,3] = QL[:,:,3]/QL[:,:,0] + pL/QL[:,:,0]
        
        PhiR[:,:,1] = QR[:,:,1]/QR[:,:,0]
        PhiR[:,:,2] = QR[:,:,2]/QR[:,:,0]
        PhiR[:,:,3] = QR[:,:,3]/QR[:,:,0] + pR/QR[:,:,0]
        
        # pressure flux vector
        P_zeta[:,:,1] = p_half_zeta * mesh.s_proj[0:-1,:,0] / mesh.s_proj[0:-1,:,4]
        P_zeta[:,:,2] = p_half_zeta * mesh.s_proj[0:-1,:,1] / mesh.s_proj[0:-1,:,4]

        P_eta[:,:,1] = p_half_eta * mesh.s_proj[:,0:-1,2] / mesh.s_proj[:,0:-1,5]
        P_eta[:,:,2] = p_half_eta * mesh.s_proj[:,0:-1,3] / mesh.s_proj[:,0:-1,5]


        flux.face_flux_muscl( mdot_half_zeta, mdot_half_eta, PhiL, PhiR, PhiB, PhiT, P_zeta, P_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, domain.M, domain.N)

        # state.residual = -np.dstack( (state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1]) ) * \
        #                           ( ( E_hat_right * np.dstack([mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4]]) \
        #                             - E_hat_left * np.dstack([mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4]]) ) + \
        #                             ( F_hat_top * np.dstack([mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5]])\
        #                             - F_hat_bot * np.dstack([mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5]]) ) )

        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj, domain.M, domain.N ) 
        state.Q[1:-1,1:-1,:] = state.Qn[1:-1,1:-1,:] + state.residual / mesh.dV4

        # L_inf-norm residual
        state.res[state.n-1,0] = np.log10( np.max(state.residual[:,:,0] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,1] = np.log10( np.max(state.residual[:,:,1] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,2] = np.log10( np.max(state.residual[:,:,2] * mesh.dV[1:-1,1:-1]) )
        state.res[state.n-1,3] = np.log10( np.max(state.residual[:,:,3] * mesh.dV[1:-1,1:-1]) ) 

        #state.res[n-1] = np.log10( np.max(state.residual * mesh.dV4) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv) )
        state.T = state.p / (gas.R_fn(gas.Cp, gas.Cv) * state.Q[:,:,0])

        # update covariant velocities
        #state = covariant(mesh, state)
        soln_vars.calc_covariant(mesh.s_proj, state.u, state.v, state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, boundary, parameters, state, gas)

        # print iteration output
        if n % 10 == 0:
            print('Iteration: ' + str(state.n) + '    ' + str(round(state.res[n-1,0],3)) + '    ' + str(round(state.res[n-1,1],3)) + \
                                           '    ' + str(round(state.res[n-1,2],3)) + '    ' + str(round(state.res[n-1,3],3)) )
            t.toc('Iteration time:')

    print('________________________________________________________________________________________________________________________________________')

    # post processing variables
    state = calc_postvars(state, gas)

    return state


# Simple Low-dissipation AUSM (SLAU) scheme {Shima and Kitamura 2009}
def SLAUmuscl( domain, mesh, boundary, parameters, state, gas ):

    print('SLAU Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # initialize speed of sound arrays
    c_half_zeta = np.zeros( (domain.M+1, domain.N+2), dtype='float', order='F' )
    c_half_eta = np.zeros( (domain.M+2, domain.N+1), dtype='float', order='F' )

    # initialize phi and p state vectors
    Phi = np.zeros( (domain.M+2, domain.N+2, 4) )
    Phi[:,:,0] = np.ones( (domain.M+2, domain.N+2) )

    PhiL = np.zeros( (domain.M+1, domain.N+2, 4) )
    PhiL[:,:,0] = np.ones( (domain.M+1, domain.N+2) )
    PhiR = np.zeros( (domain.M+1, domain.N+2, 4) )
    PhiR[:,:,0] = np.ones( (domain.M+1, domain.N+2) )
    PhiB = np.zeros( (domain.M+2, domain.N+1, 4) )
    PhiB[:,:,0] = np.ones( (domain.M+2, domain.N+1) )
    PhiT = np.zeros( (domain.M+2, domain.N+1, 4) )
    PhiT[:,:,0] = np.ones( (domain.M+2, domain.N+1) )

    P_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    P_zeta[:,:,0] = np.zeros( (domain.M+1, domain.N+2) )
    P_zeta[:,:,3] = np.zeros( (domain.M+1, domain.N+2) )

    P_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    P_eta[:,:,0] = np.zeros( (domain.M+2, domain.N+1) )
    P_eta[:,:,3] = np.zeros( (domain.M+2, domain.N+1) )

    E_hat_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    E_hat_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( [domain.M, domain.N, 4], dtype='float', order='F' )
    if state.n == 0:
        state.res = np.ones( [parameters.iterations + 1, 4] )
    else:
        state.res = np.vstack( (state.res, np.ones( [parameters.iterations + 1, 4] )) )

    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])

    t = TicToc()

    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        t.tic()

        n = n+1
        state.n = state.n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # local timestepping
        state = local_timestep( domain, mesh, state, parameters, gas )


        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        # state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / \
        #                               state.Q[:,:,0] + state.p/state.Q[:,:,0]
        state.ht = gas.ht_fn( gas, state )
        # state.Q[:,:,3] = state.Q[:,:,0]*(state.ht - state.p/state.Q[:,:,0])

        # MUSCL interpolation
        QL, QR, QB, QT = muscl.MUSCL( state.Q, parameters.epsilon, parameters.kappa, parameters.limiter )

        uL = QL[:,:,1] / QL[:,:,0]
        uR = QR[:,:,1] / QR[:,:,0]
        uB = QB[:,:,1] / QB[:,:,0]
        uT = QT[:,:,1] / QT[:,:,0]

        vL = QL[:,:,2] / QL[:,:,0]
        vR = QR[:,:,2] / QR[:,:,0]
        vB = QB[:,:,2] / QB[:,:,0]
        vT = QT[:,:,2] / QT[:,:,0]

        UL = (1 / mesh.s_proj[0:domain.M+1,:, 4]) * ( (QL[:,:,1]/QL[:,:,0])*mesh.s_proj[0:domain.M+1,:,0] + \
                                                      (QL[:,:,2]/QL[:,:,0])*mesh.s_proj[0:domain.M+1,:,1] )
        UR = (1 / mesh.s_proj[1:domain.M+2,:,4]) * ( (QR[:,:,1]/QR[:,:,0])*mesh.s_proj[1:domain.M+2,:,0] + \
                                                     (QR[:,:,2]/QR[:,:,0])*mesh.s_proj[1:domain.M+2,:,1] ) 
        VB = (1 / mesh.s_proj[:, 0:domain.N+1, 5]) * ( (QB[:,:,1]/QB[:,:,0])*mesh.s_proj[:,0:domain.N+1,2] + \
                                                       (QB[:,:,2]/QB[:,:,0])*mesh.s_proj[:,0:domain.N+1,3] )
        VT = (1 / mesh.s_proj[:, 1:domain.N+2, 5]) * ( (QT[:,:,1]/QT[:,:,0])*mesh.s_proj[:,1:domain.N+2,2] + \
                                                       (QT[:,:,2]/QT[:,:,0])*mesh.s_proj[:,1:domain.N+2,3] )   

        pL = (gas.gamma_fn(gas.Cp[0:domain.M+1,:], gas.Cv[0:domain.M+1,:])-1) * \
             ( QL[:,:,3] - (1/2)*QL[:,:,0]*((QL[:,:,1]/QL[:,:,0])**2 + (QL[:,:,2]/QL[:,:,0])**2))
        pR = (gas.gamma_fn(gas.Cp[1:domain.M+2,:], gas.Cv[1:domain.M+2,:])-1) * \
             ( QR[:,:,3] - (1/2)*QR[:,:,0]*((QR[:,:,1]/QR[:,:,0])**2 + (QR[:,:,2]/QR[:,:,0])**2))  
        pB = (gas.gamma_fn(gas.Cp[:,0:domain.N+1], gas.Cv[:,0:domain.N+1])-1) * \
             ( QB[:,:,3] - (1/2)*QB[:,:,0]*((QB[:,:,1]/QB[:,:,0])**2 + (QB[:,:,2]/QB[:,:,0])**2))
        pT = (gas.gamma_fn(gas.Cp[:,1:domain.N+2], gas.Cv[:,1:domain.N+2])-1) * \
             ( QT[:,:,3] - (1/2)*QT[:,:,0]*((QT[:,:,1]/QT[:,:,0])**2 + (QT[:,:,2]/QT[:,:,0])**2))      

        # density at cell interfaces, upwinded
        rho_half_zeta = ( QL[:,:,0] + QR[:,:,0] ) / 2
        rho_half_eta =  ( QB[:,:,0] + QT[:,:,0] ) / 2

        # density at cell centroids
        rhoL = QL[:,:,0]
        rhoR = QR[:,:,0]
        rhoD = QB[:,:,0]
        rhoU = QT[:,:,0]

        # speed of sound at cell interfaces
        # from Liou 2006 (JCP 214)
        flux.c_entropy_fixed( c_half_zeta, c_half_eta, state.ht, gas.gamma_fn(gas.Cp, gas.Cv), \
                              state.U, state.V, domain.M, domain.N )

        # cell face Mach numbers
        M_L = UL / c_half_zeta
        M_R = UR / c_half_zeta
        M_D = VB / c_half_eta
        M_U = VT / c_half_eta

        # corrective factor g in each computational direction
        g_zeta = -np.max( np.min(split.M4p( M_L ),0), -1 ) * np.min( np.max( split.M4m( M_R ), 1) )
        g_eta = -np.max( np.min(split.M4p( M_D ),0), -1 ) * np.min( np.max( split.M4m( M_U ), 1) )

        # abs(Vbar_n), density weighted covariant velocity, A8 in Shima and Kitamura 2009
        Ubar = ( rhoR*np.abs(UR) + rhoL*np.abs(UL) ) / ( rhoL + rhoR )
        Vbar = ( rhoU*np.abs(VT) + rhoD*np.abs(VB) ) / ( rhoD + rhoU )

        # corrected covariant terms
        U_bar_plus = (1-g_zeta)*Ubar + g_zeta*UR
        U_bar_minus = (1-g_zeta)*Ubar + g_zeta*UL
        V_bar_plus = (1-g_eta)*Vbar + g_eta*VT
        V_bar_minus = (1-g_eta)*Vbar + g_eta*VB

        # non-dimensional function, O(M)
        M_bar_zeta = np.minimum( 1.0, (1/c_half_zeta) * np.sqrt( ( uR**2 + vR**2 + uL**2 + vL**2 ) / 2 ) )
        M_bar_eta = np.minimum( 1.0, (1/c_half_eta) * np.sqrt( ( uT**2 + vT**2 + uB**2 + vB**2 ) / 2 ) )

        chi_zeta = ( 1 - M_bar_zeta ) ** 2
        chi_eta = ( 1 - M_bar_eta ) ** 2

        # calculate mass flux at cell interfaces
        mdot_half_zeta = (1/2) * ( rhoR*( UR+U_bar_plus ) + rhoL*(UL+U_bar_minus) - (chi_zeta/c_half_zeta)*(pR-pL) )
        mdot_half_eta = (1/2) * ( rhoU*( VT+V_bar_plus ) + rhoD*(VB+V_bar_minus) - (chi_eta/c_half_eta)*(pT-pB) )

        cr = 1e-60
        # calculate pressure flux at cell interfaces
        alpha = 0
        p_half_zeta = (pL+pR)/2 + (((split.P5p( M_L+cr, alpha )+split.P5m( M_R+cr, alpha )))/2)*(pL+pR) + \
                      (1 - chi_zeta) * (split.P5p( M_L+cr, alpha )+split.P5m( M_R+cr, alpha )-1 )*(pL+pR)/2
        p_half_eta = (pB+pT)/2 + (((split.P5p( M_D+cr, alpha )+split.P5m( M_U+cr, alpha )))/2)*(pB+pT) + \
                      (1 - chi_eta) * (split.P5p( M_D+cr, alpha )+split.P5m( M_U+cr, alpha )-1 )*(pB+pT)/2

        # update phi vector components
        PhiT[:,:,1] = QT[:,:,1]/QT[:,:,0]
        PhiT[:,:,2] = QT[:,:,2]/QT[:,:,0]
        PhiT[:,:,3] = QT[:,:,3]/QT[:,:,0] + pT/QT[:,:,0]
        
        PhiB[:,:,1] = QB[:,:,1]/QB[:,:,0]
        PhiB[:,:,2] = QB[:,:,2]/QB[:,:,0]
        PhiB[:,:,3] = QB[:,:,3]/QB[:,:,0] + pB/QB[:,:,0]
        
        PhiL[:,:,1] = QL[:,:,1]/QL[:,:,0]
        PhiL[:,:,2] = QL[:,:,2]/QL[:,:,0]
        PhiL[:,:,3] = QL[:,:,3]/QL[:,:,0] + pL/QL[:,:,0]
        
        PhiR[:,:,1] = QR[:,:,1]/QR[:,:,0]
        PhiR[:,:,2] = QR[:,:,2]/QR[:,:,0]
        PhiR[:,:,3] = QR[:,:,3]/QR[:,:,0] + pR/QR[:,:,0]
        
        # pressure flux vector
        P_zeta[:,:,1] = p_half_zeta * mesh.s_proj[0:-1,:,0] / mesh.s_proj[0:-1,:,4]
        P_zeta[:,:,2] = p_half_zeta * mesh.s_proj[0:-1,:,1] / mesh.s_proj[0:-1,:,4]

        P_eta[:,:,1] = p_half_eta * mesh.s_proj[:,0:-1,2] / mesh.s_proj[:,0:-1,5]
        P_eta[:,:,2] = p_half_eta * mesh.s_proj[:,0:-1,3] / mesh.s_proj[:,0:-1,5]

        # flux vector reconstruction
        flux.face_flux( mdot_half_zeta, mdot_half_eta, Phi, P_zeta, P_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, domain.M, domain.N)

        # state.residual = -np.dstack( (state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1]) ) * \
        #                           ( ( E_hat_right * np.dstack([mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4]]) \
        #                             - E_hat_left * np.dstack([mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4]]) ) + \
        #                             ( F_hat_top * np.dstack([mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5]])\
        #                             - F_hat_bot * np.dstack([mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5]]) ) )

        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj, domain.M, domain.N ) 
        state.Q[1:-1,1:-1,:] = state.Qn[1:-1,1:-1,:] + state.residual / mesh.dV4

        # L_inf-norm residual
        state.res[state.n-1,0] = np.log10( np.max(state.residual[:,:,0] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,1] = np.log10( np.max(state.residual[:,:,1] * mesh.dV[1:-1,1:-1]) ) 
        state.res[state.n-1,2] = np.log10( np.max(state.residual[:,:,2] * mesh.dV[1:-1,1:-1]) )
        state.res[state.n-1,3] = np.log10( np.max(state.residual[:,:,3] * mesh.dV[1:-1,1:-1]) ) 

        #state.res[n-1] = np.log10( np.max(state.residual * mesh.dV4) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv) )
        state.T = state.p / (gas.R_fn(gas.Cp, gas.Cv) * state.Q[:,:,0])

        # update covariant velocities
        #state = covariant(mesh, state)
        soln_vars.calc_covariant(mesh.s_proj, state.u, state.v, state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, boundary, parameters, state, gas)

        # print iteration output
        if n % 10 == 0:
            print('Iteration: ' + str(state.n) + '    ' + str(round(state.res[n-1,0],3)) + '    ' + str(round(state.res[n-1,1],3)) + \
                                                 '    ' + str(round(state.res[n-1,2],3)) + '    ' + str(round(state.res[n-1,3],3)) )
            t.toc('Iteration time:')

    print('________________________________________________________________________________________________________________________________________')

    # post processing variables
    state = calc_postvars(state, gas)

    return state



# calculate viscous stress components 
def calc_stress( mesh, state, gas, face ):

    if face == '-':
        u_zeta = state.u[1:-1,1:-1] - state.u[0:-2,1:-1]
        u_eta = state.u[1:-1,1:-1] - state.u[1:-1,0:-2]
        v_zeta = state.v[1:-1,1:-1] - state.v[0:-2,1:-1]
        v_eta = state.v[1:-1,1:-1] - state.v[1:-1,0:-2]
        S_zetax = mesh.s_proj[0:-2,1:-1,0]
        S_etax = mesh.s_proj[1:-1,0:-2,2]
        S_zetay = mesh.s_proj[0:-2,1:-1,1]
        S_etay = mesh.s_proj[1:-1,0:-2,3]
    elif face == '+':
        u_zeta = state.u[2:,1:-1] - state.u[1:-1,1:-1]
        u_eta = state.u[1:-1,2:] - state.u[1:-1,1:-1]
        v_zeta = state.v[2:,1:-1] - state.v[1:-1,1:-1]
        v_eta = state.v[1:-1,2:] - state.v[1:-1,1:-1]
        S_zetax = mesh.s_proj[1:-1,1:-1,0]
        S_etax = mesh.s_proj[1:-1,1:-1,2]
        S_zetay = mesh.s_proj[1:-1,1:-1,1]
        S_etay = mesh.s_proj[1:-1,1:-1,3]

    Txx = (2/3)*gas.mu_fn(state.T[1:-1,1:-1]) * ( 2*(S_zetax*u_zeta + S_etax*u_eta) - S_zetay*v_zeta + S_etay*v_eta )
    Tyy = (2/3)*gas.mu_fn(state.T[1:-1,1:-1]) * ( 2*(S_zetay*v_zeta + S_etay*v_eta) - S_zetax*u_zeta + S_etax*u_eta )
    Txy = gas.mu_fn(state.T[1:-1,1:-1]) * (S_zetay*u_zeta + S_etay*u_eta + S_zetax*v_zeta + S_etax*v_eta)

    return Txx, Tyy, Txy


def calc_postvars(state, gas):

    state.Mach = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 ) / \
                           thermo.calc_c( state.p, state.Q[:,:,0], gas.gamma_fn(gas.Cp, gas.Cv) )
    state.vel = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 )
    state.isenRatio = (1+((gas.gamma_fn(gas.Cp, gas.Cv)-1)/2)*state.Mach**2)
    state.p0Ratio = state.isenRatio**(gas.gamma_fn(gas.Cp, gas.Cv)/(gas.gamma_fn(gas.Cp, gas.Cv)-1))
    state.p0 = state.p0Ratio * state.p
    state.T0 = state.isenRatio * state.T

    return state
