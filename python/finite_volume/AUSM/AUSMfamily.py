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
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / \
                                      state.Q[:,:,0] + state.p/state.Q[:,:,0]

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

        # dissipative term (Liou_JCP_160_2000)
        #Dm_zeta = np.abs( mdot_half_zeta )
        #Dm_eta  = np.abs( mdot_half_eta )

        # flux vector reconstruction
        # E_hat_left = (1/2) * mdot_half_zeta[0:-1,1:-1] * ( Phi[0:-2,1:-1,:] + Phi[1:-1,1:-1,:] ) \
        #             -(1/2) * Dm_zeta[0:-1,1:-1] * ( Phi[1:-1,1:-1,:] - Phi[0:-2,1:-1,:] ) \
        #                    + P_zeta[0:-1,1:-1,:]
        # E_hat_right= (1/2) * mdot_half_zeta[1:,1:-1] * ( Phi[1:-1,1:-1,:] + Phi[2:,1:-1,:] ) \
        #             -(1/2) * Dm_zeta[1:,1:-1] * ( Phi[2:,1:-1,:] - Phi[1:-1,1:-1,:] ) \
        #                    + P_zeta[1:,1:-1,:]
        # F_hat_bot =  (1/2) * mdot_half_eta[1:-1,0:-1] * ( Phi[1:-1,0:-2,:] + Phi[1:-1,1:-1,:] ) \
        #             -(1/2) * Dm_eta[1:-1,0:-1] * ( Phi[1:-1,1:-1,:] - Phi[1:-1,0:-2,:] ) \
        #                    + P_eta[1:-1,0:-1,:]
        # F_hat_top =  (1/2) * mdot_half_eta[1:-1,1:] * ( Phi[1:-1,1:-1,:] + Phi[1:-1,2:,:] ) \
        #             -(1/2) * Dm_eta[1:-1,1:] * ( Phi[1:-1,2:,:] - Phi[1:-1,1:-1,:] ) \
        #                    + P_eta[1:-1,1:,:]
        flux.face_flux( mdot_half_zeta, mdot_half_eta, Phi, P_zeta, P_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, domain.M, domain.N)

        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        # state.residual = -np.dstack( (state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1]) ) * \
        #                           ( ( E_hat_right * np.dstack([mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4]]) \
        #                             - E_hat_left * np.dstack([mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4]]) ) + \
        #                             ( F_hat_top * np.dstack([mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5]])\
        #                             - F_hat_bot * np.dstack([mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5]]) ) )
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj[1:-1,1:-1,:], domain.M, domain.N ) 
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
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / state.Q[:,:,0] + state.p/state.Q[:,:,0]

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
        M_co_eta = 0.5

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

        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj[1:-1,1:-1,:], domain.M, domain.N ) 
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
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / state.Q[:,:,0] + state.p/state.Q[:,:,0]

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
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj[1:-1,1:-1,:], domain.M, domain.N ) 
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
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / state.Q[:,:,0] + state.p/state.Q[:,:,0]

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

        # c_st = thermo.calc_c_star( state.ht, gas.gamma_fn(gas.Cp, gas.Cv) )

        # c_L = c_st[0:-1,:] / np.maximum( np.sqrt(c_st[0:-1,:]), state.U[0:-1,:] )
        # c_R = c_st[1:,:] / np.maximum( np.sqrt(c_st[1:,:]), -state.U[1:,:] ) 
        # c_D = c_st[:,0:-1] / np.maximum( np.sqrt(c_st[:,0:-1]), state.V[:,0:-1] )
        # c_U = c_st[:,1:] / np.maximum( np.sqrt(c_st[:,1:]), -state.V[:,1:] )

        # c_half_zeta =  np.minimum(c_L, c_R)
        # c_half_eta =   np.minimum(c_D, c_U)

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

        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj[1:-1,1:-1,:], domain.M, domain.N ) 
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

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / \
                                      state.Q[:,:,0] + state.p/state.Q[:,:,0]

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

        c_half_zeta = np.minimum( c_L, c_R )
        c_half_eta =  np.minimum( c_B, c_T )

        # cell face Mach numbers
        M_L = UL / c_half_zeta
        M_R = UR / c_half_zeta
        M_B = VB / c_half_eta
        M_U = VT / c_half_eta

        # split interface Mach numbers in the zeta and eta directions
        M_half_zeta = split.Mvp( M_L ) + split.Mvm( M_R )
        M_half_eta = split.Mvp( M_B ) + split.Mvm( M_U )

        # calculate mass flux at cell interfaces
        mdot_half_zeta = c_half_zeta * M_half_zeta * ( np.double(M_half_zeta>0) * QL[:,:,0] + np.double(M_half_zeta<=0) * QR[:,:,0] )
        mdot_half_eta =  c_half_eta  * M_half_eta  * ( np.double(M_half_eta>0)  * QB[:,:,0] + np.double(M_half_eta<=0)  * QT[:,:,0] )

        cr = 1e-60
        # calculate pressure flux at cell interfaces
        p_half_zeta = split.P1p( M_L+cr )*pL + split.P1m( M_R+cr )*pR
        p_half_eta =  split.P1p( M_B+cr )*pB + split.P1m( M_U+cr )*pT

        # initialize Phi vector components

        PhiT[:,:,1] = QT[:,:,1]/QT[:,:,0]
        PhiT[:,:,2] = QT[:,:,2]/QT[:,:,0]
        PhiT[:,:,3] = QT[:,:,3]/QT[:,:,0]
        
        PhiB[:,:,1] = QB[:,:,1]/QB[:,:,0]
        PhiB[:,:,2] = QB[:,:,2]/QB[:,:,0]
        PhiB[:,:,3] = QB[:,:,3]/QB[:,:,0]
        
        PhiL[:,:,1] = QL[:,:,1]/QL[:,:,0]
        PhiL[:,:,2] = QL[:,:,2]/QL[:,:,0]
        PhiL[:,:,3] = QL[:,:,3]/QL[:,:,0]
        
        PhiR[:,:,1] = QR[:,:,1]/QR[:,:,0]
        PhiR[:,:,2] = QR[:,:,2]/QR[:,:,0]
        PhiR[:,:,3] = QR[:,:,3]/QR[:,:,0]
        
        # pressure flux vector
        P_zeta[:,:,1] = p_half_zeta * mesh.s_proj[0:-1,:,0] / mesh.s_proj[0:-1,:,4]
        P_zeta[:,:,2] = p_half_zeta * mesh.s_proj[0:-1,:,1] / mesh.s_proj[0:-1,:,4]

        P_eta[:,:,1] = p_half_eta * mesh.s_proj[:,0:-1,2] / mesh.s_proj[:,0:-1,5]
        P_eta[:,:,2] = p_half_eta * mesh.s_proj[:,0:-1,3] / mesh.s_proj[:,0:-1,5]

        # compute flux at each cell face
        flux.face_flux_muscl( mdot_half_zeta, mdot_half_eta, PhiL, PhiR, PhiB, PhiT, P_zeta, P_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, domain.M, domain.N)

        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj[1:-1,1:-1,:], domain.M, domain.N ) 
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

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / \
                                      state.Q[:,:,0] + state.p/state.Q[:,:,0]

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
        PhiT[:,:,3] = QT[:,:,3]/QT[:,:,0]
        
        PhiB[:,:,1] = QB[:,:,1]/QB[:,:,0]
        PhiB[:,:,2] = QB[:,:,2]/QB[:,:,0]
        PhiB[:,:,3] = QB[:,:,3]/QB[:,:,0]
        
        PhiL[:,:,1] = QL[:,:,1]/QL[:,:,0]
        PhiL[:,:,2] = QL[:,:,2]/QL[:,:,0]
        PhiL[:,:,3] = QL[:,:,3]/QL[:,:,0]
        
        PhiR[:,:,1] = QR[:,:,1]/QR[:,:,0]
        PhiR[:,:,2] = QR[:,:,2]/QR[:,:,0]
        PhiR[:,:,3] = QR[:,:,3]/QR[:,:,0]
        
        # pressure flux vector
        P_zeta[:,:,1] = p_half_zeta * mesh.s_proj[0:-1,:,0] / mesh.s_proj[0:-1,:,4]
        P_zeta[:,:,2] = p_half_zeta * mesh.s_proj[0:-1,:,1] / mesh.s_proj[0:-1,:,4]

        P_eta[:,:,1] = p_half_eta * mesh.s_proj[:,0:-1,2] / mesh.s_proj[:,0:-1,5]
        P_eta[:,:,2] = p_half_eta * mesh.s_proj[:,0:-1,3] / mesh.s_proj[:,0:-1,5]


        flux.face_flux_muscl( mdot_half_zeta, mdot_half_eta, PhiL, PhiR, PhiB, PhiT, P_zeta, P_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, domain.M, domain.N)

        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj[1:-1,1:-1,:], domain.M, domain.N ) 
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


def calc_postvars(state, gas):

    state.Mach = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 ) / \
                           thermo.calc_c( state.p, state.Q[:,:,0], gas.gamma_fn(gas.Cp, gas.Cv) )
    state.vel = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 )
    state.p0 = (1+((gas.gamma_fn(gas.Cp, gas.Cv)-1)/2)*state.Mach**2)** \
                   (gas.gamma_fn(gas.Cp, gas.Cv)/(gas.gamma_fn(gas.Cp, gas.Cv)-1)) * state.p
    state.T0 = (1+((gas.gamma_fn(gas.Cp, gas.Cv)-1)/2)*state.Mach**2) * state.T

    return state
