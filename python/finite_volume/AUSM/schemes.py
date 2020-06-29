# module importation
import numpy as np
from python.finite_volume.helper import thermo, split
from python.boundary.boundary_cond import enforce_bc, covariant
from python.finite_volume.timestepping import local_timestep
import python.finite_volume.soln_vars as soln_vars
import python.finite_volume.AUSM.flux as flux

# AUSM flux vector splitting scheme 
def AUSM( domain, mesh, parameters, state, gas ):

    n = 0

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
    state.res = np.ones( [parameters.iterations + 1, 4] )

    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])


    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        n = n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma) / state.Q[:,:,0] + state.p/state.Q[:,:,0]

        # local timestepping
        state = local_timestep( mesh, state, parameters, gas )

        # density at cell interfaces, upwinded
        rho_half_zeta = ( state.Q[0:-1,:,0] + state.Q[1:,:,0] ) / 2
        rho_half_eta =  ( state.Q[:,0:-1,0] + state.Q[:,1:,0] ) / 2

        # speed of sound at cell interfaces
        # from Liou 2006 (JCP 214)

        c_st = thermo.calc_c_star( state.ht, gas.gamma )

        c_L = c_st[0:-1,:] / np.maximum( np.sqrt(c_st[0:-1,:]), state.U[0:-1,:] )
        c_R = c_st[1:,:] / np.maximum( np.sqrt(c_st[1:,:]), -state.U[1:,:] ) 
        c_D = c_st[:,0:-1] / np.maximum( np.sqrt(c_st[:,0:-1]), state.V[:,0:-1] )
        c_U = c_st[:,1:] / np.maximum( np.sqrt(c_st[:,1:]), -state.V[:,1:] )

        c_half_zeta = np.minimum( c_L, c_R )
        c_half_eta =  np.minimum( c_D, c_U )

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
        state.res[n-1,0] = np.log10( np.max(state.residual[:,:,0] * mesh.dV[1:-1,1:-1]) ) 
        state.res[n-1,1] = np.log10( np.max(state.residual[:,:,1] * mesh.dV[1:-1,1:-1]) ) 
        state.res[n-1,2] = np.log10( np.max(state.residual[:,:,2] * mesh.dV[1:-1,1:-1]) )
        state.res[n-1,3] = np.log10( np.max(state.residual[:,:,3] * mesh.dV[1:-1,1:-1]) ) 

        #state.res[n-1] = np.log10( np.max(state.residual * mesh.dV4) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], state.u, state.v, gas.gamma )
        state.T = state.p / (gas.R * state.Q[:,:,0])

        # update covariant velocities
        #state = covariant(mesh, state)
        soln_vars.calc_covariant(mesh.s_proj, state.u, state.v, state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, parameters, state, gas)

    # post processing variables
    state.Mach = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 ) / thermo.calc_c( state.p, state.Q[:,:,0], gas.gamma )
    state.p0 = (1+((gas.gamma-1)/2)*state.Mach**2)**(gas.gamma/(gas.gamma-1)) * state.p
    state.n = n

    return state


# AUSM+up flux vector splitting scheme
def AUSMplusup( domain, mesh, parameters, state, gas ):

    n = 0

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
    state.res = np.ones( [parameters.iterations + 1, 4] )

    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])


    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        n = n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma) / state.Q[:,:,0] + state.p/state.Q[:,:,0]

        # local timestepping
        state = local_timestep( mesh, state, parameters, gas )

        # density at cell interfaces, upwinded
        rho_half_zeta = ( state.Q[0:-1,:,0] + state.Q[1:,:,0] ) / 2
        rho_half_eta =  ( state.Q[:,0:-1,0] + state.Q[:,1:,0] ) / 2

        # speed of sound at cell interfaces
        # from Liou 2006 (JCP 214)

        c_st = thermo.calc_c_star( state.ht, gas.gamma )

        c_L = c_st[0:-1,:] / np.maximum( np.sqrt(c_st[0:-1,:]), state.U[0:-1,:] )
        c_R = c_st[1:,:] / np.maximum( np.sqrt(c_st[1:,:]), -state.U[1:,:] ) 
        c_D = c_st[:,0:-1] / np.maximum( np.sqrt(c_st[:,0:-1]), state.V[:,0:-1] )
        c_U = c_st[:,1:] / np.maximum( np.sqrt(c_st[:,1:]), -state.V[:,1:] )

        c_half_zeta = np.minimum( c_L, c_R )
        c_half_eta =  np.minimum( c_D, c_U )

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

        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj[1:-1,1:-1,:], domain.M, domain.N ) 
        state.Q[1:-1,1:-1,:] = state.Qn[1:-1,1:-1,:] + state.residual / mesh.dV4

        # L_inf-norm residual
        state.res[n-1,0] = np.log10( np.max(state.residual[:,:,0] * mesh.dV[1:-1,1:-1]) ) 
        state.res[n-1,1] = np.log10( np.max(state.residual[:,:,1] * mesh.dV[1:-1,1:-1]) ) 
        state.res[n-1,2] = np.log10( np.max(state.residual[:,:,2] * mesh.dV[1:-1,1:-1]) )
        state.res[n-1,3] = np.log10( np.max(state.residual[:,:,3] * mesh.dV[1:-1,1:-1]) ) 

        #state.res[n-1] = np.log10( np.max(state.residual * mesh.dV4) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], state.u, state.v, gas.gamma )
        state.T = state.p / (gas.R * state.Q[:,:,0])

        # update covariant velocities
        #state = covariant(mesh, state)
        soln_vars.calc_covariant(mesh.s_proj, state.u, state.v, state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, parameters, state, gas)

    # post processing variables
    state.Mach = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 ) / thermo.calc_c( state.p, state.Q[:,:,0], gas.gamma )
    state.p0 = (1+((gas.gamma-1)/2)*state.Mach**2)**(gas.gamma/(gas.gamma-1)) * state.p
    state.n = n

    return state


# AUSMDV flux vector splitting scheme 
def AUSMDV( domain, mesh, parameters, state, gas ):

    n = 0

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
    state.res = np.ones( [parameters.iterations + 1, 4] )

    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])


    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        n = n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma) / state.Q[:,:,0] + state.p/state.Q[:,:,0]

        # local timestepping
        state = local_timestep( mesh, state, parameters, gas )

        # density at cell interfaces, upwinded
        rho_half_zeta = ( state.Q[0:-1,:,0] + state.Q[1:,:,0] ) / 2
        rho_half_eta =  ( state.Q[:,0:-1,0] + state.Q[:,1:,0] ) / 2

        # speed of sound at cell interfaces
        # from Liou 2006 (JCP 214)

        c_st = thermo.calc_c_star( state.ht, gas.gamma )

        c_L = c_st[0:-1,:] / np.maximum( np.sqrt(c_st[0:-1,:]), state.U[0:-1,:] )
        c_R = c_st[1:,:] / np.maximum( np.sqrt(c_st[1:,:]), -state.U[1:,:] ) 
        c_D = c_st[:,0:-1] / np.maximum( np.sqrt(c_st[:,0:-1]), state.V[:,0:-1] )
        c_U = c_st[:,1:] / np.maximum( np.sqrt(c_st[:,1:]), -state.V[:,1:] )

        c_half_zeta = np.minimum( c_L, c_R )
        c_half_eta =  np.minimum( c_D, c_U )

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

        #tic()

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
        state.res[n-1,0] = np.log10( np.max(state.residual[:,:,0] * mesh.dV[1:-1,1:-1]) ) 
        state.res[n-1,1] = np.log10( np.max(state.residual[:,:,1] * mesh.dV[1:-1,1:-1]) ) 
        state.res[n-1,2] = np.log10( np.max(state.residual[:,:,2] * mesh.dV[1:-1,1:-1]) )
        state.res[n-1,3] = np.log10( np.max(state.residual[:,:,3] * mesh.dV[1:-1,1:-1]) ) 

        #state.res[n-1] = np.log10( np.max(state.residual * mesh.dV4) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], state.u, state.v, gas.gamma )
        state.T = state.p / (gas.R * state.Q[:,:,0])

        # update covariant velocities
        #state = covariant(mesh, state)
        soln_vars.calc_covariant(mesh.s_proj, state.u, state.v, state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, parameters, state, gas)

    # post processing variables
    state.Mach = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 ) / thermo.calc_c( state.p, state.Q[:,:,0], gas.gamma )
    state.p0 = (1+((gas.gamma-1)/2)*state.Mach**2)**(gas.gamma/(gas.gamma-1)) * state.p
    state.n = n
    
    return state

