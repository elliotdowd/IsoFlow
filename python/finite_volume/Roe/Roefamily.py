# module importation
import numpy as np
from python.finite_volume.helper import thermo, split
from python.boundary.boundary_cond import enforce_bc, covariant
from python.finite_volume.timestepping import local_timestep
import python.finite_volume.soln_vars as soln_vars
import python.finite_volume.flux as flux
import python.finite_volume.muscl as muscl
from pytictoc import TicToc

# Roe finite difference scheme, work in progress
def RoeFDS( domain, mesh, boundary, parameters, state, gas ):

    print('Roe Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # initialize phi and p state vectors
    Phi = np.zeros( (domain.M+2, domain.N+2, 4) )
    Psi = np.zeros( (domain.M+2, domain.N+2, 4) )

    P_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    P_zeta[:,:,0] = np.zeros( (domain.M+1, domain.N+2) )
    P_zeta[:,:,3] = np.zeros( (domain.M+1, domain.N+2) )

    P_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    P_eta[:,:,0] = np.zeros( (domain.M+2, domain.N+1) )
    P_eta[:,:,3] = np.zeros( (domain.M+2, domain.N+1) )

    # initialize Psi vector components
    Psi_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    Psi_zeta[:,:,1] = mesh.s_proj[1:,:,0]/mesh.s_proj[1:,:,4]
    Psi_zeta[:,:,2] = mesh.s_proj[1:,:,1]/mesh.s_proj[1:,:,4]

    Psi_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    Psi_eta[:,:,1] = mesh.s_proj[:,1:,2]/mesh.s_proj[:,1:,5]
    Psi_eta[:,:,2] = mesh.s_proj[:,1:,3]/mesh.s_proj[:,1:,5]

    # normal vectors
    N_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    N_zeta[:,:,1] = mesh.s_proj[0:-1,:,0]/mesh.s_proj[0:-1,:,4]
    N_zeta[:,:,2] = mesh.s_proj[0:-1,:,1]/mesh.s_proj[0:-1,:,4]

    N_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    N_eta[:,:,1] = mesh.s_proj[:,0:-1,2]/mesh.s_proj[:,0:-1,5]
    N_eta[:,:,2] = mesh.s_proj[:,0:-1,3]/mesh.s_proj[:,0:-1,5]

    # flux jacobian
    A = np.zeros( (domain.M+2, domain.N+2, 4, 4), dtype='float', order='F' )

    E_hat_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    E_hat_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )

    Ec_half = np.zeros( (domain.M+1, domain.N+2, 4), dtype='float', order='F' )
    Fc_half = np.zeros( (domain.M+2, domain.N+1, 4), dtype='float', order='F' )

    Ed_half = np.zeros( (domain.M+1, domain.N+2, 4), dtype='float', order='F' )
    Fd_half = np.zeros( (domain.M+2, domain.N+1, 4), dtype='float', order='F' )

    state.residual = np.zeros( [domain.M, domain.N, 4], dtype='float', order='F' )
    state.res = np.ones( [parameters.iterations + 1, 4] )

    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])

    t = TicToc()

    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        t.tic()

        n = n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # local timestepping
        state = local_timestep( domain, mesh, state, parameters, gas )

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / \
                                      state.Q[:,:,0] + state.p/state.Q[:,:,0]

        # Roe averaged variables
        rho_half_zeta = np.sqrt(state.Q[0:-1,:,0]*state.Q[1:,:,0])
        rho_half_eta = np.sqrt(state.Q[:,0:-1,0]*state.Q[:,1:,0])

        # QL, QR, QB, QT = muscl.MUSCL( state.Q, parameters.epsilon, parameters.kappa, parameters.limiter )

        # QL = state.Q[0:-2,1:-1,:]
        # QR = state.Q[1:-1,1:-1,:]
        # QB = state.Q[1:-1,0:-2,:]
        # QT = state.Q[1:-1,1:-1,:]

        # cell normal vectors
        nx = mesh.s_proj[1:,:,0]/mesh.s_proj[1:,:,4]
        t.tic()

        # calculate flux jacobian matrix
        for i in range(0, domain.M+2):

            for j in range(0, domain.N+2):

                A[i,j,:,:] = flux_jacobian( state.u[i,j], state.v[i,j], \
                                            mesh.s_proj[i,j,0]/mesh.s_proj[i,j,4], mesh.s_proj[i,j,1]/mesh.s_proj[i,j,4], \
                                            gas.gamma_fn(gas.Cp[i,j], gas.Cv[i,j]), state.Q[i,j,3]/state.Q[i,j,0] )

                A[i,j,:,:] = A[i,j,:,:] / mesh.dV[i,j]
                

        for i in range(1, domain.M+1):

            for j in range(1, domain.N+1):

                ih = i-1
                jh = j-1

                E_hat_left[ih,jh,:] =  np.matmul( A[i-1,j,:,:], state.Q[i-1,j,:] )
                E_hat_right[ih,jh,:] = np.matmul( A[i,j,:,:], state.Q[i,j,:] )
                F_hat_bot[ih,jh,:]  =  np.matmul( A[i,j-1,:,:], state.Q[i,j-1,:] )
                F_hat_top[ih,jh,:]  =  np.matmul( A[i,j,:,:], state.Q[i,j,:] )


        t.toc('Flux Jacobian time:')



        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        flux.residual( state.residual, state.dt[1:-1, 1:-1] * mesh.dV[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
                          mesh.s_proj[1:-1,1:-1,:], domain.M, domain.N )
        state.Q[1:-1,1:-1,:] = state.Qn[1:-1,1:-1,:] + state.residual / mesh.dV4

        # L_inf-norm residual
        state.res[n-1,0] = np.log10( np.max(state.residual[:,:,0] * mesh.dV[1:-1,1:-1]) ) 
        state.res[n-1,1] = np.log10( np.max(state.residual[:,:,1] * mesh.dV[1:-1,1:-1]) ) 
        state.res[n-1,2] = np.log10( np.max(state.residual[:,:,2] * mesh.dV[1:-1,1:-1]) )
        state.res[n-1,3] = np.log10( np.max(state.residual[:,:,3] * mesh.dV[1:-1,1:-1]) ) 

        #state.res[n-1] = np.log10( np.max(state.residual * mesh.dV4) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv) )
        state.T = state.p / (gas.R_fn(gas.Cp, gas.Cv) * state.Q[:,:,0])

        # update covariant velocities
        soln_vars.calc_covariant(mesh.s_proj, state.u, state.v, state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, boundary, parameters, state, gas)

        # print iteration output
        if n % 10 == 0:
            print('Iteration: ' + str(n) + '    ' + str(round(state.res[n-1,0],3)) + '    ' + str(round(state.res[n-1,1],3)) + \
                                           '    ' + str(round(state.res[n-1,2],3)) + '    ' + str(round(state.res[n-1,3],3)) )
            t.toc('Iteration time:')

    print('________________________________________________________________________________________________________________________________________')

    # post processing variables
    state = calc_postvars(state, gas, n)

    return state


# Roe FV scheme
# see Ren, Gu, Li "New Roe Scheme For All Speed Flow" (2011)
def RoeFVS( domain, mesh, boundary, parameters, state, gas ):

    print('Roe Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # initialize phi and p state vectors
    Phi = np.zeros( (domain.M+2, domain.N+2, 4), dtype='float', order='F' )
    Psi = np.zeros( (domain.M+2, domain.N+2, 4), dtype='float', order='F' )

    P_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    P_zeta[:,:,0] = np.zeros( (domain.M+1, domain.N+2) )
    P_zeta[:,:,3] = np.zeros( (domain.M+1, domain.N+2) )

    P_eta = np.zeros( (domain.M+2, domain.N+1, 4), dtype='float', order='F' )
    P_eta[:,:,0] = np.zeros( (domain.M+2, domain.N+1) )
    P_eta[:,:,3] = np.zeros( (domain.M+2, domain.N+1) )

    # initialize Psi vector components
    Psi_zeta = np.zeros( (domain.M+1, domain.N+2, 4), dtype='float', order='F' )
    Psi_zeta[:,:,1] = mesh.s_proj[1:,:,0]/mesh.s_proj[1:,:,4]
    Psi_zeta[:,:,2] = mesh.s_proj[1:,:,1]/mesh.s_proj[1:,:,4]

    Psi_eta = np.zeros( (domain.M+2, domain.N+1, 4), dtype='float', order='F' )
    Psi_eta[:,:,1] = mesh.s_proj[:,1:,2]/mesh.s_proj[:,1:,5]
    Psi_eta[:,:,2] = mesh.s_proj[:,1:,3]/mesh.s_proj[:,1:,5]

    # normal vectors
    N_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    N_zeta[:,:,1] = mesh.s_proj[0:-1,:,0]/mesh.s_proj[0:-1,:,4]
    N_zeta[:,:,2] = mesh.s_proj[0:-1,:,1]/mesh.s_proj[0:-1,:,4]

    N_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    N_eta[:,:,1] = mesh.s_proj[:,0:-1,2]/mesh.s_proj[:,0:-1,5]
    N_eta[:,:,2] = mesh.s_proj[:,0:-1,3]/mesh.s_proj[:,0:-1,5]

    E_hat_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    E_hat_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )

    # convective terms

    Fbar_left = np.zeros( (domain.M+1, domain.N+2, 4), dtype='float', order='F' )
    Fbar_right = np.zeros( (domain.M+1, domain.N+2, 4), dtype='float', order='F' )
    Fbar_bot = np.zeros( (domain.M+2, domain.N+1, 4), dtype='float', order='F' )
    Fbar_top = np.zeros( (domain.M+2, domain.N+1, 4), dtype='float', order='F' )

    Ec_half = np.zeros( (domain.M+1, domain.N+2, 4), dtype='float', order='F' )
    Fc_half = np.zeros( (domain.M+2, domain.N+1, 4), dtype='float', order='F' )

    Ed_half = np.zeros( (domain.M+1, domain.N+2, 4), dtype='float', order='F' )
    Fd_half = np.zeros( (domain.M+2, domain.N+1, 4), dtype='float', order='F' )

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

        # initialize Phi vector components
        Phi[:,:,0] = state.Q[:,:,0]
        Phi[:,:,1] = state.Q[:,:,0]*state.u
        Phi[:,:,2] = state.Q[:,:,0]*state.v
        Phi[:,:,3] = state.Q[:,:,0]*state.ht

        # Roe averaged quantities
        rho_zeta = np.sqrt( state.Q[0:-1,:,0] * state.Q[1:,:,0] )
        rho_eta = np.sqrt( state.Q[:,0:-1,0] * state.Q[:,1:,0] )
        U_roe = roe_average( state.Q[0:-1,:,0], state.Q[1:,:,0], state.U[0:-1,:], state.U[1:,:] )
        V_roe = roe_average( state.Q[:,0:-1,0], state.Q[:,1:,0], state.U[:,0:-1], state.U[:,1:] )

        # [0, nx, ny, U]
        Psi_zeta[:,:,3] = U_roe
        Psi_eta[:,:,3]  = V_roe
        
        # calculate central term in eqch computational direction
        # Fbar_left = flux.roe_fbar( Fbar_left, Phi[0:-1,:,:], state.U[0:-1,:], N_zeta, state.p[0:-1,:], domain.M+1, domain.N+2 )
        # Fbar_right = flux.roe_fbar( Fbar_right, Phi[1:,:,:], state.U[1:,:], N_zeta, state.p[1:,:], domain.M+1, domain.N+2 )

        # Ec_half =  (1/2) * ( Fbar_left + Fbar_right )
        # Fc_half =  (1/2) * ( Fbar_bot + Fbar_top )

        Ec_half =  (1/2) * ( F_bar( Phi[0:-1,:,:], state.U[0:-1,:], N_zeta, state.p[0:-1,:], domain.M+1, domain.N+2 ) + \
                             F_bar( Phi[1:,:,:], state.U[1:,:], N_zeta, state.p[1:,:], domain.M+1, domain.N+2 ) )
        Fc_half =  (1/2) * ( F_bar( Phi[:,0:-1,:], state.V[:,0:-1], N_eta, state.p[:,0:-1], domain.M+2, domain.N+1 ) + \
                             F_bar( Phi[:,1:,:], state.V[:,1:], N_eta, state.p[:,1:], domain.M+2, domain.N+1 ) )

        # upwind dissipation
        D_zeta = state.eig[0,0,:,:]
        D_eta  = state.eig[1,0,:,:]

        # equations 13-16 in Ren et Al 2011.
        Pu_zeta = np.maximum( 0, state.c[1:,:] - np.abs(state.U[1:,:]) ) * rho_zeta * ( state.U[1:,:]-state.U[0:-1,:] )
        Pp_zeta = np.sign(state.U[1:,:]) * np.minimum( np.abs(state.U[1:,:]), state.c[1:,:] ) * ( state.p[1:,:]-state.p[0:-1,:] )/state.c[1:,:]
        Uu_zeta = np.sign(state.U[1:,:]) * np.minimum( np.abs(state.U[1:,:]), state.c[1:,:] ) * ( state.U[1:,:]-state.U[0:-1,:] )/state.c[1:,:]
        Up_zeta = np.maximum( 0, state.c[1:,:]-np.abs(state.U[1:,:])) * ( state.p[1:,:]-state.p[0:-1,:] ) / (rho_zeta*state.c[1:,:]**2)

        Pu_eta  = np.maximum( 0, state.c[:,1:] - np.abs(state.V[:,1:]) ) * rho_eta * ( state.V[:,1:]-state.V[:,0:-1] )
        Pp_eta  = np.sign(state.V[:,1:]) * np.minimum( np.abs(state.V[:,1:]), state.c[:,1:] ) * ( state.p[:,1:]-state.p[:,0:-1] )/state.c[:,1:]
        Uu_eta  = np.sign(state.V[:,1:]) * np.minimum( np.abs(state.V[:,1:]), state.c[:,1:] ) * ( state.V[:,1:]-state.V[:,0:-1] )/state.c[:,1:]
        Up_eta  = np.maximum( 0, state.c[:,1:]-np.abs(state.V[:,1:]) ) * ( state.p[:,1:]-state.p[:,0:-1] ) / (rho_eta*state.c[:,1:]**2)

        # numerical dissipation terms
        flux.roe_dissipation( Ed_half, Fd_half, state.p, Phi, Psi_zeta, Psi_eta, D_zeta, D_eta, Pu_zeta, Pp_zeta, Uu_zeta, Up_zeta, \
                                                              Pu_eta, Pp_eta, Uu_eta, Up_eta, domain.M, domain.N )

        # for i in range(0, 4):
        #     if i < 3:
        #         Ed_half[:,:,i] = -(1/2) * ( D_zeta[1:,:] * (Phi[1:,:,i]-Phi[0:-1,:,i]) + (Pu_zeta+Pp_zeta)*Psi_zeta[:,:,i] + \
        #                                                                                 (Uu_zeta+Up_zeta)*Phi[1:,:,i] )
        #         Fd_half[:,:,i] = -(1/2) * ( D_eta[:,1:] * (Phi[:,1:,i]-Phi[:,0:-1,i]) +  (Pu_eta+Pp_eta)*Psi_eta[:,:,i] + \
        #                                                                                 (Uu_eta+Up_eta)*Phi[:,1:,i] )
        #     else:
        #         Ed_half[:,:,i] = -(1/2) * ( D_zeta[1:,:] * ((Phi[1:,:,i]-state.p[1:,:])-(Phi[0:-1,:,i]-state.p[0:-1,:])) + (Pu_zeta+Pp_zeta)*Psi_zeta[:,:,i] + \
        #                                                                                                                    (Uu_zeta+Up_zeta)*Phi[1:,:,i] )
        #         Fd_half[:,:,i] = -(1/2) * ( D_eta[:,1:] * ((Phi[:,1:,i]-state.p[:,1:])-(Phi[:,0:-1,i]-state.p[:,0:-1])) +  (Pu_eta+Pp_eta)*Psi_eta[:,:,i] + \
        #                                                                                 (Uu_eta+Up_eta)*Phi[:,1:,i] )

        # flux summation
        E_hat_left  = Ec_half[0:-1,1:-1,:] + Ed_half[0:-1,1:-1,:]
        E_hat_right = Ec_half[1:,1:-1,:] + Ed_half[1:,1:-1,:]
        F_hat_bot   = Fc_half[1:-1,0:-1,:] + Fd_half[1:-1,0:-1,:]
        F_hat_top   = Fc_half[1:-1,1:,:] + Fd_half[1:-1,1:,:]

        state.residual = -np.dstack( (state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1], state.dt[1:-1, 1:-1]) ) * \
                                  ( ( E_hat_right * np.dstack([mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4], mesh.s_proj[1:-1,1:-1,4]]) \
                                    - E_hat_left * np.dstack([mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4], mesh.s_proj[0:-2,1:-1,4]]) ) + \
                                    ( F_hat_top * np.dstack([mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5], mesh.s_proj[1:-1,1:-1,5]])\
                                    - F_hat_bot * np.dstack([mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5], mesh.s_proj[1:-1,0:-2,5]]) ) )

        # update residuals and state vector at each interior cell, from Fortran 90 subroutine
        # flux.residual( state.residual, state.dt[1:-1, 1:-1], E_hat_left, E_hat_right, F_hat_bot, F_hat_top,\
        #                   mesh.s_proj[1:-1,1:-1,:], domain.M, domain.N ) 
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


# Improved Roe FV scheme
# see Ren, Gu, Li "New Roe Scheme For All Speed Flow" (2011)
def RoeFVSimproved( domain, mesh, boundary, parameters, state, gas ):

    print('Roe Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # initialize phi and p state vectors
    Phi = np.zeros( (domain.M+2, domain.N+2, 4) )
    Psi = np.zeros( (domain.M+2, domain.N+2, 4) )

    P_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    P_zeta[:,:,0] = np.zeros( (domain.M+1, domain.N+2) )
    P_zeta[:,:,3] = np.zeros( (domain.M+1, domain.N+2) )

    P_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    P_eta[:,:,0] = np.zeros( (domain.M+2, domain.N+1) )
    P_eta[:,:,3] = np.zeros( (domain.M+2, domain.N+1) )

    # initialize Psi vector components
    Psi_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    Psi_zeta[:,:,1] = mesh.s_proj[1:,:,0]/mesh.s_proj[1:,:,4]
    Psi_zeta[:,:,2] = mesh.s_proj[1:,:,1]/mesh.s_proj[1:,:,4]

    Psi_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    Psi_eta[:,:,1] = mesh.s_proj[:,1:,2]/mesh.s_proj[:,1:,5]
    Psi_eta[:,:,2] = mesh.s_proj[:,1:,3]/mesh.s_proj[:,1:,5]

    # normal vectors
    N_zeta = np.zeros( (domain.M+1, domain.N+2, 4) )
    N_zeta[:,:,1] = mesh.s_proj[0:-1,:,0]/mesh.s_proj[0:-1,:,4]
    N_zeta[:,:,2] = mesh.s_proj[0:-1,:,1]/mesh.s_proj[0:-1,:,4]

    N_eta = np.zeros( (domain.M+2, domain.N+1, 4) )
    N_eta[:,:,1] = mesh.s_proj[:,0:-1,2]/mesh.s_proj[:,0:-1,5]
    N_eta[:,:,2] = mesh.s_proj[:,0:-1,3]/mesh.s_proj[:,0:-1,5]

    E_hat_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    E_hat_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    F_hat_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )

    Ec_half = np.zeros( (domain.M+1, domain.N+2, 4), dtype='float', order='F' )
    Fc_half = np.zeros( (domain.M+2, domain.N+1, 4), dtype='float', order='F' )

    Ed_left = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    Ed_right = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    Fd_bot = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )
    Fd_top = np.zeros( (domain.M, domain.N, 4), dtype='float', order='F' )


    state.residual = np.zeros( [domain.M, domain.N, 4], dtype='float', order='F' )
    state.res = np.ones( [parameters.iterations + 1, 4] )

    mesh.dV4 = np.dstack([mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1], mesh.dV[1:-1,1:-1]])

    t = TicToc()

    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        t.tic()

        n = n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # local timestepping
        state = local_timestep( domain, mesh, state, parameters, gas )

        # simplify variable notation from state vector
        state.u = state.Q[:,:,1] / state.Q[:,:,0]
        state.v = state.Q[:,:,2] / state.Q[:,:,0]
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv)) / \
                                      state.Q[:,:,0] + state.p/state.Q[:,:,0]

        # initialize Phi vector components
        Phi[:,:,0] = state.Q[:,:,0]
        Phi[:,:,1] = state.Q[:,:,0]*state.u
        Phi[:,:,2] = state.Q[:,:,0]*state.v
        Phi[:,:,3] = state.Q[:,:,0]*state.ht

        # [0, nx, ny, U]
        Psi_zeta[:,:,3] = state.U[1:,:]
        Psi_eta[:,:,3]  = state.V[:,1:]

        # calculate central term in eqch computational direction
        Ec_half =  (1/2) * ( F_bar( Phi[0:-1,:,:], state.U[0:-1,:], N_zeta, state.p[0:-1,:], domain.M+1, domain.N+2 ) + \
                             F_bar( Phi[1:,:,:], state.U[1:,:], N_zeta, state.p[1:,:], domain.M+1, domain.N+2 ) )
        Fc_half =  (1/2) * ( F_bar( Phi[:,0:-1,:], state.V[:,0:-1], N_eta, state.p[:,0:-1], domain.M+2, domain.N+1 ) + \
                             F_bar( Phi[:,1:,:], state.V[:,1:], N_eta, state.p[:,1:], domain.M+2, domain.N+1 ) )

        # shock detector quantity
        b = shock_detect( state.p )

        # s1 and s2 coefficients
        M = np.sqrt( state.u**2 + state.v**2 ) / state.c
        s1 = 1 - f8(M[1:-1,1:-1])
        s2 = f8(b)

        # upwind dissipation
        D_zeta = state.eig[0,0,:,:]
        D_eta  = state.eig[1,0,:,:]

        # equations 13-16 in Ren et Al 2011.
        Pu_zeta = np.maximum( 0, state.c[1:,:] - np.abs(state.U[1:,:]) ) * state.Q[1:,:,0]*( state.U[1:,:]-state.U[0:-1,:] )
        Pp_zeta = np.sign(state.U[1:,:]) * np.minimum( np.abs(state.U[1:,:]), state.c[1:,:] ) * ( state.p[1:,:]-state.p[0:-1,:] )/state.c[1:,:]
        Uu_zeta = np.sign(state.U[1:,:]) * np.minimum( np.abs(state.U[1:,:]), state.c[1:,:] ) * ( state.U[1:,:]-state.U[0:-1,:] )/state.c[1:,:]

        Up_left = s1 * s2 * np.maximum( 0, state.c[0:-2,1:-1]-np.abs(state.U[0:-2,1:-1])) * (state.p[1:-1,1:-1]-state.p[0:-2,1:-1]) / \
                                                                                       (state.Q[0:-2,1:-1,0]*state.c[0:-2,1:-1]**2)
        Up_right = s1 * s2 * np.maximum( 0, state.c[1:-1,1:-1]-np.abs(state.U[1:-1,1:-1])) * (state.p[2:,1:-1]-state.p[1:-1,1:-1]) / \
                                                                                        (state.Q[1:-1,1:-1,0]*state.c[1:-1,1:-1]**2)

        Pu_eta  = np.maximum( 0, state.c[:,1:] - np.abs(state.V[:,1:]) ) * state.Q[:,1:,0]*( state.V[:,1:]-state.V[:,0:-1] )
        Pp_eta  = np.sign(state.V[:,1:]) * np.minimum( np.abs(state.V[:,1:]), state.c[:,1:] ) * ( state.p[:,1:]-state.p[:,0:-1] )/state.c[:,1:]
        Uu_eta  = np.sign(state.V[:,1:]) * np.minimum( np.abs(state.V[:,1:]), state.c[:,1:] ) * ( state.V[:,1:]-state.V[:,0:-1] )/state.c[:,1:]

        Up_bot  = s1 * s2 * np.maximum( 0, state.c[1:-1,0:-2]-np.abs(state.V[1:-1,0:-2]) ) * ( state.p[1:-1,1:-1]-state.p[1:-1,0:-2] ) / \
                                                                                           (state.Q[1:-1,0:-2,0]*state.c[1:-1,0:-2]**2)
        Up_top  = s1 * s2 * np.maximum( 0, state.c[1:-1,1:-1]-np.abs(state.V[1:-1,1:-1]) ) * ( state.p[1:-1,2:]-state.p[1:-1,1:-1] ) / \
                                                                                           (state.Q[1:-1,1:-1,0]*state.c[1:-1,1:-1]**2)

        # set Up values to zero for supersonic flow:
        mask = np.where(M[1:-1,1:-1]>1)
        Up_left[mask] = 0
        Up_right[mask] = 0
        Up_bot[mask] = 0
        Up_top[mask] = 0


        # numerical dissipation terms
        for i in range(0, 4):
            if i < 3:
                Ed_left[:,:,i] = -(1/2) * ( D_zeta[0:-2,1:-1] * (Phi[1:-1,1:-1,i]-Phi[0:-2,1:-1,i]) + (Pu_zeta[0:-1,1:-1]+Pp_zeta[0:-1,1:-1])*Psi_zeta[0:-1,1:-1,i] + \
                                                                                                      (Uu_zeta[0:-1,1:-1]+Up_left)*Phi[0:-2,1:-1,i] )
                Ed_right[:,:,i] = -(1/2) * ( D_zeta[1:-1,1:-1] * (Phi[2:,1:-1,i]-Phi[1:-1,1:-1,i]) + (Pu_zeta[1:,1:-1]+Pp_zeta[1:,1:-1])*Psi_zeta[1:,1:-1,i] + \
                                                                                                      (Uu_zeta[1:,1:-1]+Up_right)*Phi[1:-1,1:-1,i] )
                Fd_bot[:,:,i] = -(1/2) * ( D_eta[1:-1,0:-2] * (Phi[1:-1,1:-1,i]-Phi[1:-1,0:-2,i]) + (Pu_eta[1:-1,0:-1]+Pp_eta[1:-1,0:-1])*Psi_eta[1:-1,0:-1,i] + \
                                                                                                      (Uu_eta[1:-1,0:-1]+Up_bot)*Phi[1:-1,0:-2,i] )
                Fd_top[:,:,i] = -(1/2) * ( D_eta[1:-1,1:-1] * (Phi[1:-1,2:,i]-Phi[1:-1,1:-1,i]) + (Pu_eta[1:-1,1:]+Pp_eta[1:-1,1:])*Psi_eta[1:-1,1:,i] + \
                                                                                                      (Uu_eta[1:-1,1:]+Up_bot)*Phi[1:-1,1:-1,i] )

            else:
                Ed_left[:,:,i] = -(1/2) * ( D_zeta[0:-2,1:-1]  * (Phi[1:-1,1:-1,i]-state.p[1:-1,1:-1]-Phi[0:-2,1:-1,i]-state.p[0:-2,1:-1]) + \
                                                                 (Pu_zeta[0:-1,1:-1]+Pp_zeta[0:-1,1:-1])*Psi_zeta[0:-1,1:-1,i] + \
                                                                 (Uu_zeta[0:-1,1:-1]+Up_left)*Phi[0:-2,1:-1,i] )
                Ed_right[:,:,i] = -(1/2) * ( D_zeta[1:-1,1:-1] * (Phi[2:,1:-1,i]-state.p[2:,1:-1]-Phi[1:-1,1:-1,i]-state.p[1:-1,1:-1]) + \
                                                                 (Pu_zeta[1:,1:-1]+Pp_zeta[1:,1:-1])*Psi_zeta[1:,1:-1,i] + \
                                                                 (Uu_zeta[1:,1:-1]+Up_right)*Phi[1:-1,1:-1,i] )
                Fd_bot[:,:,i] = -(1/2) * ( D_eta[1:-1,0:-2] * (Phi[1:-1,1:-1,i]-state.p[1:-1,1:-1]-Phi[1:-1,0:-2,i]-state.p[1:-1,0:-2]) + \
                                                              (Pu_eta[1:-1,0:-1]+Pp_eta[1:-1,0:-1])*Psi_eta[1:-1,0:-1,i] + \
                                                              (Uu_eta[1:-1,0:-1]+Up_bot)*Phi[1:-1,0:-2,i] )
                Fd_top[:,:,i] = -(1/2) * ( D_eta[1:-1,1:-1] * (Phi[1:-1,2:,i]-state.p[1:-1,2:]-Phi[1:-1,1:-1,i]-state.p[1:-1,1:-1]) + \
                                                              (Pu_eta[1:-1,1:]+Pp_eta[1:-1,1:])*Psi_eta[1:-1,1:,i] + \
                                                              (Uu_eta[1:-1,1:]+Up_bot)*Phi[1:-1,1:-1,i] )


        # flux summation
        # E_hat_left  = Ec_half[0:-1,1:-1,:] + Ed_left[0:-1,1:-1,:]
        # E_hat_right = Ec_half[1:,1:-1,:] + Ed_right[1:,1:-1,:]
        # F_hat_bot   = Fc_half[1:-1,0:-1,:] + Fd_bot[1:-1,0:-1,:]
        # F_hat_top   = Fc_half[1:-1,1:,:] + Fd_top[1:-1,1:,:]
        E_hat_left  = Ec_half[0:-1,1:-1,:] + Ed_left
        E_hat_right = Ec_half[1:,1:-1,:] + Ed_right
        F_hat_bot   = Fc_half[1:-1,0:-1,:] + Fd_bot
        F_hat_top   = Fc_half[1:-1,1:,:] + Fd_top

        # flux vector reconstruction
        # flux.face_flux( mdot_half_zeta, mdot_half_eta, Phi, P_zeta, P_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, domain.M, domain.N)

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
        state.p = thermo.calc_p( state.Q[:,:,0], state.Q[:,:,3], state.u, state.v, gas.gamma_fn(gas.Cp, gas.Cv) )
        state.T = state.p / (gas.R_fn(gas.Cp, gas.Cv) * state.Q[:,:,0])

        # update covariant velocities
        soln_vars.calc_covariant(mesh.s_proj, state.u, state.v, state.U, state.V, domain.M+2, domain.N+2)

        # enforce boundary conditions
        state = enforce_bc(domain, mesh, boundary, parameters, state, gas)

        # print iteration output
        if n % 10 == 0:
            print('Iteration: ' + str(n) + '    ' + str(round(state.res[n-1,0],3)) + '    ' + str(round(state.res[n-1,1],3)) + \
                                           '    ' + str(round(state.res[n-1,2],3)) + '    ' + str(round(state.res[n-1,3],3)) )
            t.toc('Iteration time:')

    print('________________________________________________________________________________________________________________________________________')

    # post processing variables
    state = calc_postvars(state, gas, n)

    return state


# Roe averaging function on quantity phi
def roe_average( rhoL, rhoR, phiL, phiR ):

    avg = ( phiL*np.sqrt(rhoL) + phiR*np.sqrt(rhoR) ) / ( np.sqrt(rhoL) + np.sqrt(rhoR) )

    return avg

# see appendix A.7 in Computational Fluid Dynamics: Principles and Applications (Blazek, 2004)
def flux_jacobian( u, v, nx, ny, gam, E ):

    phi = (1/2) * (gam-1) * (u**2 + v**2)
    a1 = gam*E - phi
    a2 = gam - 1
    a3 = gam - 2
    V =  nx*u + ny*v

    A = np.array( ( 0,              nx,                 ny,                 0,              \
                    nx*phi-u*V,     V-a3*nx*u,          ny*u-a2*nx*v,       a2*nx,          \
                    ny*phi*v*V,     nx*v-a2*ny*u,       V-a3*ny*v,          a2*ny,          \
                    V*(phi-a1),     nx*a1-a2*u*V,       ny*a1-a2*v*V,       gam*V           ) ) . \
                    reshape( ( 4, 4 ) )

    return A

# equation 20 in Ren et Al. 2011
def f8( phi ):
    
    func = (np.minimum( phi * np.sqrt(4+(1-phi**2)**2)/(1+phi**2), 1 ))**8

    return func

# compute face quantity for determining F_c,1/2
def F_bar( Phi, U, normal, p, M, N ):

    Fbar = np.zeros( (M, N, 4) )
    for i in range(0, 4):
        Fbar[:,:,i] = U * Phi[:,:,i] + normal[:,:,i] * p

    return Fbar

# shock detector function
def shock_detect( p ):

    P = lambda p, pp: np.minimum(p/pp, pp/p)

    b = np.minimum( np.minimum( np.minimum( P(p[1:-1,1:-1],p[2:,1:-1]), P(p[1:-1,1:-1],p[2:,0:-2]) ), \
                    np.minimum( P(p[1:-1,1:-1],p[2:,2:]), P(p[1:-1,1:-1],p[1:-1,0:-2]) ) ), P(p[1:-1,1:-1],p[1:-1,2:]) )

    return b


def calc_postvars(state, gas):

    state.Mach = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 ) / \
                           thermo.calc_c( state.p, state.Q[:,:,0], gas.gamma_fn(gas.Cp, gas.Cv) )
    state.vel = np.sqrt( (state.Q[:,:,1]/state.Q[:,:,0])**2 + (state.Q[:,:,2]/state.Q[:,:,0])**2 )
    state.p0 = (1+((gas.gamma_fn(gas.Cp, gas.Cv)-1)/2)*state.Mach**2)** \
                   (gas.gamma_fn(gas.Cp, gas.Cv)/(gas.gamma_fn(gas.Cp, gas.Cv)-1)) * state.p
    state.T0 = (1+((gas.gamma_fn(gas.Cp, gas.Cv)-1)/2)*state.Mach**2) * state.T

    return state