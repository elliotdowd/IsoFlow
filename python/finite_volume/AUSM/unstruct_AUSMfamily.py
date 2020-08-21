import numpy as np
from python.finite_volume.helper import split, thermo
from python.boundary.unstruct_boundary_cond import unstruct_covariant, unstruct_boundary_cond
from python.finite_volume.timestepping import local_unstruct_timestep
from pytictoc import TicToc

def unstruct_AUSM( domain, mesh, state, parameters, gas ):

    print('AUSM Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # number of elements and faces
    mesh.M = len( mesh.elements )
    mesh.N = len( mesh.faces )

    # initialize speed of sound arrays
    c_L = np.zeros( mesh.N, dtype='float', order='F' )
    c_R = np.zeros( mesh.N, dtype='float', order='F' )
    c_half = np.zeros( mesh.N, dtype='float', order='F' )

    # left and right side Mach number arrays
    M_L = np.zeros( mesh.N, dtype='float', order='F' )
    M_R = np.zeros( mesh.N, dtype='float', order='F' )
    M_half = np.zeros( mesh.N, dtype='float', order='F' )
    mdot_half = np.zeros( mesh.N, dtype='float', order='F' )
    p_half = np.zeros( mesh.N, dtype='float', order='F' )

    # initialize phi and p state vectors
    Phi = np.zeros( (mesh.M, 4) )
    Phi[:,0] = np.ones( mesh.M )

    Phibound = np.zeros( (len(mesh.bdry_ind[0]), 4) )
    Phibound[:,0] = np.ones( len(mesh.bdry_ind[0]) )

    P_half = np.zeros( (mesh.N, 4) )
    P_half[:,0] = np.zeros( mesh.N )
    P_half[:,3] = np.zeros( mesh.N )

    F_hat = np.zeros( (mesh.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( (mesh.M, 4), dtype='float', order='F' )
    if state.n == 0:
        state.res = np.ones( [parameters.iterations + 1, 4] )
    else:
        state.res = np.vstack( (state.res, np.ones( [parameters.iterations + 1, 4] )) )

    mesh.dV4 = np.transpose( np.vstack( [mesh.dV, mesh.dV, mesh.dV, mesh.dV] ) )

    t = TicToc()

    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        t.tic()

        n = n+1
        state.n = state.n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # local timestepping
        state = local_unstruct_timestep( mesh, state, parameters, gas )
        dt = np.transpose( np.vstack( [state.dt, state.dt, state.dt, state.dt] ) )

        # simplify variable notation from state vector
        state.u = state.Q[:,1] / state.Q[:,0]
        state.v = state.Q[:,2] / state.Q[:,0]
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,0], state.u, state.v, gas.gamma ) \
                                    / state.Q[:,0] + state.p / state.Q[:,0]


        # initialize Phi vector components
        Phi[:,1] = state.u
        Phi[:,2] = state.v
        Phi[:,3] = state.ht
        Phibound[:,1] = state.Qbound[:,1]/state.Qbound[:,0]
        Phibound[:,2] = state.Qbound[:,2]/state.Qbound[:,0]
        Phibound[:,3] = thermo.calc_rho_et(state.pbound, state.Qbound[:,0], state.ubound, state.vbound, gas.gamma ) / state.Qbound[:,0] + \
                                           state.pbound / state.Qbound[:,0]

        # characteristic sound speed at each element centroid
        c_st = thermo.calc_c_star( state.ht, gas.gamma )
        c_stbound = thermo.calc_c_star( state.Qbound[:,3]/state.Qbound[:,0] + state.pbound/state.Qbound[:,0], gas.gamma )

        cr = 1e-60

        # loop through cell faces for finite volume integration
        for i, face in enumerate( mesh.faces ):

            # pointers from faces to left (interior) and right (exterior) elements
            L = int( mesh.face_pairs[i,0] )
            R = int( mesh.face_pairs[i,1] )

            # sound speeds at each cell face (account for boundary cells)
            if R < 0:
                bound_i = int(np.where(np.isin(mesh.bdry_ind,i))[1])
                c_L[i] = c_st[L] / np.maximum( np.sqrt(c_st[L]), state.U[i,0] )
                c_R[i] = c_stbound[bound_i] / np.maximum( np.sqrt(c_stbound[bound_i]), -state.Ubound[bound_i,1] ) 

                c_half[i] = np.min( [ c_L[i], c_R[i] ] )

                M_L[i] = state.U[i,0] / c_half[i]
                M_R[i] = state.Ubound[bound_i,1] / c_half[i]

                # interface Mach numbers
                M_half[i] = split.Mvp( M_L[i] ) + split.Mvm( M_R[i] )

                # interface mass flux
                mdot_half[i] = c_half[i] * M_half[i] * ( np.double(M_half[i]>0) * state.Q[L,0] + np.double(M_half[i]<=0) * state.Qbound[bound_i,0] )

                # calculate pressure flux at cell interfaces
                p_half[i] = split.P1p( M_L[i]+cr )*state.p[L] + split.P1m( M_R[i]+cr )*state.pbound[bound_i]

                # pressure flux vector
                P_half[i,1] = p_half[i] * mesh.n[i,0]
                P_half[i,2] = p_half[i] * mesh.n[i,1]
            
                F_hat[i,:] = (1/2) * mdot_half[i] * ( Phi[L,:] + Phibound[bound_i,:] ) - \
                             (1/2) * abs(mdot_half[i]) * ( Phibound[bound_i,:] - Phi[L,:] ) + \
                                   + P_half[i,:]

                state.residual[L] = state.residual[L] - ( F_hat[i,:] * mesh.dS[i] )

            else:
                c_L[i] = c_st[L] / np.maximum( np.sqrt(c_st[L]), state.U[i,0] )
                c_R[i] = c_st[R] / np.maximum( np.sqrt(c_st[R]), -state.U[i,1] )

                c_half[i] = np.min( [ c_L[i], c_R[i] ] )

                M_L[i] = state.U[i,0] / c_half[i]
                M_R[i] = state.U[i,1] / c_half[i]

                # interface Mach numbers
                M_half[i] = split.Mvp( M_L[i] ) + split.Mvm( M_R[i] )

                # interface mass flux
                mdot_half[i] = c_half[i] * M_half[i] * ( np.double(M_half[i]>0) * state.Q[L,0] + np.double(M_half[i]<=0) * state.Q[R,0] )

                # calculate pressure flux at cell interfaces
                p_half[i] = split.P1p( M_L[i]+cr )*state.p[L] + split.P1m( M_R[i]+cr )*state.p[R]
            
                # pressure flux vector
                P_half[i,1] = p_half[i] * mesh.n[i,0]
                P_half[i,2] = p_half[i] * mesh.n[i,1]
            
                F_hat[i,:] = (1/2) * mdot_half[i] * ( Phi[L,:] + Phi[R,:] ) - \
                             (1/2) * abs(mdot_half[i]) * ( Phi[R,:] - Phi[L,:] ) + \
                                   + P_half[i,:]

                state.residual[L] = state.residual[L] - ( F_hat[i,:] * mesh.dS[i] )
                state.residual[R] = state.residual[R] + ( F_hat[i,:] * mesh.dS[i] )

        # update state vector
        # state.Q = state.Qn + state.residual / mesh.dV4
        state.Q = state.Qn + dt * state.residual

        # L_inf-norm residual
        state.res[state.n-1,0] = np.log10( np.max(state.residual[:,0] * mesh.dV) ) 
        state.res[state.n-1,1] = np.log10( np.max(state.residual[:,1] * mesh.dV) ) 
        state.res[state.n-1,2] = np.log10( np.max(state.residual[:,2] * mesh.dV) )
        state.res[state.n-1,3] = np.log10( np.max(state.residual[:,3] * mesh.dV) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,0], state.Q[:,3], state.u, state.v, gas.gamma )
        state.T = state.p / (gas.R * state.Q[:,0])

        # update covariant velocities
        state = unstruct_covariant( mesh, state )

        # enforce boundary conditions
        state = unstruct_boundary_cond( domain, mesh, state, parameters, gas )

        # print iteration output
        print('Iteration: ' + str(state.n) + '    ' + str(round(state.res[n-1,0],3)) + '    ' + str(round(state.res[n-1,1],3)) + \
                                                '    ' + str(round(state.res[n-1,2],3)) + '    ' + str(round(state.res[n-1,3],3)) )
        t.toc('Iteration time:')


    print('________________________________________________________________________________________________________________________________________')

    state = calc_postvars(state, gas)

    return state


def avg_vars( domain, mesh, state, parameters, gas ):
    print('AUSM Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # number of elements and faces
    mesh.M = len( mesh.elements )
    mesh.N = len( mesh.faces )

    # initialize speed of sound arrays
    c_L = np.zeros( mesh.N, dtype='float', order='F' )
    c_R = np.zeros( mesh.N, dtype='float', order='F' )
    c_half = np.zeros( mesh.N, dtype='float', order='F' )

    # left and right side Mach number arrays
    M_L = np.zeros( mesh.N, dtype='float', order='F' )
    M_R = np.zeros( mesh.N, dtype='float', order='F' )
    M_half = np.zeros( mesh.N, dtype='float', order='F' )
    mdot_half = np.zeros( mesh.N, dtype='float', order='F' )
    p_half = np.zeros( mesh.N, dtype='float', order='F' )

    # initialize phi and p state vectors
    Phi = np.zeros( (mesh.M, 4) )
    Phi[:,0] = np.ones( mesh.M )

    Phi_half = np.zeros( (mesh.N, 4) )

    Phibound = np.zeros( (len(mesh.bdry_ind[0]), 4) )
    Phibound[:,0] = np.ones( len(mesh.bdry_ind[0]) )

    P_half = np.zeros( (mesh.N, 4) )
    P_half[:,0] = np.zeros( mesh.N )
    P_half[:,3] = np.zeros( mesh.N )

    # initialize averaged quantities
    rho_half = np.zeros( mesh.N )
    U_half = np.zeros( mesh.N )
    p_half = np.zeros( mesh.N )
    u_half = np.zeros( mesh.N )
    v_half = np.zeros( mesh.N )
    ht_half = np.zeros( mesh.N )

    F_hat = np.zeros( (mesh.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( (mesh.M, 4), dtype='float', order='F' )
    if state.n == 0:
        state.res = np.ones( [parameters.iterations + 1, 4] )
    else:
        state.res = np.vstack( (state.res, np.ones( [parameters.iterations + 1, 4] )) )

    mesh.dV4 = np.transpose( np.vstack( [mesh.dV, mesh.dV, mesh.dV, mesh.dV] ) )

    t = TicToc()

    while n <= parameters.iterations and max(state.res[n-1+np.int(n<10),:]) > parameters.tolerance: 

        t.tic()

        n = n+1
        state.n = state.n+1

        # state at previous timestep, use for outflow BCs
        state.Qn = state.Q

        # local timestepping
        state = local_unstruct_timestep( mesh, state, parameters, gas )
        dt = np.transpose( np.vstack( [state.dt, state.dt, state.dt, state.dt] ) )

        # simplify variable notation from state vector
        state.u = state.Q[:,1] / state.Q[:,0]
        state.v = state.Q[:,2] / state.Q[:,0]
        state.ht = thermo.calc_rho_et(state.p, state.Q[:,0], state.u, state.v, gas.gamma ) \
                                    / state.Q[:,0] + state.p / state.Q[:,0]

        ht_bound = thermo.calc_rho_et(state.pbound, state.Qbound[:,0], state.ubound, state.vbound, gas.gamma) / \
                                                    state.Qbound[:,0] + \
                                                    state.pbound / state.Qbound[:,0]

        # loop through cell faces for finite volume integration
        for i, face in enumerate( mesh.faces ):

            # pointers from faces to left (interior) and right (exterior) elements
            L = int( mesh.face_pairs[i,0] )
            R = int( mesh.face_pairs[i,1] )

            # sound speeds at each cell face (account for boundary cells)
            if R < 0:
                bound_i = int(np.where(np.isin(mesh.bdry_ind,i))[1])

                rho_half[i] = (1/2) * ( state.Q[L,0] + state.Qbound[bound_i,0] )
                U_half[i] = (1/2) * ( state.U[i,0] + state.Ubound[bound_i,1] )
                p_half[i] = (1/2) * ( state.p[L] + state.pbound[bound_i] )

                u_half[i] = (1/2) * ( state.u[L] + state.ubound[bound_i] )
                v_half[i] = (1/2) * ( state.v[L] + state.vbound[bound_i] )

                ht_half[i] = (1/2) * ( state.ht[L] + ht_bound[bound_i] )

            else:

                rho_half[i] = (1/2) * ( state.Q[L,0] + state.Q[R,0] )
                U_half[i] = (1/2) * ( state.U[i,0] + state.U[i,1] )
                p_half[i] = (1/2) * ( state.p[L] + state.p[R] )

                u_half[i] = (1/2) * ( state.u[L] + state.u[R] )
                v_half[i] = (1/2) * ( state.v[L] + state.v[R] )

                ht_half[i] = (1/2) * ( state.ht[L] + state.ht[R] )

            Phi_half[i,:] = np.array( [1, u_half[i], v_half[i], ht_half[i]] )
            P_half[i,:] = np.array( [0, mesh.n[i,0]*p_half[i], mesh.n[i,1]*p_half[i], 0] )

            F_hat[i,:] = np.array( [rho_half[i], rho_half[i], rho_half[i], rho_half[i]] ) * \
                         np.array( [U_half[i], U_half[i], U_half[i], U_half[i]] ) * \
                         Phi_half[i,:] + P_half[i,:]

            state.residual[L] = state.residual[L] + ( F_hat[i,:] * mesh.dS[i] )
            if R >= 0:
                state.residual[R] = state.residual[R] - ( F_hat[i,:] * mesh.dS[i] )


        # update state vector
        # state.Q = state.Qn + state.residual / mesh.dV4
        state.Q = state.Qn + dt * state.residual

        # L_inf-norm residual
        state.res[state.n-1,0] = np.log10( np.max(state.residual[:,0] * mesh.dV) ) 
        state.res[state.n-1,1] = np.log10( np.max(state.residual[:,1] * mesh.dV) ) 
        state.res[state.n-1,2] = np.log10( np.max(state.residual[:,2] * mesh.dV) )
        state.res[state.n-1,3] = np.log10( np.max(state.residual[:,3] * mesh.dV) ) 

        # update cell temperatures and pressures
        state.p = thermo.calc_p( state.Q[:,0], state.Q[:,3], state.u, state.v, gas.gamma )
        state.T = state.p / (gas.R * state.Q[:,0])

        # update covariant velocities
        state = unstruct_covariant( mesh, state )

        # enforce boundary conditions
        state = unstruct_boundary_cond( domain, mesh, state, parameters, gas )

        # print iteration output
        print('Iteration: ' + str(state.n) + '    ' + str(round(state.res[n-1,0],3)) + '    ' + str(round(state.res[n-1,1],3)) + \
                                                '    ' + str(round(state.res[n-1,2],3)) + '    ' + str(round(state.res[n-1,3],3)) )
        t.toc('Iteration time:')


    print('________________________________________________________________________________________________________________________________________')

    state = calc_postvars(state, gas)

    return state


# state variable post processing
def calc_postvars(state, gas):

    state.Mach = np.sqrt( (state.Q[:,1]/state.Q[:,0])**2 + (state.Q[:,2]/state.Q[:,0])**2 ) / \
                           thermo.calc_c( state.p, state.Q[:,0], gas.gamma )
    state.vel = np.sqrt( (state.Q[:,1]/state.Q[:,0])**2 + (state.Q[:,2]/state.Q[:,0])**2 )
    state.isenRatio = (1+((gas.gamma-1)/2)*state.Mach**2)
    state.p0Ratio = state.isenRatio**(gas.gamma/(gas.gamma-1))
    state.p0 = state.p0Ratio * state.p
    state.T0 = state.isenRatio * state.T

    return state