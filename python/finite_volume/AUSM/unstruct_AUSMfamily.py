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
    Phi_half = np.zeros( (mesh.N, 4) )

    Phibound = np.zeros( (len(mesh.bdry_ind[0]), 4) )
    Phibound[:,0] = np.ones( len(mesh.bdry_ind[0]) )

    P_half = np.zeros( (mesh.N, 4) )
    P_half[:,0] = np.zeros( mesh.N )
    P_half[:,3] = np.zeros( mesh.N )

    Fc_half = np.zeros( (mesh.N, 4), dtype='float', order='F' )
    F_half = np.zeros( (mesh.N, 4), dtype='float', order='F' )

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
        state.ht = ( thermo.calc_rho_et(state.p, state.Q[:,0], state.u, state.v, gas.gamma ) \
                                    / state.Q[:,0] + state.p ) / state.Q[:,0]


        # initialize Phi vector components
        Phi[:,0] = state.Q[:,0]
        Phi[:,1] = state.u * state.Q[:,0]
        Phi[:,2] = state.v * state.Q[:,0]
        Phi[:,3] = state.ht
        Phibound[:,0] = state.Qbound[:,0]
        Phibound[:,1] = state.Qbound[:,1]
        Phibound[:,2] = state.Qbound[:,2]
        Phibound[:,3] = (thermo.calc_rho_et(state.pbound, state.Qbound[:,0], state.ubound, state.vbound, gas.gamma ) / state.Qbound[:,0] + state.pbound) / state.Qbound[:,0]

        # characteristic sound speed at each element centroid
        # c_st = thermo.calc_c_star( state.ht, gas.gamma )
        # c_stbound = thermo.calc_c_star( state.Qbound[:,3]/state.Qbound[:,0] + state.pbound/state.Qbound[:,0], gas.gamma )

        cr = 1e-60

        # loop through control volumes
        for L, faces in enumerate( mesh.elem_to_face ):

            # loop through faces belonging to each CV
            for j, face in enumerate( faces ):

                # pointers from faces to left (interior) and right (exterior) elements
                # L = int( mesh.face_pairs[face,0] )
                R = int( mesh.face_pairs[face,1] )

                if R < 0:
                    bound_i = int(np.where(np.isin(mesh.bdry_ind,face))[1])
                else:
                    bound_i = 0

                # interface sound speed
                c_L[face] = state.c[L]
                c_R[face] = state.c[R]*np.double(R>-1) + thermo.calc_c( state.pbound[bound_i], state.Qbound[bound_i,0], gas.gamma )*np.double(R<0)
                c_half[face] = (1/2) * ( c_L[face] + c_R[face] )

                # left and right convective velocities
                Vl = state.U[face,0]
                Vr = state.U[face,1]*np.double(R>-1) + state.Ubound[bound_i,1]*np.double(R<0)

                if mesh.ccw[L][j] == 1:
                    N = mesh.n[face]
                else:
                    Vl = -Vl
                    Vr = -Vr
                    N = -mesh.n[face]


                M_L[face] = Vl / c_half[face]
                M_R[face] = Vr / c_half[face]

                # interface Mach numbers
                M_half[face] = split.Mvp( M_L[face] ) + split.Mvm( M_R[face] )

                # interface mass flux
                mdot_half[face] = c_half[face] * M_half[face] * ( np.double(M_half[face]>0)*state.Q[L,0] + \
                                                                  np.double(M_half[face]<=0)*state.Q[R,0] * np.double(R>-1) + \
                                                                  np.double(M_half[face]<=0)*state.Qbound[bound_i,0] * np.double(R<0) )

                # calculate pressure flux at cell interfaces
                p_half[face] = split.P1p( M_L[face]+cr )*state.p[L] + split.P1m( M_R[face]+cr )*state.p[R]*np.double(R>-1) + \
                                                                      split.P1m( M_R[face]+cr )*state.pbound[bound_i]*np.double(R<0)
            
                # pressure flux vector
                P_half[face,1] = p_half[face] * N[0]
                P_half[face,2] = p_half[face] * N[1]

                PhiL = Phi[L,:]
                PhiR = Phi[R,:]*np.double(R>-1) + Phibound[bound_i,:]*np.double(R<0)

                Fc_half[face,:] = (1/2) * mdot_half[face] * (PhiL + PhiR) - \
                                  (1/2) * np.abs(mdot_half[face]) * (PhiR - PhiL)

                # F_half[face,:] = ( Fc_half[face,:] + P_half[face,:] ) * mesh.dS[face]
                F_half[face,:] = ( Fc_half[face,:] ) * mesh.dS[face]

            state.residual[L] = np.sum( F_half[faces,:], axis=0 )


        # update state vector
        state.Q = state.Qn - dt * state.residual / mesh.dV4

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


def centered_scheme( domain, mesh, state, parameters, gas ):

    print('AUSM Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # number of elements and faces
    mesh.M = len( mesh.elements )
    mesh.N = len( mesh.faces )

    # cell centered fluxes
    C = np.zeros( [ mesh.M, 4 ], dtype='float' )

    # E = np.zeros( (mesh.N, 4), dtype='float', order='F' )
    # F = np.zeros( (mesh.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( (mesh.M, 4), dtype='float' )
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
        state.ht = ( thermo.calc_rho_et(state.p, state.Q[:,0], state.u, state.v, gas.gamma ) \
                                    / state.Q[:,0] + state.p ) / state.Q[:,0]

        cr = 1e-60

        # loop through control volumes
        for L, faces in enumerate( mesh.elem_to_face ):

            # loop through faces belonging to each CV
            for j, face in enumerate( faces ):

                # pointers from faces to left (interior) and right (exterior) elements
                R = int( mesh.face_pairs[face,1] )

                QL = state.Q[L,:]
                pL = state.p[L]

                if L == 258:
                    L

                if R < 0:
                    bound_i = int(np.where(np.isin(mesh.bdry_ind,face))[1])
                    QR = state.Qbound[bound_i,:]
                    pR = state.pbound[bound_i]
                else:
                    bound_i = 0
                    QR = state.Q[R,:]
                    pR = state.p[R]

                if mesh.ccw[L][j] == 1:
                    N = mesh.n[face]
                    sign = 1
                else:
                    N = -mesh.n[face]
                    sign = -1

                Q_half = (1/2) * ( QL + QR )
                p_half = (1/2) * ( pL + pR )

                E = calc_E( Q_half, p_half )
                F = calc_F( Q_half, p_half )

                C[L] = C[L] + E*sign*mesh.dy[face] - F*sign*mesh.dx[face]

            state.residual[L] = C[L]


        # update state vector
        state.Q = state.Qn - dt * state.residual / mesh.dV4

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


def centered_face( domain, mesh, state, parameters, gas ):

    print('AUSM Scheme: ' + 'CFL = ' + str(parameters.CFL))
    print('________________________________________________________________________________________________________________________________________')
    print('                          mass          u            v        energy')
    print('________________________________________________________________________________________________________________________________________')

    n = 0

    # number of elements and faces
    mesh.M = len( mesh.elements )
    mesh.N = len( mesh.faces )

    # cell centered fluxes
    C = np.zeros( [ mesh.M, 4 ], dtype='float' )

    # E = np.zeros( (mesh.N, 4), dtype='float', order='F' )
    # F = np.zeros( (mesh.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( (mesh.M, 4), dtype='float' )
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
        state.ht = ( thermo.calc_rho_et(state.p, state.Q[:,0], state.u, state.v, gas.gamma ) \
                                    / state.Q[:,0] + state.p ) / state.Q[:,0]

        cr = 1e-60

        # loop through faces
        for i, elems in enumerate( mesh.face_pairs ):

            # pointers from faces to left (interior) and right (exterior) elements
            L = elems[0]
            R = elems[1]

            QL = state.Q[L,:]
            pL = state.p[L]

            if R < 0:
                bound_i = int(np.where(np.isin(mesh.bdry_ind,i))[1])
                QR = state.Qbound[bound_i,:]
                pR = state.pbound[bound_i]
            else:
                bound_i = 0
                QR = state.Q[R,:]
                pR = state.p[R]

            Q_half = (1/2) * ( QL + QR )
            p_half = (1/2) * ( pL + pR )

            E = calc_E( Q_half, p_half )
            F = calc_F( Q_half, p_half )

            C[L] = C[L] + E*mesh.dy[i] - F*mesh.dx[i]

            state.residual[L] = C[L]


        # update state vector
        state.Q = state.Qn - dt * state.residual / mesh.dV4

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


def unstruct_AUSM_facebased( domain, mesh, state, parameters, gas ):

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

    Fc_half = np.zeros( (mesh.N, 4), dtype='float', order='F' )
    F_half = np.zeros( (mesh.N, 4), dtype='float', order='F' )

    state.residual = np.zeros( (mesh.M, 4), dtype='float', order='F' )
    if state.n == 0:
        state.res = np.ones( [parameters.iterations + 1, 4] )
    else:
        state.res = np.vstack( (state.res, np.ones( [parameters.iterations + 1, 4] )) )

    mesh.dV4 = np.transpose( np.vstack( [mesh.dV, mesh.dV, mesh.dV, mesh.dV] ) )

    cr = 1e-30

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
        state.ht = ( thermo.calc_rho_et(state.p, state.Q[:,0], state.u, state.v, gas.gamma ) \
                                    / state.Q[:,0] + state.p ) / state.Q[:,0]

        # initialize Phi vector components
        Phi[:,0] = np.ones( mesh.M )
        Phi[:,1] = state.u
        Phi[:,2] = state.v
        Phi[:,3] = state.ht

        Phibound[:,0] = np.ones( len(state.Qbound[:,0]) )
        Phibound[:,1] = state.Qbound[:,1] / state.Qbound[:,0]
        Phibound[:,2] = state.Qbound[:,2] / state.Qbound[:,0]
        Phibound[:,3] = (thermo.calc_rho_et(state.pbound, state.Qbound[:,0], state.ubound, state.vbound, gas.gamma ) / state.Qbound[:,0] + state.pbound) / state.Qbound[:,0]


        # loop through cell faces for finite volume integration
        for i, face in enumerate( mesh.faces ):

            # pointers from faces to left (interior) and right (exterior) elements
            L = int( mesh.face_pairs[i,0] )
            R = int( mesh.face_pairs[i,1] )

            # sound speeds at each cell face (account for boundary cells)
            if R < 0:
                bound_i = int(np.where(np.isin(mesh.bdry_ind, i))[1])
            else:
                bound_i = 0

            # c_L[i] = c_st[L] / np.maximum( np.sqrt(c_st[L]), state.U[i,0] )
            # c_R[i] = c_st[R] / np.maximum( np.sqrt(c_st[R]), -state.U[i,1] )
            c_L[i] = state.c[L]
            c_R[i] = state.c[R]

            c_half[i] = (1/2) * ( c_L[i] + c_R[i] )

            # left and right convective velocities
            Vl = state.U[i,0]
            Vr = state.U[i,1]*np.double(R>-1) + state.Ubound[bound_i,1]*np.double(R<0)

            M_L[i] = Vl / c_half[i]
            M_R[i] = Vr / c_half[i]

            # interface Mach numbers
            M_half[i] = split.Mvp( M_L[i] ) + split.Mvm( M_R[i] )

            # interface mass flux
            mdot_half[i] = c_half[i] * M_half[i] * ( np.double(M_half[i]>0)*state.Q[L,0] + \
                                                     np.double(M_half[i]<=0)*state.Q[R,0] * np.double(R>-1) + \
                                                     np.double(M_half[i]<=0)*state.Qbound[bound_i,0] * np.double(R<0) )

            # calculate pressure flux at cell interfaces
            p_half[i] = split.P1p( M_L[i]+cr )*state.p[L] + split.P1m( M_R[i]+cr )*state.p[R]*np.double(R>-1) + \
                                                            split.P1m( M_R[i]+cr )*state.pbound[bound_i]*np.double(R<0)
        
            # pressure flux vector
            P_half[i,1] = p_half[i] * mesh.n[i,0]
            P_half[i,2] = p_half[i] * mesh.n[i,1]

            PhiL = Phi[L,:]
            PhiR = Phi[R,:]*np.double(R>-1) + Phibound[bound_i,:]*np.double(R<0)

            Fc_half[i,:] = (1/2)*mdot_half[i]*(PhiL+PhiR) - \
                           (1/2)*np.abs(mdot_half[i])*(PhiR-PhiL)

            F_half[i,:] = Fc_half[i,:] + P_half[i,:]

            state.residual[L] = state.residual[L] + ( F_half[i,:] * mesh.dS[i] )
            if R < -1:
                state.residual[R] = state.residual[R] - ( F_half[i,:] * mesh.dS[i] )

        # update state vector
        state.Q = state.Qn - dt * state.residual / mesh.dV4

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

    Fc_hat = np.zeros( (mesh.N, 4), dtype='float', order='F' )
    F_half = np.zeros( (mesh.N, 4), dtype='float', order='F' )

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
        state.ht = (thermo.calc_rho_et(state.p, state.Q[:,0], state.u, state.v, gas.gamma )  + state.p) / state.Q[:,0]

        ht_bound = (thermo.calc_rho_et(state.pbound, state.Qbound[:,0], state.ubound, state.vbound, gas.gamma) + state.pbound) / state.Qbound[:,0]

        # loop through cell faces for finite volume integration
        for i, face in enumerate( mesh.faces ):

            # pointers from faces to left (interior) and right (exterior) elements
            L = int( mesh.face_pairs[i,0] )
            R = int( mesh.face_pairs[i,1] )

            # sound speeds at each cell face (account for boundary cells)
            if R < 0:
                bound_i = int(np.where(np.isin(mesh.bdry_ind,i))[1])

                rhoL = state.Q[L,0]
                rhoR = state.Qbound[bound_i,0]

                UL = state.U[i,0]
                UR = state.Ubound[bound_i,1]

                pL = state.p[L]
                pR = state.pbound[bound_i]

                uL = state.u[L]
                uR = state.ubound[bound_i]
                vL = state.v[L]
                vR = state.vbound[bound_i]

                htL = state.ht[L]
                htR = ht_bound[bound_i]

                QL = state.Q[L,:]
                QR = state.Qbound[bound_i,:]

            else:

                rhoL = state.Q[L,0]
                rhoR = state.Q[R,0]

                UL = state.U[i,0]
                UR = state.U[i,1]

                pL = state.p[L]
                pR = state.p[R]

                uL = state.u[L]
                uR = state.u[R]
                vL = state.v[L]
                vR = state.v[R]

                htL = state.ht[L]
                htR = state.ht[R]

                QL = state.Q[L,:]
                QR = state.Q[R,:]

            PhiL = np.array( [1, uL, vL, htL ])
            PhiR = np.array( [1, uR, vR, htR ])

            # mdot4 = np.array( [rho_half[i], rho_half[i], rho_half[i], rho_half[i]] ) * \
            #         np.array( [U_half[i], U_half[i], U_half[i], U_half[i]] )
            # Phi_half[i,:] = np.array( [1, u_half[i], v_half[i], ht_half[i]] )
            # P_half[i,:] = np.array( [0, mesh.n[i,0]*p_half[i], mesh.n[i,1]*p_half[i], 0] )
            ubar = u_half*mesh.dy[i]/mesh.dS[i] - v_half*mesh.dx[i]/mesh.dS[i]
            vbar = u_half*mesh.dx[i]/mesh.dS[i] + v_half*mesh.dy[i]/mesh.dS[i]

            # HdS = np.array( [ rho_half*ubar, rho_half*ubar**2, ubar*vbar, ubar*(ht_half) ] )

            # F_half[i,:] = np.matmul( np.linalg.inv( transform( mesh.dx[i], mesh.dy[i], mesh.dS[i] ) ), HdS )

            P_half[i,1] = (1/2) * (pL+pR) * mesh.n[i,0]
            P_half[i,2] = (1/2) * (pL+pR) * mesh.n[i,1]

            F_half[i,:] = (1/8) * (rhoL+rhoR) * (UL+UR) * (PhiL+PhiR) + P_half[i,:]

            state.residual[L] = state.residual[L] + dt[L] * ( F_half[i,:] * mesh.dS[i] )
            if R >= 0:
                state.residual[R] = state.residual[R] - dt[R] * ( F_half[i,:] * mesh.dS[i] )

            # F_hat[i,:] = mdot4 * Phi_half[i,:]# + P_half[i,:]

            # state.residual[L] = state.residual[L] + dt[L] * ( F_hat[i,:] * mesh.dS[i] )
            # if R >= 0:
            #     state.residual[R] = state.residual[R] - dt[R] * ( F_hat[i,:] * mesh.dS[i] )


        # update state vector
        # state.Q = state.Qn - state.residual / mesh.dV4
        state.Q = state.Q - state.residual

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


# rotation matrix for flux vectors
def transform( dx, dy, ds ):

    T = np.array( [ 1,      0,      0,      0, \
                    0,      dy/ds, -dx/ds,  0, \
                    0,      dx/ds,  dy/ds,  0, \
                    0,      0,      0,      1 ] ).reshape( [4, 4] )

    return T


def calc_E( Q, p ):

    E = np.zeros( 4 )

    E[0] = Q[1]
    E[1] = Q[1]**2/Q[0] + p
    E[2] = Q[1]*Q[2] / Q[0]
    E[3] = Q[3]*Q[1]

    return E

def calc_F( Q, p ):

    F = np.zeros( 4 )

    F[0] = Q[2]
    F[1] = Q[1]*Q[2] / Q[0]
    F[2] = Q[2]**2/Q[0] + p
    F[3] = Q[3]*Q[2]

    return F


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