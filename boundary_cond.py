

# set boundary conditions (inviscid walls at bottom, outlet at top and right)

def enforce_bc(domain, mesh, state, gas):

    state.p[:, 0] = state.p[:, 1]
    state.T[:, 0] = 300
    state.Q[:, 0:2, :] = invisc_wall(state.Q[:, 0:2, :], state.p[:, 0], state.T[:, 0], mesh.s_proj[:, 0:2, :], domain.M+2, gas)

    state.Q[domain.M+1, :, :] = state.Qn[domain.M+1, :, :]
    state.Q[:, domain.N+1, :] = state.Qn[:, domain.N+1, :]

    return state

# input Qwall[M+2, 2, 4]

def invisc_wall(Qwall, pwall, Twall, s_proj, M, gas):

    from helper import thermo
    import boundary

    u0 = Qwall[:, 0, 1] / Qwall[:, 0, 0]
    v0 = Qwall[:, 0, 2] / Qwall[:, 0, 0]

    u1 = Qwall[:, 1, 1] / Qwall[:, 1, 0]
    v1 = Qwall[:, 1, 2] / Qwall[:, 1, 0]

    boundary.slip(u0, v0, u1, v1, s_proj, M)

    Qwall[:, 0, 0] = pwall / (gas.R * Twall)
    Qwall[:, 0, 1] = u0 * Qwall[:, 0, 0]
    Qwall[:, 0, 2] = v0 * Qwall[:, 0, 0]
    Qwall[:, 0, 3] = thermo.calc_rho_et( pwall, Qwall[:, 0, 0], Qwall[:, 0, 1]/Qwall[:, 0, 0], Qwall[:, 0, 2]/Qwall[:, 0, 0], gas.gamma )
    
    return Qwall