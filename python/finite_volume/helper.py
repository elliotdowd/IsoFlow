# thermodynamics-related functions

class thermo: 
    def calc_rho_et( p, rho, u, v, gam ): 
        rho_et = (p/(gam-1)) + (1/2)*rho*(u**2 + v**2)
        return rho_et

    def calc_p( rho, Q4, u, v, gam ): 
        p = (gam-1) * ( Q4 - (1/2)*rho*(u**2 + v**2))
        return p

    def calc_c( p, rho, gam ):
        c = ( (gam * p) / rho )**0.5
        return c
    
    def calc_c_star( ht, gam ):
        cs = 2 * ht * ( (gam-1) / (gam+1) )
        return cs

class split:
    # van Leer Mach number splitting polynomials
    def Mvm( M ):
        import numpy as np
        Msplit = np.double(abs(M)<=1) * -0.25*(M - 1)**2 + \
                 np.double(abs(M)>1) * 0.5*(M-abs(M))
        return Msplit
    def Mvp( M ):
        import numpy as np
        Msplit = np.double(abs(M)<=1) * 0.25*(M + 1)**2 + \
                 np.double(abs(M)>1) * 0.5*(M+abs(M))
        return Msplit

    # first degree Mach number splitting polynomials
    def M1m( M ): 
        Msplit = 0.5 * ( M-abs(M) )
        return Msplit

    def M1p( M ):
        Msplit = 0.5 * ( M+abs(M) )
        return Msplit

    # second degree Mach number splitting polynomials
    def M2m( M ):
        import numpy as np
        M2split = np.double(np.abs(M)>=1)*split.M1m(M) + \
                  np.double(np.abs(M)<1)*-0.25*(M-1)**2
        return M2split
    def M2p( M ):
        import numpy as np
        M2split = np.double(np.abs(M)>=1)*split.M1p(M) + \
                  np.double(np.abs(M)<1)*0.25*(M+1)**2
        return M2split

    # fourth order splitting
    def M4m( M ): 
        import numpy as np
        M4minus =  np.double(np.abs(M)>=1) * split.M1m(M) + \
                   np.double(np.abs(M)<1) * (-(1/4) * (M-1)**2 - (1/8)*(M**2-1)**2); 
        return M4minus

    def M4p( M ):
        import numpy as np
        M4minus =  np.double(np.abs(M)>=1) * split.M1p(M) + \
                   np.double(np.abs(M)<1) * ((1/4) * (M+1)**2 + (1/8)*(M**2-1)**2); 
        return M4minus

    # pressure diffusion term
    def Mp( Mbar, M0, pL, pR, rho_half, c_half ):
        import numpy as np
        M = -(0.25/split.fa(M0))*np.maximum(1-Mbar**2, 0)*(pR-pL)/(rho_half*c_half**2)
        return M

    # first degree pressure splitting polynomials
    def P1m( M ):
        import numpy as np
        Psplit = np.double(np.abs(M)<=1) * 0.5*(1-M) + np.double(np.abs(M)>1) * 0.5*(M-np.abs(M))/M
        return Psplit

    def P1p( M ):
        import numpy as np
        Psplit = np.double(np.abs(M)<=1) * 0.5*(1+M) + np.double(np.abs(M)>1) * 0.5*(M+np.abs(M))/M
        return Psplit

    # third degree pressure splitting polynomials
    def P3m( M ):
        P3split = -(2-M)*split.M2m(M)
        return P3split

    def P3p( M ):
        P3split = -(2-M)*split.M2p(M)
        return P3split

    # fifth degree pressure splitting polynomials
    def P5m( M, a ):
        import numpy as np
        P5split = np.double(abs(M)>=1)*split.M1m(M)/M + \
                  np.double(abs(M)<1)*(split.M2m(M)*((-2-M)-(16*a*M*split.M2m(M))))
        return P5split
    def P5p( M, a ):
        import numpy as np
        P5split = np.double(abs(M)>=1)*split.M1p(M)/M + \
                  np.double(abs(M)<1)*(split.M2p(M)*((2-M)+(16*a*M*split.M2p(M))))
        return P5split

    # velocity diffusion term
    def Pu(ML, MR, UL, UR, rhoL, rhoR, c_half, a):
        P = -0.75*split.P5p(ML, a )*split.P5m(MR, a )*(rhoL+rhoR)*c_half*(UR-UL)
        return P

    # reference Mach number function
    def M0( Mbar, Minf ):
        import numpy as np
        M = np.minimum(1, np.maximum(Mbar**2, Minf**2))
        return M

    # Mach shaping function
    def fa( M0 ):
        f = M0*(2-M0)
        return f

    # alpha coefficient for fifth degree polynomials
    def alpha( fa ):
        a = (3/16) * (-4 + 5*fa**2)
        return a



# gradient function

def grad(x, y, f):
    import numpy as np
    # center divided difference approximation for gradient magnitude
    # note to leave one halo cell on either side of mesh in each dimension

    # central difference approximation

    fx = (f[2:,1:-1] - f[0:-2,1:-1]) / (2*(x[2:,1:-1]-x[1:-1,1:-1]))
    fy = (f[1:-1,2:] - f[1:-1,0:-2]) / (2*(y[1:-1,2:]-y[1:-1,1:-1]))
    
    df = np.sqrt(fx**2 + fy**2)

    return df
