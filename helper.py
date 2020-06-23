# thermodynamics-related functions

class thermo: 
    def calc_rho_et( p, rho, u, v, gam ): 
        rho_et = (p/(gam-1)) + 0.5*rho*(u**2 + v**2)
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

    # first degree pressure splitting polynomials
    def P1m( M ):
        import numpy as np
        Psplit = np.double(abs(M)<=1) * 0.5*(1-M) + np.double(abs(M)>1) * 0.5*(M-abs(M))/M
        return Psplit

    def P1p( M ):
        import numpy as np
        Psplit = np.double(abs(M)<=1) * 0.5*(1+M) + np.double(abs(M)>1) * 0.5*(M+abs(M))/M
        return Psplit

    # third degree pressure splitting polynomials
    def P3m( M ):
        P3split = -(2-M)*split.M2m(M)
        return P3split

    def P3p( M ):
        P3split = -(2-M)*split.M2p(M)
        return P3split



# matrix functions
class matrix:
    def m22mul( A, B ):
        import numpy as np
        result = np.zeros((2,2))
        for i in range(len(A)):
            # iterate through columns of A
            for j in range(len(B[0])):
                # iterate through rows of B
                for k in range(len(B)):
                    result[i][j] += A[i][k] * B[k][j]
        return result
