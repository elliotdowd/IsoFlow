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
    def Mm( M ): 
        import numpy as np
        Msplit = 0.5 * ( M-np.abs(M) )
        return Msplit

    def Mp( M ):
        import numpy as np
        Msplit = 0.5 * ( M+np.abs(M) )
        return Msplit

    def P1m( M ):
        import numpy as np
        Psplit = np.double(abs(M)<=1) * 0.5*(1-M) + np.double(abs(M)>1) * 0.5*(M-np.abs(M))/M
        return Psplit

    def P1p( M ):
        import numpy as np
        Psplit = np.double(abs(M)<=1) * 0.5*(1+M) + np.double(abs(M)>1) * 0.5*(M+np.abs(M))/M
        return Psplit
