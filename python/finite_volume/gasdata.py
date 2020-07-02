###########################################################################
## Calorically perfect gases (constant Cp, Cv)
###########################################################################	

# calorically perfect air
class air_cpg:
    gamma_p = 1
    Cp_p = 1006
    Cv_p = 718
    R_p = Cp_p - Cv_p
    theta = 3055.556

    def Cp_fn( k_p, Cp_p, th, T ):
        import numpy as np
        Cp = Cp_p * (1+(k_p-1)/k_p * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cp

    def Cv_fn( k_p, Cv_p, th, T ):
        import numpy as np
        Cv = Cv_p * (1+(k_p-1)/k_p * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cv

    def R_fn( Cp, Cv ):
        R = Cp - Cv
        return R

    def gamma_fn( Cp, Cv ):
        gamma = Cp / Cv
        return gamma


###########################################################################
## Thermally perfect gases (constant Cp, Cv)
###########################################################################	

# thermally perfect air
class air_tpg:
    gamma_p = 1.4
    Cp_p = 1006
    Cv_p = 718
    R_p = Cp_p - Cv_p
    theta = 3055.556

    def Cp_fn( k_p, Cp_p, th, T ):
        import numpy as np
        Cp = Cp_p * (1+(k_p-1)/k_p * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cp

    def Cv_fn( k_p, Cv_p, th, T ):
        import numpy as np
        Cv = Cv_p * (1+(k_p-1)/k_p * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cv

    def R_fn( Cp, Cv ):
        R = Cp - Cv
        return R

    def gamma_fn( Cp, Cv ):
        gamma = Cp / Cv
        return gamma