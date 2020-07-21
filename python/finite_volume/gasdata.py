###########################################################################
## References
###########################################################################	

# Cp and Cv thermally perfect relations from NACA Report 1135

###########################################################################
## Units
###########################################################################	

# Cp, Cv, R ~ J/kg*K
# Theta ~ K (molecular vibrational energy constant)

###########################################################################
## Calorically perfect gases (constant Cp, Cv)
###########################################################################	

# calorically perfect air
class air_cpg:
    gamma_p = 1
    Cp_p = 1006
    Cv_p = 718
    R_p = 287
    theta = 3055.556

    mu_p = 1.803e-5

    def Cp_fn( k_p, Cp_p, th, T ):
        import numpy as np
        Cp = Cp_p * (1+(k_p-1)/k_p * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cp

    def Cv_fn( k_p, Cv_p, th, T ):
        import numpy as np
        Cv = Cv_p * (1+(k_p-1) * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cv

    def R_fn( Cp, Cv ):
        R = Cp - Cv
        return R

    def gamma_fn( Cp, Cv ):
        gamma = Cp / Cv
        return gamma

    def mu_fn( T ):

        T0 = 273.11
        S = 110.56
        mu = air_cpg.mu_p * ( (T/T0)**(3/2) ) * ( (T0+S)/(T+S) )

        return mu

# calorically perfect CO2
class C02_cpg:
    gamma_p = 1
    Cp_p = 849
    Cv_p = 658
    R_p = Cp_p - Cv_p
    theta = 3055.556

    def Cp_fn( k_p, Cp_p, th, T ):
        import numpy as np
        Cp = Cp_p * (1+(k_p-1)/k_p * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cp

    def Cv_fn( k_p, Cv_p, th, T ):
        import numpy as np
        Cv = Cv_p * (1+(k_p-1) * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cv

    def R_fn( Cp, Cv ):
        R = Cp - Cv
        return R

    def gamma_fn( Cp, Cv ):
        gamma = Cp / Cv
        return gamma

# calorically perfect H2
class H2_cpg:
    gamma_p = 1
    Cp_p = 14310
    Cv_p = 10184
    R_p = Cp_p - Cv_p
    theta = 3055.556

    def Cp_fn( k_p, Cp_p, th, T ):
        import numpy as np
        Cp = Cp_p * (1+(k_p-1)/k_p * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cp

    def Cv_fn( k_p, Cv_p, th, T ):
        import numpy as np
        Cv = Cv_p * (1+(k_p-1) * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
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
    R_p = 287
    theta = 3055.556

    def Cp_fn( k_p, Cp_p, th, T ):
        import numpy as np
        Cp = Cp_p * (1+(k_p-1)/k_p * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cp

    def Cv_fn( k_p, Cv_p, th, T ):
        import numpy as np
        Cv = Cv_p * (1+(k_p-1) * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cv

    def R_fn( Cp, Cv ):
        R = Cp - Cv
        return R

    def gamma_fn( Cp, Cv ):
        gamma = Cp / Cv
        return gamma

# thermally perfect CO2
class C02_tpg:
    gamma_p = 1.29
    Cp_p = 849
    Cv_p = 658
    R_p = Cp_p - Cv_p
    theta = 2340

    def Cp_fn( k_p, Cp_p, th, T ):
        import numpy as np
        Cp = Cp_p * (1+(k_p-1)/k_p * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cp

    def Cv_fn( k_p, Cv_p, th, T ):
        import numpy as np
        Cv = Cv_p * (1+(k_p-1) * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cv

    def R_fn( Cp, Cv ):
        R = Cp - Cv
        return R

    def gamma_fn( Cp, Cv ):
        gamma = Cp / Cv
        return gamma

# thermally perfect H2
class H2_tpg:
    gamma_p = 1.405
    Cp_p = 14310
    Cv_p = 10184
    R_p = Cp_p - Cv_p
    theta = 4960

    def Cp_fn( k_p, Cp_p, th, T ):
        import numpy as np
        Cp = Cp_p * (1+(k_p-1)/k_p * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cp

    def Cv_fn( k_p, Cv_p, th, T ):
        import numpy as np
        Cv = Cv_p * (1+(k_p-1) * ((th/T)**2*np.exp(th/T)/(np.exp(th/T)-1)**2))
        return Cv

    def R_fn( Cp, Cv ):
        R = Cp - Cv
        return R

    def gamma_fn( Cp, Cv ):
        gamma = Cp / Cv
        return gamma
