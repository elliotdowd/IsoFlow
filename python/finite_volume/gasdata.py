from python.finite_volume.helper import thermo
import numpy as np

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

###########################################################################
## mixture component class, input chemical formula, fifth order fit coefficient vector for Kp/p, mole fraction of equilibrium mixture
##                                characteristic vibration temperature, characteristic temperature of dissociation
###########################################################################	

class mix_component:

    def __init__( self, formula, poly_fit, molefrac, theta, Dp ):

        self.formula = formula
        self.poly_fit = poly_fit
        self.b = molefrac
        self.tc = theta
        self.D = Dp

    # calculate equilibrium coefficient
    def Kp_fn( self, T ):

        poly = self.poly_fit
        N = len(poly)

        Kp = 0

        for i in range(0, len(poly)):
            Kp = Kp + poly[i]*T**((N-1)-i)

        # fit is logarithmic, return arithmetic
        Kp = 10 ** Kp

        return Kp

    # calculate fraction of dissociated gas
    def diss_frac( self, Kp, p ):

        # fitting functions for finding molar fraction of dissociated gas
        alpha = lambda x: (4762*np.sqrt(x))/(1881*np.sqrt(x) + np.sqrt(19048000 + 8300161*x))
        beta = lambda x:  (7658*np.sqrt(x))/(1329*np.sqrt(x) + np.sqrt(153160000 + 40056241*x))
        
        # convert pressure to atm
        p = p * 9.86923e-6

        if self.formula == 'O2':
            a = alpha(Kp/p)
        elif self.formula == 'N2':
            a = beta(Kp/p)

        return a


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

    def ht_fn( self, state ):

        ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, self.gamma_fn(self.Cp, self.Cv)) / \
                                    state.Q[:,:,0] + state.p/state.Q[:,:,0]

        return ht

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

    def ht_fn( self, state ):

        ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, self.gamma_fn(self.Cp, self.Cv)) / \
                                    state.Q[:,:,0] + state.p/state.Q[:,:,0]

        return ht

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

    def ht_fn( self, state ):

        ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, self.gamma_fn(self.Cp, self.Cv)) / \
                                    state.Q[:,:,0] + state.p/state.Q[:,:,0]

        return ht

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

    mu_p = 1.803e-5

    # set up mixture class list
    mix = []
    mix.append( mix_component( 'O2', np.array([1.731E-10, -2.678E-06, 1.472E-02, -2.663E+01]), 0.21, 2270, 59000 ) )
    mix.append( mix_component( 'N2', np.array([1.388E-18, -0.00000000000005464, 0.0000000008668, -0.000006996, 0.02992, -56.08]), 0.79, 3390, 113200 ) )

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

    def ht_fn( self, state ):
        # self.option = 'dissociation'
        if self.option == 'perfect':
            ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, self.gamma_fn(self.Cp, self.Cv)) / \
                                        state.Q[:,:,0] + state.p/state.Q[:,:,0]
        elif self.option == 'dissociation':
            O2 = self.mix[0]
            N2 = self.mix[1]
            # determine equilibrium constants divided by pressure, account for exponent overflow
            Kpalpha_p = O2.Kp_fn(state.T)
            Kpalpha_p[np.isinf(Kpalpha_p)] = 1E10
            Kpbeta_p = N2.Kp_fn(state.T)
            Kpbeta_p[np.isinf(Kpbeta_p)] = 1E10

            # determine mole fractions of dissociated gas
            state.O2alpha = O2.diss_frac(Kpalpha_p, state.p)
            state.N2alpha = N2.diss_frac(Kpbeta_p, state.p)

            # modified version of equation 3.43 from Hermann 1965 (UARI report no. 30)
            # third line replaces (7/2)T term with pressure based calculation for better match with rho*et enthalpy calculation
            e = self.R_fn(self.Cp,self.Cv)*( O2.b*state.O2alpha*O2.D + N2.b*state.N2alpha*N2.D + \
                                             (1-state.O2alpha)*O2.b*(O2.tc/((np.exp(O2.tc/state.T)-1))) + \
                                             (N2.b)*(1-state.N2alpha)*(N2.tc/((np.exp(N2.tc/state.T)-1))) )

            et = e + thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, self.gamma_fn(self.Cp, self.Cv))/state.Q[:,:,0]
                                              

            # h = self.R_fn(self.Cp,self.Cv)*( O2.b*state.O2alpha*O2.D + N2.b*state.N2alpha*N2.D + \
            #                                 (3/2)*O2.b*state.O2alpha*state.T + (3/2)*N2.b*state.N2alpha*state.T + \
            #                                 # state.p/((self.gamma_fn(self.Cp,self.Cv)-1)*state.Q[:,:,0]) / self.R_fn(self.Cp,self.Cv) + \
            #                                 (7/2)*state.T + \
            #                                 (1-state.O2alpha)*O2.b*(O2.tc/((np.exp(O2.tc/state.T)-1))) + \
            #                                 (N2.b)*(1-state.N2alpha)*(N2.tc/((np.exp(N2.tc/state.T)-1))) )

            h = e + self.R_fn(self.Cp,self.Cv)*( (3/2)*O2.b*state.O2alpha*state.T + (3/2)*N2.b*state.N2alpha*state.T + \
                                                 (7/2)*state.T )

            ht = h + (state.u**2 + state.v**2)/2
            # state.Q[:,:,3] = state.Q[:,:,0] * et

        return ht

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

    def ht_fn( self, state ):

        ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, self.gamma_fn(self.Cp, self.Cv)) / \
                                    state.Q[:,:,0] + state.p/state.Q[:,:,0]

        return ht

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

    def ht_fn( self, state ):

        ht = thermo.calc_rho_et(state.p, state.Q[:,:,0], state.u, state.v, self.gamma_fn(self.Cp, self.Cv)) / \
                                    state.Q[:,:,0] + state.p/state.Q[:,:,0]

        return ht
        
