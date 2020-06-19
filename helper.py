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


# measure cpu time w/ tic toc
class cpu_time:
    def TicTocGenerator():
        # Generator that returns time differences
        import time
        ti = 0           # initial time
        tf = time.time() # final time
        while True:
            ti = tf
            tf = time.time()
            yield tf-ti # returns the time difference

    TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

    # This will be the main function through which we define both tic() and toc()
    def toc(tempBool=True):
        # Prints the time difference yielded by generator instance TicToc
        tempTimeInterval = next(TicToc)
        if tempBool:
            print( "init_state time: %f seconds.\n" %tempTimeInterval )

    def tic():
        # Records a time in TicToc, marks the beginning of a time interval
        toc(False)