eps = 1
kap = 1

def QL( Qm2, Qm1, Qi, Qp1, Qp2, kap, eps, limiter ):

    rL = lambda Qm2, Qm1, Qi: (Qi-Qm1) / (Qm1-Qm2)
    QL = Qm1 + (eps/4) * ( (1-kap)*(Qm1-Qm2)*limiter(rL(Qm2, Qm1, Qi)) + \
                           (1+kap)*(Qi-Qm1)*limiter(1/rL(Qm2, Qm1, Qi)) )
    return QL

def QR( Qm1, Qi, Qp1, kap, eps, limiter ):

    rR = lambda Qm1, Qi, Qp1: (Qi-Qm1) / (Qp1-Qi)
    QR = Qi - (eps/4) * ( (1+kap)*(Qi-Qm1)*limiter(1/rR(Qm1, Qi, Qp1)) + \
                          (1-kap)*(Qp1-Qi)*limiter(rR(Qm1, Qi, Qp1)) )
    return QR

