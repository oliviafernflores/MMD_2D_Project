import dadi
import numpy as np
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

def IM(params, ns, pts):
    """
    ns = (n1,n2)
    params = (s,nu1,nu2,T,m12,m21, F)

    Isolation-with-migration model with exponential pop growth.

    s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    s,nu1,nu2,T,m12,m21, F = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F, F), (2, 2))
    return fs