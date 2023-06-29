from contextlib import AsyncExitStack
import dadi
import pickle
import nlopt
import matplotlib.pyplot as plt
import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum
def bottlegrowth_split_mig(params, ns, pts):
    nuB,nuF,m,T,Ts, F = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    if T >= Ts:
        nu_func = lambda t: nuB*numpy.exp(numpy.log(nuF/nuB) * t/T)
        phi = Integration.one_pop(phi, xx, T-Ts, nu_func)

        phi = PhiManip.phi_1D_to_2D(xx, phi)
        nu0 = nu_func(T-Ts)
        nu_func = lambda t: nu0*numpy.exp(numpy.log(nuF/nu0) * t/Ts)
        phi = Integration.two_pops(phi, xx, Ts, nu_func, nu_func, m12=m, m21=m)
    else:
        phi = PhiManip.phi_1D_to_2D(xx, phi)
        phi = Integration.two_pops(phi, xx, Ts-T, 1, 1, m12=m, m21=m)
        nu_func = lambda t: nuB*numpy.exp(numpy.log(nuF/nuB) * t/T)
        phi = Integration.two_pops(phi, xx, T, nu_func, nu_func, m12=m, m21=m)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F, 1e-10), (2, 2))
    return fs
def bottlegrowth_split(params, ns, pts):
    nuB,nuF,T,Ts, F = params
    return bottlegrowth_split_mig((nuB,nuF,0,T,Ts, F), ns, pts)
def IM(params, ns, pts):
    s,nu1,nu2,T,m12,m21, F = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F, 1e-10), (2, 2))
    return fs
def IM_pre(params, ns, pts):
    nuPre,TPre,s,nu1,nu2,T,m12,m21, F = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, TPre, nu=nuPre)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1_0 = nuPre*s
    nu2_0 = nuPre*(1-s)
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F, 1e-10), (2, 2))
    return fs
def split_asym_mig(params, ns, pts):
    nu1,nu2,T,m12,m21, F = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F, 1e-10), (2, 2))
    return fs
def split_delay_mig(params, ns, pts):
    nu1,nu2,Tpre,Tmig,m12,m21, F = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Tpre, nu1, nu2, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, Tmig, nu1, nu2, m12=m12, m21=m21)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F, 1e-10), (2, 2))
    return fs
def split_mig(params, ns, pts):
    nu1,nu2,T,m, F = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F, 1e-10), (2, 2))
    return fs