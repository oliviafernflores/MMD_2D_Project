import dadi
import matplotlib.pyplot as plt
import numpy as np


def plot_im_demography(fs, popt):
    demo_model = dadi.Demographics2D.IM
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    