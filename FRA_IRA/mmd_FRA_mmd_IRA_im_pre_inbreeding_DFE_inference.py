import pickle, random
import numpy as np
import nlopt
import dadi
import dadi.DFE as DFE
import matplotlib.pyplot as plt
from dadi.DFE import *

def dfe_1():
    demo_params = [2.3892200757761812, 0.09865800019618395, 0.059546476404449086, 0.29863880148384436, 2.9956042751467993, 0.1765401593467742, 2.1701288461162407, 0.19321473756551447]
    theta_0 = 26681.214011656422
    theta_ns = theta_0 * 2.4
    ns = [16, 16]