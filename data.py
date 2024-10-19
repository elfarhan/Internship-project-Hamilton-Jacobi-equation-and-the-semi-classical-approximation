# coding: utf-8
from library import *
import numpy as np
import sympy as s
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib import interactive
from matplotlib.animation import FuncAnimation
from scipy.stats import norm
import math
import multiprocessing as multiproc
import datetime
from joblib import Parallel, delayed
from numba import jit
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
from IPython.display import display, Math
# symmetric two wells
b = 30
c = 1600
m = -0
a_4 = 1*.2
a_3 = -5*.2
a_2 = -3*.2
a_1 = -5*.2
a_0 = -4*.2

modelParameters = [1,1,0.,1,1.,1,a_1,a_2,a_3,a_4, a_0]
parameterNames = ["dimension", "m",'\mu', '\lambda', 'p_0','hbar','a_1','a_2','a_3','a_4','a_0']

one_zero_params = parameters(modelParameters, parameterNames)
initConds = initialConditions(W,p_0, one_zero_params)
one_zero = PhysicalSystem(parameters=one_zero_params,
                                     initialConditions=initConds,
                                     Hamiltonian= asymmetric_Oscillator_H,
                                    HamiltonEquations=asymmetric_Oscillator_H_eq)
from library import *
import numpy as np
import sympy as s
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib import interactive
from matplotlib.animation import FuncAnimation
from scipy.stats import norm
import math
import multiprocessing as multiproc
import datetime
from joblib import Parallel, delayed
from numba import jit
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
from IPython.display import display, Math
b = 30
c = 1600
m = -0
a_4 = 1*.2
a_3 = -5*.2
a_2 = -3*.2
a_1 = -5*.2
a_0 = -4*.2

modelParameters = [1,1,0.,1,1.,1,a_1,a_2,a_3,a_4, a_0]
parameterNames = ["dimension", "m",'\mu', '\lambda', 'p_0','hbar','a_1','a_2','a_3','a_4','a_0']

one_zero_params = parameters(modelParameters, parameterNames)
initConds = initialConditions(W,p_0, one_zero_params)
one_zero = PhysicalSystem(parameters=one_zero_params,
                                     initialConditions=initConds,
                                     Hamiltonian= asymmetric_Oscillator_H,
                                    HamiltonEquations=asymmetric_Oscillator_H_eq)
ONEZERO = WKB(one_zero, N = 1000, method = 'count', dx = 0.2)
one_zerotrajectories, one_zeropsi_WKB = ONEZERO.psi_eval()
# ASHO
minimum_location = 2

a_4 = 1
a_2 = -2*minimum_location**2
a_1 = 0
a_3 = 0
a_0 = minimum_location**4

modelParameters = [1,1.,0.,1,1.,1,a_1,a_2,a_3,a_4, a_0]
parameterNames = ["dimension", "m",'\mu', '\lambda', 'p_0','hbar','a_1','a_2','a_3','a_4','a_0']
symmetric_two_wells_params = parameters(modelParameters, parameterNames)
symmetric_two_wells = PhysicalSystem(parameters=symmetric_two_wells_params,
                                     initialConditions=initConds,
                                     Hamiltonian= asymmetric_Oscillator_H,
                                    HamiltonEquations=asymmetric_Oscillator_H_eq)
ASHO = WKB(symmetric_two_wells, N = 1000, method = 'count', dx = 0.2)
trajectoriesASHO, ASHOpsi_WKB = ASHO.psi_eval()
minimum_location = 2.5

a_4 = 1
a_3 = -2*minimum_location
a_2 = 0
a_1 = 2*minimum_location**3
a_0 = -minimum_location**4

modelParameters = [1,1,0.,1,1.,1,a_1,a_2,a_3,a_4, a_0]
parameterNames = ["dimension", "m",'\mu', '\lambda', 'p_0','hbar','a_1','a_2','a_3','a_4','a_0']
deeper_well_params = parameters(modelParameters, parameterNames)
deeper_well = PhysicalSystem(parameters=deeper_well_params,
                                     initialConditions=initConds,
                                     Hamiltonian= asymmetric_Oscillator_H,
                                    HamiltonEquations=asymmetric_Oscillator_H_eq)
DW = WKB(deeper_well,N = 1000, method = 'count', dx = 0.2)
deeper_welltrajectories, deeper_wellpsi_WKB = DW.psi_eval()
