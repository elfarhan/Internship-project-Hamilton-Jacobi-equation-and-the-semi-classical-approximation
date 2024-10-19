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

class numpy_to_sympy:
    def __enter__(self):
        self.original_np = np.exp
        np.exp = s.exp
        # Add more mappings as needed
        return self
    def __exit__(self, type, value, traceback):
        np.exp = self.original_np
        # Restore other mappings as needed


class parameters():
    def __init__(self, modelParameters, parameterNames):
        self.modelParameters = modelParameters
        self.parameterNames = parameterNames
        return
    def show(self):
        print("Parameters:")
        nameSpace = [s.symbols(x) for x in self.parameterNames]
        for i in range(len(nameSpace)):
            display(Math(f'{s.latex(nameSpace[i])}'+" = "+str(self.modelParameters[i])))
        return
    def parseParameter(self, p):
        for i in range(len(self.parameterNames)):
            if self.parameterNames[i]==p:
                return self.modelParameters[i]
        print("parameter "+ p+ "is not defined")
        retrun 

class initialConditions():
    def __init__(self, W, p_0, parameters, method = "track"):
        self.W = W
        self.p_0 = p_0
        self.parameters = parameters
        self.dimension = parameters.modelParameters[0]
        #with numpy_to_sympy():
        #    self.show()
        self.W = lambda q: W(self.parameters, q)
        self.p_0 = lambda q: p_0(self.parameters, q)
        self.method = method
        #self.q_0, self.p_0 = self.initilze_particles()
        return 
    def initilze_particles(self):
        ## assume initial distribution of particles is a gaussian
        self.l = self.parameters.modelParameters[4]
        if self.method == "track":
            if self.dimension == 1:
                self.x_max = 3*1/(2*self.l) # interval boundary of the simulation
                x_max = 3*self.x_max
                q_0 = np.concatenate((np.linspace(-x_max,-self.x_max,int(1/10*self.N)),\
                                  np.linspace(-self.x_max,self.x_max, int(4/5*self.N)),\
                                 np.linspace(self.x_max,x_max, int(1/10*self.N))))
                self.q_0 = q_0
            
            elif self.dimension == 2:
                self.x_max = 3*1/(2*max(self.l))
                x_max = 3*self.x_max
                N = int(np.sqrt(N))
                x_arr = np.concatenate((np.linspace(-x_max,-self.x_max,int(1/10*N)),\
                                  np.linspace(-self.x_max,self.x_max, int(4/5*N)),\
                                 np.linspace(self.x_max,x_max, int(1/10*N))))
                q_0 = []
                for x in x_arr:
                    for y in x_arr:
                        q_0.append(np.array([x,y]))
        if self.method == "count":
            if self.dimension == 1:
                self.x_max = 30*1/(2*self.l)
                n = 0
                q_0 = [] 
            
                while n<self.N:
                
                    # use Von Neumann rejection sampling
                    x = np.random.uniform(-self.x_max*15, self.x_max*15)  # Generate a random value in the desired range
                    y = np.random.uniform(0, 1)  # Generate a random value between 0 and the maximum probability
                
                    if y <= self.W(x):
                        n +=1
                        q_0.append(x)
            elif self.dimension == 2:
                self.x1_max = 3*1/(2*self.l[0])
                self.x2_max = 3*1/(2*self.l[1])
                n = 0
                q_0 = [] 
                while n<self.N:
                    x1 = np.random.uniform(-self.x1_max*15, self.x1_max*15)
                    x2 = np.random.uniform(-self.x2_max*15, self.x2_max*15)
                    y = np.random.uniform(0, 1)
                    if y <= self.W(np.array([x1,x2])):
                            n +=1
                            q_0.append(np.array([x1,x2]))
                    
            p_0 = [self.p_0(q) for q in q_0]
        return np.array(q_0), np.array(p_0)
    
    def show(self):
        nameSpace = [s.symbols(x) for x in self.parameters.parameterNames]
        show_parameters = parameters(nameSpace, self.parameters.parameterNames)
        q = s.symbols("q")
        show_parameters.modelParameters[0]=1
        print("Initial Cnditions:")
        display(Math(f'{s.latex("W(q)")}' + ' = ' + f'{s.latex(self.W(show_parameters, q))}'))
        display(Math("P_0(q)" + ' = ' + f'{s.latex(self.p_0(show_parameters, q))}'))

class PhysicalSystem():
    """ a class to handle the Hamiltonian, the Hamilton equations
       as well as the model parameters"""
    def __init__(self, Hamiltonian, HamiltonEquations, initialConditions, parameters):
        """Args: 
            - Hamiltonian: model Hamiltonian args: Hamiltonian(modelParameters, p, q)
            - HamiltonEquations: classical equations args: HamiltonEquations(modelParameters, p, q)
            both functions take model_params as only parameter
            """
        self.parameters = parameters
        self.initialConditions = initialConditions
        self.HamiltonEquations = HamiltonEquations
        self.Hamiltonian = Hamiltonian
        self.show()
        self.HamiltonEquations = lambda p,q: HamiltonEquations(self.parameters, p, q)
        self.Hamiltonian = lambda p,q: Hamiltonian(self.parameters, p, q)
        #self.plot_potential()
        return
    def plot_potential(self):
        q = np.linspace(-6,6)
        V = np.array([self.Hamiltonian(p=0,q=x) for x in q])
        plt.close('all')
        plt.plot(q,V)
        plt.title("The potential V(q)")
        plt.xlabel(r'$q$')
        plt.ylabel(r'$V(q)$')
        plt.grid()
        plt.show()
        return
    def show(self):
        aa = s.symbols("dt")
        aaa = s.symbols("dq")
        aaaa = s.symbols("dp")
        frac_p = aaaa/aa
        frac_q = aaa/aa
        p = s.symbols("p")
        q = s.symbols("q")
        nameSpace = [s.symbols(x) for x in self.parameters.parameterNames]
        show_parameters = parameters(nameSpace, self.parameters.parameterNames)
        show_parameters.modelParameters[0]=1
        print("Hamiltonian:")
        display(Math('H(p,q) = '+f'{s.latex(self.Hamiltonian(show_parameters, p, q))}'))
        print("Hamilton equations:")
        display(Math(f'{s.latex(frac_q)}' + ' = ' + f'{s.latex(self.HamiltonEquations(show_parameters, p, q)[0])}'))      
        display(Math(f'{s.latex(frac_p)}' + ' = ' + f'{s.latex(self.HamiltonEquations(show_parameters, p, q)[1])}')) 
        self.parameters.show()
    

class WKB():
    def __init__(self, PhysicalSystem, N = 1000, dt = 0.01, t_max = 20, dx = 0.01, method = "track"):
        
        # inherit PhysicalSystem
        self.PhysicalSystem = PhysicalSystem
        self.parameters = PhysicalSystem.parameters
        self.initialConditions = PhysicalSystem.initialConditions
        self.hamilton_equations = PhysicalSystem.HamiltonEquations
        self.H = PhysicalSystem.Hamiltonian

        # define integration parameters
        self.dx = dx
        self.dimension = self.parameters.modelParameters[0]
        self.h = self.parameters.parseParameter("hbar") # define hbar
        self.dt = dt
        self.num_steps = int(t_max/dt)
        self.t = [i*dt for i in range(self.num_steps)]
        
        # initialize a set of particles
        PhysicalSystem.initialConditions.N = N
        self.N = N
        PhysicalSystem.initialConditions.method = method
        self.q_0, self.p_0 = PhysicalSystem.initialConditions.initilze_particles()
        self.trajectories = []
        return 
    
    
    
    def verlet_solver(self, q_0, p_0, dt, num_steps):
        """
        Verlet integrator to solve Hamilton's equations of motion.
        Args:
            q_0: Initial position
            p_0: Initial momentum
            dt: Time step size
            num_steps: Number of integration steps
        Returns:
            q_traj: Array of position values over time
            p_traj: Array of momentum values over time
            S: classical action of the trajectory
        """
        # set the mass 
        m = self.parameters.parseParameter("m")
        
        # Initialize arrays to store trajectory and action
        S = np.zeros(num_steps)
        if self.dimension == 1:
            
            q = np.zeros(num_steps)
            p = np.zeros(num_steps)
            # Set initial conditions
            q[0] = q_0
            p[0] = p_0 
            S[0] = p_0*q[0]
            
            # Perform the Verlet integration
            for i in range(1, num_steps):
        
                # Compute the half-step momentum
                p_half = p[i-1] + 0.5 * dt * self.hamilton_equations(p[i-1], q[i-1])[1]

                # Compute the new position
                q[i] = q[i-1] + dt * p_half/m

                # Compute the new momentum
                p[i] = p_half + 0.5 * dt * self.hamilton_equations(p_half, q[i])[1]
        
                # Compute the action
                S[i] = S[i-1] + (p[i]**2/m - self.H(p[i],q[i]))*dt
                
        elif self.dimension == 2:
            
            q = np.array([np.array([0.0,0.0]) for _ in range(num_steps)])
            p = np.array([np.array([0.0,0.0]) for _ in range(num_steps)])
            # Set initial conditions
            q[0] = q_0
            p[0] = p_0 
            S[0] = p_0@q[0]
            # Perform the Verlet integration
            for i in range(1, num_steps):
        
                # Compute the half-step momentum
                p_half = p[i-1] + 0.5 * dt * self.hamilton_equations(p[i-1],q[i-1])[1]

                # Compute the new position
                q[i] = q[i-1] + dt * p_half/m

                # Compute the new momentum
                p[i] = p_half + 0.5 * dt * self.hamilton_equations(p_half, q[i])[1]
        
                # Compute the action
                S[i] = S[i-1] + (p[i]@p[i]/m - self.H(p[i],q[i]))*dt
        
            
        return q, p, S
    
    
    def compute_a_trajectory(self, i):
        x_0 = self.q_0[i]
        p_0 = self.p_0[i]
        q, p, s = self.verlet_solver(x_0, p_0, self.dt, self.num_steps)
        return q, s
        
    def find_trajectories(self):
        """
            a function to calculate the trajectories of the particles starting at
            the ensemble self.q_0 over the time t
            Args:
                self.q_0: intial positions
            Returns:
                self.trajectories: trajectories[position,time]
                self.S: action of the trajectories S[position,time]
        """

        # Parallel computation
        results = Parallel(n_jobs=multiproc.cpu_count()-1)(delayed(self.compute_a_trajectory)(i) for i in range(len(self.q_0)))
            
        # Extract results
        trajectories = [res[0] for res in results]
        S = [res[1] for res in results]
        self.S = np.array(S)
        self.trajectories = np.array(trajectories)
        return
    
    def find_densities(self):
        """
        a function to calculate the wave function in the WKB approximation
        Args:
            self.trajectories: list of q[x][t]
        Returns:
            rho: rho[t,q[t]]
        """
        if self.PhysicalSystem.initialConditions.method == "track":
            
            self.valid_index  = []
            N = self.N
            for j in range(int(np.sqrt(N))**2):
                if (j+1)%int(np.sqrt(N)) != 0 and (j)%int(np.sqrt(N)) != 0 and j-int(np.sqrt(N)) >= 0\
                    and (j+int(np.sqrt(N))) < int(np.sqrt(N))**2 and self.dimension ==2:
                    self.valid_index.append(j)
            if self.dimension==1:
                self.valid_index = [j for j in range(1,len(self.trajectories[:,0])-1)]
            
            densities = Parallel(n_jobs=multiproc.cpu_count()-1)(delayed(rho_track)\
                                (j, time, self.trajectories, self.initialConditions, self.dimension)
                                for time in range(self.num_steps) for j in self.valid_index )  
                
            self.densities = np.array(densities)
            self.densities.shape = (self.num_steps,len(self.valid_index))
            
        elif self.PhysicalSystem.initialConditions.method == "count":
            
            self.valid_index  = [j for j in range(len(self.trajectories[:,0]))]
            densities = Parallel(n_jobs=multiproc.cpu_count()-1)(delayed(density)\
                              (position, time, self.dx, self.trajectories)
                              for time in range(self.num_steps) for position in self.trajectories[:,time] ) 
        
            self.densities = np.array(densities)
            self.densities.shape = (self.num_steps,len(self.trajectories[:,0]))
        return

    def psi_eval(self):
        """
        a function to calculate the wave function in the WKB approximation
        Args:
            x: position
            t: time is an intger
            self.dx: volume element size (side length)
            self.trajectories: list of q[x][t]
        Returns:
            psi_WKB: psi[q[t],t]
            trajectories[q[t],t]
            
        """
        if self.PhysicalSystem.initialConditions.method == "count":
            #solve e.o.m. and find particle density
            self.find_trajectories()
            self.find_densities()
        
            #intialize output and used variables
            S = self.S
            trajectories = self.trajectories
            densities = self.densities
            psi_WKB = np.zeros((self.num_steps,len(trajectories[:])), dtype = complex)
        
            # compute the wave function
            for time in range(self.num_steps):
            
                for position in range(len(trajectories[:,time])):
                
                    psi_WKB[time,position] = np.sqrt(densities[time,position])*np.exp(1j/self.h*S[position,time])
            
            self.psi_WKB = psi_WKB
        
            return self.trajectories, psi_WKB
        elif self.PhysicalSystem.initialConditions.method == "track":
            #solve e.o.m. and find particle density
                
            self.find_trajectories()
            self.find_densities()
            
            #intialize output and used variables
            S = np.array([self.S[j] for j in self.valid_index])
            trajectories = np.array([self.trajectories[j] for j in self.valid_index])
            densities = self.densities
            psi_WKB = np.zeros((self.num_steps,len(trajectories[:])), dtype = complex)
            
            # compute the wave function
            for time in range(self.num_steps):
                
                for position in range(len(trajectories[:,time])):
                    
                    psi_WKB[time,position] = densities[time,position]#np.sqrt(densities[time,position])*np.exp(1j/self.h*S[position,time])
            return trajectories, psi_WKB
    

def rho_track(j, time, trajectories, initialConditions, dimension):
    # find W
    x = trajectories[j,0]
    w = initialConditions.W(x)
    #find area at t = 0
    q = trajectories[:,0]
    if dimension == 2:
        x1,y1 = q[above(j,N),0],q[above(j,N),1]
        x2,y2 = q[left(j,N),0],q[left(j,N),1]
        x3,y3 = q[below(j,N),0],q[below(j,N),1]
        x4,y4 = q[right(j,N),0],q[right(j,N),1]
        A1 =(x1*y2 + x2*y3 + x3*y4 + x4*y1 )  
        A2 = x2*y1 + x3*y2 + x4*y3 + x1*y4
        area_0 = 1/2 *(A1-A2)
    if dimension ==1:
        area_0 = abs(q[j+1]-q[j-1])

    #find area at t = time
    q = trajectories[:,time]
    if dimension == 2:
        x1,y1 = q[above(j,N),0],q[above(j,N),1]
        x2,y2 = q[left(j,N),0],q[left(j,N),1]
        x3,y3 = q[below(j,N),0],q[below(j,N),1]
        x4,y4 = q[right(j,N),0],q[right(j,N),1]
        A1 =(x1*y2 + x2*y3 + x3*y4 + x4*y1 )  
        A2 = x2*y1 + x3*y2 + x4*y3 + x1*y4
        area_t = 1/2 *(A1-A2)
    if dimension ==1:
        area_t = abs(q[j+1]-q[j-1])
                     
    return w * area_0 / area_t
@jit
def density(x, t, dx, trajectories):
        """
        method ='count'
        a function to calculate the number of particles density
        Args:
            x: position
            t: time is an intger
            dx: volume element size (side length)
            trajectories: list of q[x][t]
        Returns:
            rho: particle number density
        """
        # set initial number of particles
        N = 0
    
        for q in trajectories:
        
            #compute the number of particles in a box with side-length dx centred around x 
            position = q[t]  
            
            if np.linalg.norm(position-x)<=dx:
                
                N +=1
                
        # calculate the density
        rho = N/len(trajectories[:])/(2*dx)
    
        return rho

def psi_analytical(x,t, parameters):
    """
    calculates the wkb approximation of the wave function acccording
    to the analytical derivartion
    Args:
         x: position
         t: time
         rest are model parameters
    Returns:
         psi_WKB(x,t)
    """
    p_0=parameters.parseParameter("p_0")
    m=parameters.parseParameter("m")
    w=parameters.parseParameter("\omega")
    l=parameters.parseParameter("\lambda")
    mu=parameters.parseParameter("\mu")
    h=parameters.parseParameter("hbar")
    dimension = parameters.modelParameters[0]
    if dimension ==1:
        Z = np.sqrt(1/np.sqrt(np.pi/2) * l/abs(np.cos(w*t)))
        re = -l**2*(x/np.cos(w*t)-p_0/(m*w)*np.tan(w*t)-mu)**2
        im = 1j/h *(p_0*(x/np.cos(w*t)-p_0/m/w*np.tan(w*t)) + np.tan(w*t)*(p_0**2/(2*m*w)-1/2*m*w*x**2))
        psi = Z * np.exp(re+im)
    if dimension >1:
        Z = 1
        im = 1j*0.0
        re = 0.0
        for i in range(len(w)):
            Z = Z * np.sqrt(l[i]/np.sqrt(np.pi/2) /abs(np.cos(w[i]*t)))
            im += 1j/h *(p_0[i]*(x[i]/np.cos(w[i]*t)-p_0[i]/m/w[i]*np.tan(w[i]*t))\
                         + np.tan(w[i]*t)*(p_0[i]**2/(2*m*w[i])-1/2*m*w[i]*x[i]**2))
            re -= l[i]**2*(x[i]/np.cos(w[i]*t)-p_0[i]/(m*w[i])*np.tan(w[i]*t)-mu[i])**2
        psi = Z * np.exp(re+im)
    return psi

def psi_analytical_eval(x_arr, t_arr, parameters):
    """
    calculates the wkb approximation of the wave function acccording
    to the analytical derivartion
    Args:
         x: list of positions
         t: list of time
    Returns:
         psi_WKB[time,position]
    """
    dimension = parameters.modelParameters[0]
    
    if dimension==1:
        
        PSI = []
    
        for xx in x_arr:
        
            dummy = []
            for tt in t_arr:
            
                dummy.append(psi_analytical(xx,tt, parameters))
            
            PSI.append(np.array(dummy))
        
        PSI = np.array(PSI)
        
        return PSI
    if dimension == 2:
        
        resolution = len(x_arr)
        PSI_2D = np.zeros((resolution,resolution, len(t_arr)), dtype=complex)
        for k in range(len(t_arr)):
            for i in range(resolution):
                for j in range(resolution):
                    PSI_2D[i,j,k] = psi_analytical([x_arr[i],x_arr[j]],t_arr[k],parameters)
        
        return PSI_2D
        

@jit
def H(n, x):
    """
    Hermite polynomial of degree n
    Args:
        n: degree of polynomial
        x: where the polynomial is evaluated
    Returns:
        H_n(x)
     """
        
    if n == 0:
        return 1
        
    elif n == 1:
        return 2.0 * x
        
    else:
            
        H0 = 1
        H1 = 2.0 * x
            
        for k in range(2, n + 1):
            Hk = 2.0 * x * H1 - 2.0 * (k - 1) * H0
            H0 = H1
            H1 = Hk
        return Hk
class Harmonic_Oscillator_QM():
    
    def __init__(self, l, w, m, h, p_0, mu, n_cutoff, zero_energey = True):
        
        # parameters
        self.l = l
        self.w = w
        self.m = m
        self.h = h
        self.p_0 = p_0
        self.mu = mu
        self.n_cutoff = n_cutoff
        
        self.a = self.l**2+1/2*self.m*self.w/self.h
        self.b = 2*self.l*self.mu+1j*self.p_0/self.h
        self.c = self.l**2*self.mu**2
        self.alpha =  np.sqrt(self.m*self.w/self.h)
        self.s = (self.alpha**2-self.a)/self.a
        self.k = self.alpha*self.b/self.a
        self.zero_energey = zero_energey
        self.G(n_cutoff)
        
        return
    
        
    def G(self, n_cutoff):
        """
        Helping function to evaluate scalar product <n|psi>
        n_cutoff: degree of H_n at which the infite sum is truncated
        """
        s = self.s
        k = self.k
        t = 0
        self.G_t = []
        self.G_t.append(np.sqrt(np.pi/self.a)*np.exp(s*t**2+k*t))
        self.G_t.append((2*s*t+k) * (np.sqrt(np.pi/self.a)*np.exp(s*t**2+k*t)))
        
        for n in range(2, n_cutoff+1):
            
            self.G_t.append(2*s*(n-1)*self.G_t[n-2] + (2*s*t+k)*self.G_t[n-1])
            
        return
    
    def Z(self, n):
        """
        normalization factor for |psi_n>
        Args:
            n: nth eigenstate of the Hamiltonian that need to be normalized
        """
        return ( self.m*self.w/(np.pi*self.h) )**(1/4)* 1/np.sqrt(2.0**n * math.factorial(n))
    
    def psi_dot_n(self, n):
        """
        scalar product <n|psi>
        Args:
            n: nth eigenstate of the Hamiltonian
        """
        return self.Z(n)*np.exp(-self.c + self.b**2/(4*self.a)) * self.G_t[n] * np.sqrt(self.l/np.sqrt(np.pi/2))
    
    def E(self,n):
        """
        Energy of the nth eigenstate
        Args:
            n: nth eigenstate of the Hamiltonian
        """
        if self.zero_energey:
            return self.h*self.w*(n+.5)
        else:
            return self.h*self.w*(n)#+.5)
    
    def psi_n(self, n, x):
        """
        nth eigenstate in real space representation
        Args:
            n: nth eigenstate of the Hamiltonian
            x: where wave function is evaluated
        """
        return self.Z(n)* H(n,np.sqrt(self.m*self.w/self.h)*x)*np.exp(-1/2 * self.m*self.w/self.h*x**2)
    
    def psi(self,x,t):
        """
        calculates exact time evolution of the wave function in quantum mechanics
        Args:
            x: position
            t: time
        Returns:
            psi_QM(x,t)
        """
        Sum = 0
        
        for n in range(0,self.n_cutoff):
            
            Sum += np.exp(-1j*self.E(n)*t)* self.psi_dot_n(n) * self.psi_n(n,x)
            
        return Sum
    
    def psi_eval(self, x_arr, t_arr):
        """
        calculates the wkb approximation of the wave function acccording
        to the analytical derivartion
        Args:
            x: list of positions
            t: list of time
        Returns:
            psi_QM[time,position]
        """
        PSI = Parallel(n_jobs=multiproc.cpu_count()-1)(delayed(self.psi)\
                              (xx,tt)
                              for xx in x_arr for tt in t_arr)
        PSI = np.array(PSI)
        PSI.shape = (len(x_arr),len(t_arr))
    
        return PSI

    

def harmonic_oscillator_H(Parameters, p, q):
    # modelParameters: 
    dimension = Parameters.modelParameters[0]
    m = Parameters.parseParameter('m')
    w = Parameters.parseParameter('\omega')
     
    p = np.array(p)
    q = np.array(q)
    if dimension == 1:
        return 1/(2*m)*p*p + 1/2*m*w**2*q**2
    else:
        return 1/(2*m)*np.dot(p, p) + 1/2*m*w**2@q**2
    
def harmonic_oscillator_H_eq(Parameters, p, q):
    """
        Hamilton's equations of motion for a simple harmonic oscillator.
        Args:
            q: Position coordinate
            p: Momentum coordinate
        Returns:
            dq_dt: Derivative of position
            dp_dt: Derivative of momentum
        """
    
    # modelParameters: 
    dimension = Parameters.modelParameters[0]
    m = Parameters.parseParameter('m')
    w = Parameters.parseParameter('\omega')
    
    # Compute the derivatives
    dq_dt = p/m
    dp_dt = -m*w**2*q
    return dq_dt, dp_dt

def W(Parameters, q):
    mu = Parameters.parseParameter('\mu')
    l = Parameters.parseParameter('\lambda')
    dimension = Parameters.modelParameters[0]
    if dimension == 2:
        return np.exp(-2*(q-mu)@np.diag(l)**2@(q-mu))*l[0]*l[1]/np.pi/2
    if dimension ==1:
        return np.exp(-2*l**2*(q-mu)**2)*l/np.sqrt(np.pi/2)

def p_0(Parameters, q):
    return Parameters.parseParameter('p_0')

def Harmonic_oscillator1D(p, mu,l, w, m, h, n_cutoff = 170, resolution = 300, zero_energey = True):
    x_max = np.sqrt((p/(m*w))**2+10**2)
    # implement WKB:
    modelParameters = [1,m,w,mu,l,p,h]
    parameterNames = ["dimension", "m",'\omega','\mu', '\lambda', 'p_0','hbar']
    params = parameters(modelParameters, parameterNames)
    initConds = initialConditions(W,p_0, params)
    HO2D = PhysicalSystem(harmonic_oscillator_H, harmonic_oscillator_H_eq, initConds, params)
    a = WKB(HO2D, N = 500, method = 'count', dx = 0.2)
    trajectories, psi_WKB = a.psi_eval()
    # parameters of integration 
    x_max = np.sqrt((p/(m*w))**2+10**2)
    x_arr = np.linspace(-x_max,x_max,resolution)
    t_arr = np.array(a.t)
    # implement WKB analytical
    psi_WKB_analytical = psi_analytical_eval(x_arr, t_arr, p, m, w, h, l, mu, 1)
    # implement QM 
    QM_ = Harmonic_Oscillator_QM(l, w, m, h, p, mu, n_cutoff, zero_energey)
    psi_QM = QM_.psi_eval(x_arr, t_arr)
    return trajectories, x_arr, psi_WKB, psi_WKB_analytical, psi_QM

def Harmonic_oscillator2D(p, mu,l, w, m, h, n_cutoff = 170, resolution = 300, zero_energey = True):
    # implement WKB:
    modelParameters = [2,m,w,mu,l,p,h]
    parameterNames = ["dimension", "m",'\omega','\mu', '\lambda', 'p_0','hbar']
    params = parameters(modelParameters, parameterNames)
    initConds = initialConditions(W,p_0, params)
    HO2D = PhysicalSystem(harmonic_oscillator_H, harmonic_oscillator_H_eq, initConds, params)
    a = WKB(HO2D, N = 500, method = 'count', dx = 0.2)
    trajectories, psi_WKB_2D = a.psi_eval()
    # implement QM
    x_max = np.sqrt((max(p)/(m*w[0]))**2+10**2)
    x_arr = np.linspace(-x_max,x_max,resolution)
    t_arr = np.array(a.t)
    QM_0 = Harmonic_Oscillator_QM(l[0], w[0], m, h, p[0], mu[0], n_cutoff, zero_energey)
    QM_1 = Harmonic_Oscillator_QM(l[1], w[1], m, h, p[1], mu[1], n_cutoff, zero_energey)
    psi_QM_0 = QM_0.psi_eval(x_arr, t_arr)
    psi_QM_1 = QM_1.psi_eval(x_arr, t_arr)
    PSI_QM_2D = np.zeros((resolution,resolution, len(t_arr)), dtype=complex)
    X, Y = np.meshgrid(x_arr, x_arr)       
    for i in range(len(psi_QM_0)):
        for j in range(len(psi_QM_1)):
                for k in range(len(t_arr)):
                    PSI_QM_2D[i,j,k] = psi_QM_0[i,k]*psi_QM_1[j,k]
    
    return trajectories, psi_WKB_2D, X, Y, PSI_QM_2D

def asymmetric_Oscillator_H(Parameters, p, q):
    # modelParameters: 
    dimension = Parameters.modelParameters[0]
    m = Parameters.parseParameter('m')
    a_1 = Parameters.parseParameter('a_1')
    a_2 = Parameters.parseParameter('a_2')
    a_3 = Parameters.parseParameter('a_3')
    a_4 = Parameters.parseParameter('a_4')
    a_0 = Parameters.parseParameter('a_0')
    p = np.array(p)
    q = np.array(q)
    if dimension == 1:
        return 1/(2*m)*p*p + 1/2*m*(a_0 + a_1*q + a_2*q**2 + a_3*q**3 + a_4*q**4)
    else:
        return 1/(2*m)*np.dot(p, p) + 1/2*m*(a_0 + a_1*q + a_2*q**2 + a_3*q**3 + a_4*q**4)
    
def asymmetric_Oscillator_H_eq(Parameters, p, q):
    """
        Hamilton's equations of motion for a simple harmonic oscillator.
        Args:
            q: Position coordinate
            p: Momentum coordinate
        Returns:
            dq_dt: Derivative of position
            dp_dt: Derivative of momentum
        """
    
    # modelParameters: 
    dimension = Parameters.modelParameters[0]
    m = Parameters.parseParameter('m')
    a_1 = Parameters.parseParameter('a_1')
    a_2 = Parameters.parseParameter('a_2')
    a_3 = Parameters.parseParameter('a_3')
    a_4 = Parameters.parseParameter('a_4')
    
    # Compute the derivatives
    dq_dt = p/m
    dp_dt = -1/2*m*(a_1 + 2*a_2*q+ 3*a_3*q**2 + 4*a_4*q**3)
    return dq_dt, dp_dt