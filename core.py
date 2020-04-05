import numpy as np
from scipy.optimize import brentq
from math import sqrt, exp
from scipy.optimize import curve_fit, basinhopping
import matplotlib.pyplot as plt

class CubicEq:
    def __init__(self, T, R, a, b):
        
        self.T = T
        self.R = R
        self.a = a
        self.b = b
        
    def set_p(self, p):
        self.p = p
        
    def obj_funct(self, v):
        p = self.p
        T = self.T 
        R = self.R 
        a = self.a
        b = self.b 
        return p*sqrt(T)*v**3-R*T**(3.0/2)*v**2+(a-R*T**(3/2)*b-p*b**2*sqrt(T))*v-a*b        


def select_valid_root(r):
    if r[0].imag == 0 and r[0].real > 0.0:
        return r[0].real
    if r[1].imag == 0 and r[1].real > 0.0:
        return r[1].real    
    if r[2].imag == 0 and r[2].real > 0.0:
        return r[2].real
    print("Problem with the roots of cubic volumen equation")
    return None     

class YCurve:
    def __init__(self, TC = -240.17+273.15, PC=12.77*101325, T=75.0):
        self.R = 8.314472
        self.a = (0.42748*(self.R**2)*(TC**(2.5)))/PC
        self.b = 0.08664*self.R*TC/PC
        self.T = T        

    def cubicEq(self, v):
        p = self.p
        T = self.T 
        R = self.R 
        a = self.a
        b = self.b 
        
        return p*np.sqrt(T)*v**3-R*T**(3.0/2)*v**2+(a-R*T**(3/2)*b-p*b**2*np.sqrt(T))*v-a*b


    def findRoot(self):
        T = self.T
        R = self.R
        a, b = self.a, self.b 

        if type(self.p) is np.ndarray:
            v = np.empty(shape=self.p.shape)
            for i in range(v.size):
                p = self.p[i]
                c3 = p*sqrt(T)                    
                c2 = -R*T**(3.0/2)
                c1 = (a-R*T**(3/2)*b-p*b**2*sqrt(T))
                c0 = -a*b
                v[i] = select_valid_root(np.roots([c3, c2, c1, c0])) 
            return v        
        else:
            p = self.p
            c3 = p*sqrt(T)                    
            c2 = -R*T**(3.0/2)
            c1 = (a-R*T**(3/2)*b-p*b**2*sqrt(T))
            c0 = -a*b
            return select_valid_root(np.roots([c3, c2, c1, c0])) 

    def curve(self, p, n0, q, n, vad):
        self.p = p
        v = self.findRoot()        
        return n0*(1.0-np.exp(-q*p**n))-vad/v
    
class ObjetiveFunct:
    def __init__(self, xn, yn, curve):
        self.xn = xn
        self.yn = yn
        self.curve = curve
        
    def obj_funct(self, params):
        n0, q, n, vad = params
        yp = self.curve.curve(self.xn, n0=n0, q=q, n=n, vad=vad)
        return ((self.yn-yp)**2).sum()

def zhou_funct(p, n0, q, n, c1, c2):
    return n0*(1.0-np.exp(-q*p**n)) - c1*1.0e-4*c2*1.0e-6*p**2

class ZhouObjetiveFunct:
    def __init__(self, xn, yn):
        self.xn = xn
        self.yn = yn
        
    def obj_funct(self, params):
        n0, q, n, c1, c2 = params
        yp = zhou_funct(self.xn, n0, q, n, c1, c2)
        return ((self.yn-yp)**2).sum()

def bondad(yn, fn):
    ss_tot = ((yn - yn.mean())**2).sum()
    ss_res = ((yn-fn)**2).sum()
    return 1 - ss_res/ss_tot

def rel_error(yn, yp):
    return 100.0* ( np.abs(yn-yp)/yn ).mean()

if __name__ == "__main__":
    pass 
