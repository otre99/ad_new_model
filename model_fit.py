import argparse
from core import YCurve
import numpy as np 
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt 

plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)
plt.rc('font', size=18, weight="bold")
plt.rc('figure', autolayout=True, figsize=(9,6))
plt.rc('axes', titlesize=16, labelsize=18, labelweight="bold", linewidth=2.0)
plt.rc('lines', linewidth=2, markersize=6)
plt.rc('legend',fontsize=16)
plt.rc('mathtext', fontset='stix')
plt.rc('font', family='STIXGeneral')
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

if __name__ == "__main__":
    from math import isinf

    def inside_range(a,b, x, vname):
        if not a<x<b:
            print("Warning: Initial value for parameter `{}` " 
                  "must be in ({},{}), however you set {}={}".format(vname, a,b,vname,x))
            return False
        return True 

    parser = argparse.ArgumentParser()
    parser.add_argument("data_file", help="Experimental data")
    parser.add_argument("--xsf", help="Scale factor for x experimental data", type=float, default=1.0)
    parser.add_argument("--ysf", help="Scale factor for y experimental data", type=float, default=1.0)
    
    parser.add_argument("--pc", help="Critical pressure [Pa]", type=float, default=1293920.25)
    parser.add_argument("--tc", help="Critical temperature [K]", type=float, default=32.979999)
    parser.add_argument("--temp", help="Temperature [K]", type=float, default=75.0)    


    parser.add_argument("--n0",  help="n0 (Starting Value)", type=float, default=3.728490)
    parser.add_argument("--q",   help="q (Starting Value)", type=float, default=0.010016)
    parser.add_argument("--n",   help="n (Starting Value)", type=float, default=0.414666)    
    parser.add_argument("--vad", help="vad (Starting Value)", type=float, default=0.000288)
    parser.add_argument('--fit', action='store_true')

    
    FLAGS, _ = parser.parse_known_args()    
    data = np.loadtxt(FLAGS.data_file)    
    model = YCurve(TC=FLAGS.tc, PC=FLAGS.pc, T=FLAGS.temp)

    # Y model 
    n0  = FLAGS.n0  
    q   = FLAGS.q   
    n   = FLAGS.n   
    vad = FLAGS.vad 

    xn, yn = data[:,0], data[:,1]
    xn*=FLAGS.xsf
    yn*=FLAGS.ysf
        
    if FLAGS.fit:
        n0Min = -np.inf
        n0Max = +np.inf 
        f1 = inside_range(n0Min, n0Max,n0,"n0")

        qMin = -np.inf
        qMax = +np.inf 
        f2 = inside_range(qMin, qMax,n0,"q")
        
        nMin = -np.inf 
        nMax = +np.inf 
        f3 = inside_range(nMin, nMax,n,"n")
        
        vadMin = -np.inf
        vadMax = +np.inf
        f4 = inside_range(vadMin, vadMax,vad,"vad")
        
        # ensure a feasible staring point         
        if f1&f2&f3&f4:        
            bounds = ([n0Min, qMin, nMin, vadMin], [n0Max, qMax, nMax, vadMax])
            yp0, results = curve_fit(model.curve, xn, yn, 
                                 p0=[n0, q, n, vad], bounds=bounds)
            n0, q, n, vad = yp0
        else:
            print("Error. Not fitted was done!!")
    
    xx = np.linspace(xn.min(), xn.max(), 12000)    
    yy = model.curve(xx, n0=n0, q=q, n=n, vad=vad)
        
    plt.scatter(xn, yn, label="data", color="black")
    plt.plot(xx, yy, label="model", color="black")
    plt.xlabel("Pa")
    plt.ylabel(r"mol . mol$^{-1}$")
    plt.grid()
    plt.legend()
    tt = r"n0={:.6f} q={:.6f} n={:.6f} Vad={:.6f}".format(n0, q, n, vad)
    plt.title(tt)
    plt.show()
    
 
