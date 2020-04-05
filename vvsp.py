from core import YCurve

curve = YCurve(TC=-240.17+273.15, PC=12.77*101325, T=273.15)
eq = curve.getCubicEq()

p = 101000.325
print("{:>16s} -> {:>16s}".format("Pressure", "Volume"))
print("{:16.8f} -> {:16.8f}".format(p, eq.get_root(p)))

