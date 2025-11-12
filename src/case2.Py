import ply
import numpy as np
import matplotlib.pyplot as plt

# S3 Example with foam core

e1 = 77e9
e2 = 77e9
g12 = 10e9
nu12 = 0.
xc = 666e6
yc = 666e6
sc = 94e6
# Param foam
ef = 27e6
nuf = 0.1
gf = ef/(2*(1+nuf))
xf = .4e6
yf = .4e6
sf = .15e6
theta1 = 0.*np.pi/180.
theta2 = 0.*np.pi/180.
theta3 = -0.*np.pi/180.
z1 = 4e-3
z2 = 3e-3
z3 = -3e-3
z4 = -4e-3
ply1 = ply.Ply(e1, e2, g12, nu12, z1, z2, theta1, x=xc, y=yc, s=sc)
ply2 = ply.Ply(ef, ef, gf, nuf, z2, z3, theta2, x=xf, y=yf, s=sf)
ply3 = ply.Ply(e1, e2, g12, nu12, z3, z4, theta3, x=xc, y=yc, s=sc)
ply_list = [ply1, ply2, ply3]
laminate = ply.Laminate(ply_list)
laminate.abd_global()
eps_top = 1e-3*np.array([[1.75],
                            [-6.37],
                            [0.0]])
eps_bot = 1e-3*np.array([[-2.46],
                            [6.6],
                            [0.0]])
laminate.set_deformation(eps_top, eps_bot)
print("Laminate eps0:", laminate.eps0.flatten())
print("Laminate phi:", laminate.phi.flatten())
laminate.update_tsai_hill()
print("Tsai-Hill criterion:", laminate.ply_list[0].tscrit, laminate.ply_list[1].tscrit, laminate.ply_list[2].tscrit)
print("Max Tsai-Hill criterion;", laminate.tscrit)
print("Ply 1 stress", laminate.ply_list[0].stress(np.linspace(ply1.z_bot,ply1.z_top,10)))
laminate.plot_laminate()
laminate.plot_stress_distribution()