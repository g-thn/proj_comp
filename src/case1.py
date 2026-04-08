import ply
import numpy as np
import matplotlib.pyplot as plt


 # S2 example
e1 = 25e9
e2 = 5e9
g12 = 10e9
nu12 = 0.4
xc = 666e6
yc = 666e6
sc = 94e6
theta1 = 30.*np.pi/180.
theta2 = 0.*np.pi/180.
theta3 = -30.*np.pi/180.
z1 =6e-3
z2 = 3e-3
z3 = -3e-3
z4 = -6e-3
ply1 = ply.Ply(e1, e2, g12, nu12, z1, z2, theta1, x=xc, y=yc, s=sc)
ply2 = ply.Ply(e1, e2, g12, nu12, z2, z3, theta2, x=xc, y=yc, s=sc)
ply3 = ply.Ply(e1, e2, g12, nu12, z3, z4, theta3, x=xc, y=yc, s=sc)
ply_list = [ply1, ply2, ply3]
print("Ply 1 stiffness matrix global:\n", ply1.q_global)
print("Ply 2 stiffness matrix global:\n", ply2.q_global)
print("Ply 2 stiffness matrix global:\n", ply2.q_global)
laminate = ply.Laminate(ply_list)
laminate.abd_global()
print("Laminate ABD matrix:\n", laminate.abd)
print("Laminate a:", laminate.a)
print("Laminate b:", laminate.b)
print("Laminate d:", laminate.d)
laminate.set_load(np.array([[1.e6],
                            [0.],
                            [0.]]),
                   np.array([[0.],
                             [0.],
                             [0.]]))
laminate.comp_deformation()
print("Laminate eps0:", laminate.eps0.flatten())
print("Laminate phi:", laminate.phi.flatten())
laminate.update_tsai_hill()
laminate.plot_laminate()
laminate.plot_stress_distribution()
laminate.plot_tsai_hill()