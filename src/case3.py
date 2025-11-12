import ply
import numpy as np
import matplotlib.pyplot as plt

# 8 plies example

e1 = 25e9
e2 = 5e9
g12 = 10e9
nu12 = 0.4
xc = 666e6
yc = 666e6
sc = 94e6
thickness = 20e-3
nb_layers = 8
z_all = np.linspace(thickness/2, -thickness/2, nb_layers+1)
z_tops = z_all[:-1]
z_bots = z_all[1:]
angles_deg = [0.0, 45.0, -45.0, 90.0, 90.0, -45.0, 45.0, 0.0]
ply_list = []
for i in range(nb_layers):
    theta = angles_deg[i]*np.pi/180.
    ply_tmp = ply.Ply(e1, e2, g12, nu12, z_tops[i], z_bots[i], theta, x=xc, y=yc, s=sc)
    ply_list.append(ply_tmp)
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
laminate.update_tsai_hill()
print("Tsai-Hill criterion overall:", laminate.tscrit)
print("Elastic energy of the laminate:", laminate.comp_elastic_energy())