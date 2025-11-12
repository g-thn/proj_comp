import ply
import numpy as np
import matplotlib.pyplot as plt

# Brute force example

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

angle_configurations = [0.0,30.,-30, 45.0, -45.0, 90.]
energy_max = 0.
best_config = None
progress = 0
nb_valid_config = 0
for a1 in angle_configurations:
    for a2 in angle_configurations:
        for a3 in angle_configurations:
            for a4 in angle_configurations:
                angles_deg = [a1, a2, a3, a4, a4, a3, a2, a1]
                ply_list = []
                for i in range(nb_layers):
                    theta = angles_deg[i]*np.pi/180.
                    ply_tmp = ply.Ply(e1, e2, g12, nu12, z_tops[i], z_bots[i], theta, x=xc, y=yc, s=sc)
                    ply_list.append(ply_tmp)
                laminate = ply.Laminate(ply_list)
                laminate.abd_global()
                laminate.set_load(np.array([[1.e4],
                                            [-5.e4],
                                            [0.]]),
                                  np.array([[1.e4],
                                            [1e4],
                                            [0.]]))
                laminate.comp_deformation()
                laminate.update_tsai_hill()
                
                if laminate.tscrit < 1.:
                    nb_valid_config += 1
                    energy_max = max(energy_max, laminate.comp_elastic_energy())
                    if energy_max == laminate.elastic_energy:
                        best_config = [a1, a2, a3, a4]
                progress += 1
                if progress % 10 == 0:
                    print(f"Progress: {100*progress/4**len(angle_configurations)} % evaluated.", end='\r')
                del laminate

print("Number of valid configurations:", nb_valid_config)
print("Best angle configuration (deg):", best_config)
print("Maximum elastic energy over all configurations:", energy_max)