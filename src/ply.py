import numpy as np
import matplotlib.pyplot as plt

class Ply:
    def __init__(self,e1, e2, g12, nu12,z_top,z_bot,theta, x = 0., y = 0., s = 0.):
        self.e1 = e1
        self.e2 = e2
        self.g12 = g12
        self.nu12 = nu12
        self.nu21 = self.nu12 * self.e2 / self.e1
        self.z_top = z_top
        self.z_bot = z_bot
        self.thickness = z_top - z_bot
        self.theta = theta
        self.x = x
        self.y = y
        self.s = s
        self.compute_invmatrot()
        self.matrot = np.linalg.inv(self.invmatrot)
        self.q_loc()
        self.q_glob()
        self.a_loc()
        self.b_loc()
        self.d_loc()
        self.abd_loc()
        self.eps0 = np.array([[0],
                              [0],
                              [0]])
        self.phi = np.array([[0],
                             [0],
                             [0]])


    def compute_invmatrot(self):
        c = np.cos(self.theta)
        s = np.sin(self.theta)
        self.invmatrot = np.array([[c**2, s**2,      -2*s*c],
                                   [s**2, c**2,       2*s*c],
                                   [ c*s, -c*s, c**2 - s**2]])
        self.matrot = np.linalg.inv(self.invmatrot)
        return self.invmatrot

    def q_loc(self):
        self.q_local = np.array([[self.e1/(1 - self.nu12*self.nu21)        , self.nu21*self.e1/(1 - self.nu12*self.nu21), 0       ],
                                 [self.nu12*self.e2/(1 - self.nu12*self.nu21), self.e2/(1 - self.nu12*self.nu21)          , 0       ],
                                 [0                                        , 0                                          , self.g12]])
        return self.q_local
     
    def q_glob(self):
        self.q_global = np.dot(np.dot(self.invmatrot, self.q_local), self.invmatrot.T)
        return self.q_global
    
    def a_loc(self):
        self.a = self.q_global*(self.z_top - self.z_bot)
        return self.a
    
    def b_loc(self):
        self.b = 0.5*self.q_global*(self.z_top**2 - self.z_bot**2)
        return self.b
    
    def d_loc(self):
        self.d = self.q_global*(self.z_top**3 - self.z_bot**3)/3
        return self.d
    
    def abd_loc(self):
        a = self.a_loc()
        b = self.b_loc()
        d = self.d_loc()
        abd = np.block([[a, b],
                        [b, d]])
        self.abd = abd
        return abd
    
    def set_eps_phi(self, eps, phi):
        self.eps = eps
        self.phi = phi
    
    def strain(self,z):
        return self.eps + z*self.phi

    def stress(self, z):
        stress = np.dot(self.q_global, self.strain(z))
        return stress
    
    def set_param_Ts(self, x, y, s):
        self.x = x
        self.y = y
        self.s = s

    def tsai_hill_crit(self):
        crit_values = []
        for z in [self.z_top, self.z_bot]:
            sig_glob = self.stress(z)
            sig_loc = np.dot(self.matrot, sig_glob)
            sig1 = sig_loc[0,0]
            sig2 = sig_loc[1,0]
            tau12 = sig_loc[2,0]
            crit = (sig1/self.x)**2 - (sig1*sig2)/(self.x**2) + (sig2/self.y)**2 + (tau12/self.s)**2
            crit_values.append(crit)
        self.tscrit =  np.max(crit_values)
        return self.tscrit

class Laminate:
    def __init__(self, ply_list):
        self.ply_list = ply_list
        self.abd_global()
        self.abd_inv = np.linalg.inv(self.abd)
        self.load = np.zeros((6,1))
        self.deformation = np.zeros((6,1))
        self.eps0 = np.zeros((3,1))
        self.phi = np.zeros((3,1))
        self.force = np.zeros((3,1))
        self.moment = np.zeros((3,1))
        self.comp_thickness()
        self.z_top = .5*self.thickness
        self.z_bot = -.5*self.thickness
        self.tscrit = 0.

    
    def comp_thickness(self):
        self.thickness = 0.
        for ply in self.ply_list:
            self.thickness += ply.thickness

    def set_load(self, force, moment):
        self.force = force
        self.moment = moment
        self.load = np.vstack((force, moment))
    
    def set_deformation(self, eps_top, eps_bot):
        self.eps0 = 0.5*(eps_top + eps_bot)
        self.phi = (eps_top - eps_bot)/self.thickness
        self.deformation = np.vstack((self.eps0, self.phi))
        self.update_plies()
        
    def update_plies(self):
        for ply in self.ply_list:
            ply.set_eps_phi(self.eps0, self.phi)

    def abd_global(self):
        a_tot = np.zeros((3,3))
        b_tot = np.zeros((3,3))
        d_tot = np.zeros((3,3))
        for ply in self.ply_list:
            a_tot += ply.a_loc()
            b_tot += ply.b_loc()
            d_tot += ply.d_loc()
        abd_tot = np.block([[a_tot, b_tot],
                            [b_tot, d_tot]])
        self.abd = abd_tot
        return abd_tot
    
    def comp_deformation(self):
        self.deformation = np.dot(self.abd_inv, self.load)
        self.eps0 = self.deformation[0:3,:]
        self.phi = self.deformation[3:6,:]
        self.update_plies()
        return self.deformation
    
    def comp_load(self):
        self.load = np.dot(self.abd, self.deformation)
        self.force = self.load[0:3,:]
        self.moment = self.load[3:6,:]
        return self.load

    def update_tsai_hill(self):
        tscrit_list = []
        for ply in self.ply_list:
            ply.set_eps_phi(self.eps0, self.phi)
            tscrit = ply.tsai_hill_crit()
            tscrit_list.append(tscrit)
        self.tscrit = np.max(tscrit_list)

    def plot_laminate(self):
        pass
    
    def plot_stress_distribution(self, num_points=100):
        pass
    
if __name__ == "__main__":
    
    # Param ply
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
    theta3 = 0.*np.pi/180.
    z1 = 4-3
    z2 = 3e-3
    z3 = -3e-3
    z4 = -4e-3
    ply1 = Ply(e1, e2, g12, nu12, z1, z2, theta1, x=xc, y=yc, s=sc)
    ply2 = Ply(ef, ef, gf, nuf, z2, z3, theta2, x=xf, y=yf, s=sf)
    ply3 = Ply(e1, e2, g12, nu12, z3, z4, theta3, x=xc, y=yc, s=sc)
    ply_list = [ply1, ply2, ply3]
    laminate = Laminate(ply_list)
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
    laminate.plot_laminate()
    laminate.plot_stress_distribution()