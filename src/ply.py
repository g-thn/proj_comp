import numpy as np
import matplotlib.pyplot as plt

class Ply_comp:
    def __init__(self,e1, e2, g12, nu12,z_top,z_bot,theta):
        self.e1 = e1
        self.e2 = e2
        self.g12 = g12
        self.nu12 = nu12
        self.nu21 = self.nu12 * self.e2 / self.e1
        self.z_top = z_top
        self.z_bot = z_bot
        self.thickness = z_top - z_bot
        self.theta = theta
        self.invmatrot = self.compute_invmatrot()
        self.q_local = self.q_loc()
        self.q_global = self.q_glob()
        self.a = self.a_loc()
        self.b = self.b_loc()
        self.d = self.d_loc()
        self.abd = self.abd_loc()
        self.eps0 = np.array([[0],
                              [0],
                              [0]])
        self.phi = np.array([[0],
                             [0],
                             [0]])


    def compute_invmatrot(self):
        c = np.cos(self.theta)
        s = np.sin(self.theta)
        return np.array([[c**2, s**2,      -2*s*c],
                         [s**2, c**2,       2*s*c],
                         [ c*s, -c*s, c**2 - s**2]])

    def q_loc(self):
        return np.array([[self.e1/(1 - self.nu12*self.nu21)        , self.nu21*self.e1/(1 - self.nu12*self.nu21), 0       ],
                         [self.nu12*self.e2/(1 - self.nu12*self.nu21), self.e2/(1 - self.nu12*self.nu21)          , 0       ],
                         [0                                        , 0                                          , self.g12]])
    def q_glob(self):
        return np.dot(np.dot(self.invmatrot, self.q_local), self.invmatrot.T)
    
    def a_loc(self):
        return self.q_global*(self.z_top - self.z_bot)
    
    def b_loc(self):
        return 0.5*self.q_global*(self.z_top**2 - self.z_bot**2)
    
    def d_loc(self):
        return self.q_global*(self.z_top**3 - self.z_bot**3)/3
    
    def abd_loc(self):
        a = self.a_loc()
        b = self.b_loc()
        d = self.d_loc()
        abd = np.block([[a, b],
                        [b, d]])
        return abd
    
    def set_eps_phi(self, eps, phi):
        self.eps = eps
        self.phi = phi
    
    def strain(self,z):
        return self.eps + z*self.phi

    def stress(self, z):
        stress = np.dot(self.q_global, self.strain(z))
        return stress
    
    def set_param_TS(self, X, Y, S):
        self.X = X
        self.Y = Y
        self.S = S

    def tsai_hill_crit(self):
        sig_glob_top = self.stress(self.z_top)
        sig_glob_bot = self.stress(self.z_bot)
        sig_loc_top = np.dot(self.invmatrot, sig_glob_top)
        sig_loc_bot = np.dot(self.invmatrot, sig_glob_bot)
        pass
        
        

if __name__ == "__main__":
    e1 = 25e9
    e2 = 5e9
    g12 = 10e9
    nu12 = 0.4
    theta1 = 30.*np.pi/180.
    theta2 = 0.*np.pi/180.
    theta3 = -30.*np.pi/180.
    z1 = 6e-3
    z2 = 3e-3
    z3 = -3e-3
    z4 = -6e-3
    ply1 = Ply_comp(e1, e2, g12, nu12, z1, z2, theta1)
    ply2 = Ply_comp(e1, e2, g12, nu12, z2, z3, theta2)
    ply3 = Ply_comp(e1, e2, g12, nu12, z3, z4, theta3)
    a_tot = ply1.a_loc() + ply2.a_loc() + ply3.a_loc()
    b_tot = ply1.b_loc() + ply2.b_loc() + ply3.b_loc()
    d_tot = ply1.d_loc() + ply2.d_loc() + ply3.d_loc()
    print("A matrix:")
    print(a_tot)
    print("B matrix:")
    print(b_tot)
    print("D matrix:")
    print(d_tot)
    eps0 = np.array([[0.001],
                     [0.002],
                     [0.005]])
    phi = np.array([[0.03],
                    [0.01],
                    [0.004]])
    ply1.set_eps_phi(eps0, phi)
    strain_fun = ply1.strain()
    z = np.linspace(ply1.z_bot, ply1.z_top, )
    print("Strain through the thickness of ply 1:")
    print(strain_fun(z))
    plt.plot(z, strain_fun(z)[0,:], label='epsilon_x')
    plt.show()