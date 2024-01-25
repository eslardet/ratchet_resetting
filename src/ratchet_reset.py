import numpy as np



class Particle:
    def __init__(self, x0, dt, D, L, a, h, r):
        self.x = x0
        self.phase = 'd' ## Can be 'd' for diffusion or 'r' for ratchet/return
        self.t = 0
        self.dt = dt
        self.D = D
        self.L = L
        self.a = a
        self.h = h
        self.r = r
    
    def pbc_wrap(self):
        self.x = self.x % self.L

    def get_force(self):
        if self.x < self.a:
            return -self.h/self.a
        elif self.x > self.a:
            return self.h/(self.L - self.a)
        else:
            return 0

    def move(self):
        if self.phase == 'd':
            prob_r = np.random.uniform(0,1)
            if prob_r < self.r*self.dt:
                self.phase = 'r'
                self.x += self.get_force()*self.dt + np.sqrt(2*self.D*self.dt)*np.random.normal()
            else:
                self.x += np.sqrt(2*self.D*self.dt)*np.random.normal()
            self.pbc_wrap()

        elif self.phase == 'r':
            self.x += self.get_force()*self.dt + np.sqrt(2*self.D*self.dt)*np.random.normal()
            # print(self.get_force()*self.dt, np.sqrt(2*self.D*self.dt)*np.random.normal())
            if self.x <= 0 or self.x >= self.L:
                self.phase = 'd'
                self.pbc_wrap()
        else:
            raise ValueError('Particle phase must be "d" or "r"')
        
        self.t += self.dt

