import numpy as np
import matplotlib.pyplot as plt
import time
from particle import Particle

seed = 1
np.random.seed(seed)

# Simulation parameters
L = 1
a = 0.2
h = 1
D = 0.1
dt = 0.001
r = 0.2
x0 = 0

samples = 10**3
total_t = 10

# num_steps = 20000
num_steps = int(total_t/dt)

final_loc = []
t0 = time.time()
for s in range(samples):
    # Create particle
    particle = Particle(x0, dt, D, L, a, h, r)
    for i in range(num_steps):
        particle.move()
    final_loc.append(particle.x)
print('Time taken: ', time.time() - t0)
fig, ax = plt.subplots()
ax.hist(final_loc, bins=100)
plt.show()