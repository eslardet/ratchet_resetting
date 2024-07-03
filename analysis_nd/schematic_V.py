import matplotlib.pyplot as plt
import numpy as np
import os


fig, ax = plt.subplots(figsize=(6, 2))


def V(x, a, h, L):
    x = x % L
    if x <= a:
        return h/a * x
    else:
        return h/(L-a) * (L-x)

a = 0.2
h = 0.1
L = 1.0

x_plot = np.linspace(-0.2, 2.1, 1000)
y_plot = np.array([-V(x, a, h, L) for x in x_plot])
# ax.plot(x_plot, y_plot, color='black')
ax.plot(x_plot, y_plot, color='tab:red', linewidth=2)
ax.plot(x_plot, y_plot*0, color='tab:blue', linewidth=2)

y_min = np.min(y_plot) - 0.01
y_max = np.max(y_plot) + 0.01

ax.vlines(0, y_min, y_max, color='black', linestyle='--', alpha=0.2, linewidth=5)
ax.vlines(L, y_min, y_max, color='black', linestyle='--', alpha=0.2, linewidth=5)
ax.vlines(2*L, y_min, y_max, color='black', linestyle='--', alpha=0.2, linewidth=5)

ax.set_xlim(-0.2, 2.1)
ax.set_ylim(y_min, y_max)


# remove x ticks
ax.set_xticks([])
ax.set_yticks([])

# remove x and y axis
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

folder = os.path.abspath('./plots/trajectory/')
if not os.path.exists(folder):
    os.makedirs(folder)
filename = 'potential_color'
plt.savefig(os.path.join(folder, filename + '.svg'), bbox_inches='tight')


# plt.show()