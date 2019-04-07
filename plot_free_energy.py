import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import itertools

x, y, z = np.genfromtxt(r'fel_teste.txt', unpack=True)

# Generate fake data
#x = np.random.normal(size=1000)
#y = x * 3 + np.random.normal(size=1000)
cm = plt.cm.get_cmap('hot')
# Calculate the point density
xy = np.vstack([x,y])

idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

xmin = 5
xmax = 45
ymin = 5
ymax = 50

fig, ax = plt.subplots()
p = ax.scatter(x, y, c=z, s=100, cmap=cm, edgecolors='none')
plt.colorbar(p)
plt.xlabel("X")
plt.ylabel("Y")
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.show()
