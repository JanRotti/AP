import numpy as np
import matplotlib.pyplot as plt
from random import random

x = np.linspace(0,1,26)
y =  4 * np.sqrt(x) - 3*x - 1 * x**2 + 2 * x**3 - 2 * x**4


plt.plot(x,y,'o')
plt.show()

with open('data/wing.dat','w+') as f:
    for x,y in zip(x,y):
        f.write("{} {}\n".format(x,y))