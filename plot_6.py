import matplotlib.pyplot as plt
import numpy as np

with open('data/jump.dat','r') as f: 
    lines = [[float(v) for v in l.strip().split("\t")] for l in f]

t = [v[0] for v in lines]
x = [v[1] for v in lines]
c = [v[2] for v in lines]


fig,ax = plt.subplots()
ax.plot(t,x,label="height",color="orange")
ax.set_ylabel("height",color="orange",fontsize=12)
ax2=ax.twinx()
ax2.plot(t,c,label="speed",color="blue")
ax2.set_ylabel("velocity",color="blue",fontsize=12)
plt.show()