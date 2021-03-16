import pandas as pd
import matplotlib.pyplot as plt
import os

output_files = ["canal_explicit.dat","canal_implicit.dat"];

for on in output_files:
    ts_name = on.split("_")[1]
    df = pd.read_csv("data/" + on, names=["y","u"])
    plt.plot(df["u"],df["y"], label = ts_name)


plt.legend()
plt.show()