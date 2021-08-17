import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv("iters.csv")
ovi1=df["OVI_1"]
ovi1opt=df["OVI_1_opt"]
ovi100=df["OVI_100"]
ovi100opt=df["OVI_100_opt"]

fig = plt.figure()
ax = fig.add_axes([0,0,1,1])

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_title("OVI opt comparison")
ax.set_xlabel("Without opt")
ax.set_ylabel("With opt")

ax.scatter(ovi1,ovi1opt,color='r')
ax.scatter(ovi100,ovi100opt,color='b')
plt.show()

