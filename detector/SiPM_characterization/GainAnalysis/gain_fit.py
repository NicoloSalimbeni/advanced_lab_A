import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.odr import *

filepath="txt_files/C_SX"
df=pd.read_csv(filepath+".txt",sep="\t")

x_data=df["V[dV]"]
y_data=df["Gain[mV]"]

err_x=df["SigmaV[dV]"]
err_y=df["SigmaGain[mv]"]

def linear(p,x):
    m,q=p
    return m*x+q


# Create a model for fitting.
linear_model = Model(linear)

# Create a RealData object using our initiated data from above.
data = RealData(x_data, y_data, sx=err_x, sy=err_y)

# Set up ODR with the model and data.
odr = ODR(data, linear_model, beta0=[0., -50.])

# Run the regression.
out = odr.run()

# Use the in-built pprint method to give us results.
out.pprint()

fig,ax=plt.subplots(1,1,figsize=(13,12))
x_fit=np.linspace(280,x_data[len(x_data)-1],1000)
V_b=-out.beta[1]/out.beta[0]
ax.errorbar(x=x_data,y=y_data,xerr=err_x,yerr=err_y,linestyle="none",capsize=8)
ax.plot(x_fit,linear(out.beta,x_fit),linestyle="dashed")
ax.axvline(V_b, color='red', linestyle='dashed', linewidth=1,label= "$V_b$: {} [dV]".format(round(V_b,3)))
ax.set_title("Gain vs $V_{bias}$",fontsize=22)
ax.set_xlabel("$V_{bias}$ [dV]",fontsize=16)
ax.set_ylabel("Gain [mV]",fontsize=16)
ax.grid()
ax.legend(loc="upper left",fontsize=14)
plt.show()
fig.savefig("filepath.pdf")

