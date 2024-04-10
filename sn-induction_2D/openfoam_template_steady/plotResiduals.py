import pandas as pd
import matplotlib.pyplot as plt

data = {}

values = ["p_rgh_0", "Ux_0", "Uy_0", "Uz_0", "T_0"]
f, ax = plt.subplots(figsize=(15, 5))

for key in values:
    data[key] = pd.read_csv(
        "logs/" + key, delimiter="\t", header=None, names=["iter", "residual"]
    )
    plot = ax.plot(data[key]["iter"], data[key]["residual"])

ax.set_yscale("log")
ax = plt.gca()
ax.legend(data.keys(), loc="upper right")
ax.set_xlabel("Iterations")
ax.set_ylabel("Residuals")
plt.savefig("residuals.png", dpi=300)
