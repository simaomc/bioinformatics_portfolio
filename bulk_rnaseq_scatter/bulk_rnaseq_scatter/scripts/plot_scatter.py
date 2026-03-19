import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
df = pd.read_csv("normalized_counts.csv")

# Example: average expression
df["mit_mean"] = df[["mit1","mit2","mit3","mit4"]].mean(axis=1)
df["pm_mean"] = df[["pm1","pm2","pm3"]].mean(axis=1)

# Log scale
x = np.log10(df["mit_mean"] + 1)
y = np.log10(df["pm_mean"] + 1)

# Plot
plt.scatter(x, y, alpha=0.5)
plt.xlabel("Mitotic (log10)")
plt.ylabel("Post-meiotic (log10)")
plt.title("Gene Expression Scatter Plot")

plt.savefig("../results/scatter_plot.png")
plt.show()
