import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("../output/clearance_times.csv")
df.columns = df.columns.str.strip()

df_grouped = df.groupby("number_best", as_index=False)["clearance_time"].mean()

plt.plot(df_grouped["number_best"], df_grouped["clearance_time"], marker='o')

plt.xlabel("Number of Best Cells")
plt.ylabel("Clearance Time")
plt.title("Clearance Time Depending on the Number of Best Cells.")
plt.grid(True)

for x, y in zip(df_grouped["number_best"], df_grouped["clearance_time"]):
    plt.text(x, y-5, f"{y:.1f}", fontsize=8, ha='left', va='bottom')

plt.savefig("../output/clearance_plot.png", dpi=300)

plt.show()