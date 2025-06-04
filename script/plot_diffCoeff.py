import matplotlib.pyplot as plt
import pandas as pd

#plot D vs persistence time (fixed speed)
df_tau = pd.read_csv("../output/diffusion_vs_tau.txt", delim_whitespace=True)

plt.figure(figsize=(6, 5))
plt.plot(df_tau["tau"], df_tau["D"], marker='o', color='blue')
plt.xlabel("Persistence time (τ)")
plt.ylabel("Diffusion Coefficient (D)")
plt.title("D vs τ (fixed speed)")
plt.grid(True)
plt.tight_layout()
plt.savefig("../output/D_vs_tau_fixed_speed.png")
plt.show()

#plot D vs speed (fixed persistence)
df_speed = pd.read_csv("../output/diffusion_vs_speed.txt", delim_whitespace=True)

plt.figure(figsize=(6, 5))
plt.plot(df_speed["speed"], df_speed["D"], marker='o', color='green')
plt.xlabel("Speed (v)")
plt.ylabel("Diffusion Coefficient (D)")
plt.title("D vs v (fixed persistence time)")
plt.grid(True)
plt.tight_layout()
plt.savefig("../output/D_vs_speed_fixed_tau.png")
plt.show()
