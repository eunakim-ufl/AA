# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use("Agg")                     # save to file without needing a display
import pandas as pd
import matplotlib.pyplot as plt

r = pd.read_csv("results_p1.csv"); r["cells"] = r.n * r.m
mem = pd.read_csv("mem_p1.csv");   mem["cells"] = mem.n * mem.m
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))

# (1) running time
for exp in ["S1", "S2size"]:
    g = r[(r.experiment == exp) & (r.dist == "english")].sort_values("cells")
    ax[0].loglog(g.cells, g.time_lin_us, "o-", label="{} two-row".format(exp))
    gf = g[g.time_full_us > 0]
    ax[0].loglog(gf.cells, gf.time_full_us, "s--", label="{} full table".format(exp))
ax[0].set_xlabel("n*m (cells)"); ax[0].set_ylabel("time (us)")
ax[0].set_title("Running time"); ax[0].legend()

# (2) memory  (dRSS only for the full table; two-row arrays are below RSS resolution)
for mode, g in mem.groupby("mode"):
    g = g.sort_values("cells")
    ax[1].loglog(g.cells, g.theory_bytes / 1024.0, "o-", label="{} theory".format(mode))
    ax[1].loglog(g.cells, g.heap_peak_bytes / 1024.0, "x--", label="{} heap peak".format(mode))
    if mode == "full":
        ax[1].loglog(g.cells, (g.rss_peak_after_kb - g.rss_cur_before_kb).clip(lower=1),
                     "s:", label="{} dRSS".format(mode))
ax[1].set_xlabel("n*m (cells)"); ax[1].set_ylabel("memory (KB)")
ax[1].set_title("Memory"); ax[1].legend(fontsize=8)

# (3) Scenario 2 delta sweep
s = r[r.experiment == "S2sweep"]
for d, g in s.groupby("dist"):
    g = g.sort_values("delta")
    ax[2].semilogy(g.delta, g.len, "o-", label="{}: length".format(d))
    ax[2].semilogy(g.delta, g.mismatches.clip(lower=0.5), "x--", label="{}: mismatches".format(d))
ax[2].set_xlabel("delta (units of 0.001%)")
ax[2].set_title("Scenario 2: delta sweep (n=m=2000)"); ax[2].legend()

plt.tight_layout(); plt.savefig("p1_plots.png", dpi=150)
print("saved p1_plots.png")
