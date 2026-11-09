import pandas as pd, matplotlib.pyplot as plt
t = pd.read_csv("results_p2.csv"); t["cells"] = t.m * t.n
mem = pd.read_csv("mem_p2.csv");   mem["cells"] = mem.m * mem.n

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for p, g in t.groupby("p"):
    ax[0].loglog(g.cells, g.avg_time_us, "o-", label=f"p={p}")
ax[0].set(xlabel="m·n (cells)", ylabel="time (µs)", title="Running time"); ax[0].legend()

g = mem[mem.p == 0.1]
ax[1].loglog(g.cells, g.theory_bytes / 1024, "o-", label="theoretical")
ax[1].loglog(g.cells, g.heap_peak_bytes / 1024, "x--", label="measured heap peak")
ax[1].loglog(g.cells, (g.rss_after_kb - g.rss_before_kb).clip(lower=1), "s:", label="measured ΔRSS (4KB pages)")
ax[1].set(xlabel="m·n (cells)", ylabel="memory (KB)", title="Memory usage"); ax[1].legend()

plt.tight_layout(); plt.savefig("p2_plots.png", dpi=150)
