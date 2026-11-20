#!/usr/bin/env python3
"""Parse VPR logs for ECE1387 Exercise 1 and produce CSV tables + plots."""
import re, os, glob, json, statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/claude/e1/runs"
OUT = "/home/claude/e1/results"
os.makedirs(OUT, exist_ok=True)
CIRCUITS = ["or1200", "clma", "pdc", "riscv_core_lut6"]
ROW = re.compile(r"^\s*(\d+)\s+[\d.]+\s+([\d.]+e[+-]\d+)\s+([\d.]+)\s+([\d.]+)\s")


def read(path):
    try:
        with open(path) as f:
            return f.read()
    except OSError:
        return None


def place_info(log):
    m = re.findall(r"Placement cost: [\d.]+, bb_cost: ([\d.]+)", log)
    w = re.findall(r"BB estimate of min-dist \(placement\) wire length: (\d+)", log)
    return (float(m[-1]) if m else None, int(w[-1]) if w else None)


def anneal_curve(log):
    sec = log.split("Tnum", 1)[1].split("Placement Quench", 1)[0]
    pts = []
    for line in sec.splitlines():
        # first data row is glued onto the "(sec) (ns)" header line
        line = re.sub(r"^\s*\(sec\).*?\(ns\)\s*", "", line)
        m = ROW.match(line)
        if m:
            pts.append((int(m.group(1)), float(m.group(2)), float(m.group(4))))
    return pts


def ms(xs):
    xs = [x for x in xs if x is not None]
    if not xs:
        return (None, None, 0)
    return (st.mean(xs), st.stdev(xs) if len(xs) > 1 else 0.0, len(xs))


def seeds(d):
    return sorted(glob.glob(os.path.join(d, "s*")), key=lambda p: int(p.rsplit("s", 1)[1]))


summary = {}

# ---------- Part A ----------
# A1: cost vs temperature (seed 1 default auto schedule)
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
for ax, c in zip(axes.flat, CIRCUITS):
    log = read(f"{ROOT}/A2/{c}/s1/vpr_stdout.log")
    if not log:
        ax.set_title(f"{c} (no data)")
        continue
    pts = [p for p in anneal_curve(log) if p[1] > 0]
    T = [p[1] for p in pts]; C = [p[2] for p in pts]
    ax.plot(T, C, marker=".", lw=1)
    ax.set_xscale("log"); ax.invert_xaxis()
    ax.set_xlabel("Temperature T (log scale, decreasing →)")
    ax.set_ylabel("Avg bounding-box cost (HPWL)")
    ax.set_title(c)
    ax.grid(alpha=.3)
    with open(f"{OUT}/A1_{c}_cost_vs_T.csv", "w") as f:
        f.write("tnum,T,av_bb_cost\n")
        for p in pts:
            f.write(f"{p[0]},{p[1]},{p[2]}\n")
fig.suptitle("Exercise A.1 — HPWL cost vs. temperature (default auto schedule, seed 1)")
fig.tight_layout(); fig.savefig(f"{OUT}/A1_cost_vs_T.png", dpi=150); plt.close(fig)


def a_table(sub):
    res = {}
    for c in CIRCUITS:
        vals = [place_info(read(f"{d}/vpr_stdout.log") or "") for d in seeds(f"{ROOT}/{sub}/{c}")]
        res[c] = {"bb": ms([v[0] for v in vals]), "wl": ms([v[1] for v in vals]),
                  "raw": [v[0] for v in vals]}
    return res

summary["A2"] = a_table("A2")
summary["A4"] = a_table("A4")
summary["A3"] = {}
TEMPS = ["1e-2", "1e-3", "1e-4", "1e-5", "1e-6"]
for c in CIRCUITS:
    summary["A3"][c] = {}
    for T in TEMPS:
        vals = [place_info(read(f"{d}/vpr_stdout.log") or "") for d in seeds(f"{ROOT}/A3/{c}/T{T}")]
        summary["A3"][c][T] = {"bb": ms([v[0] for v in vals]), "raw": [v[0] for v in vals]}

# A3 plot
fig, axes = plt.subplots(1, 4, figsize=(15, 3.8))
for ax, c in zip(axes, CIRCUITS):
    xs, m, s = [], [], []
    for T in TEMPS:
        r = summary["A3"][c][T]["bb"]
        if r[0] is not None:
            xs.append(float(T)); m.append(r[0]); s.append(r[1])
    if xs:
        ax.errorbar(xs, m, yerr=s, marker="o", capsize=3)
        d = summary["A2"][c]["bb"]
        if d[0]:
            ax.axhline(d[0], ls="--", color="gray", label="default auto schedule")
            ax.legend(fontsize=7)
    ax.set_xscale("log"); ax.set_title(c); ax.set_xlabel("init_t"); ax.grid(alpha=.3)
axes[0].set_ylabel("final bb_cost (mean ± sd, 5 seeds)")
fig.suptitle("Exercise A.3 — final HPWL vs. starting temperature (exit_t=1e-8, alpha_t=0.9)")
fig.tight_layout(); fig.savefig(f"{OUT}/A3_start_temp.png", dpi=150); plt.close(fig)

# ---------- Part B ----------
def route_info(log):
    if not log:
        return None, None, False, None
    cpd = re.findall(r"Final critical path delay \(least slack\): ([\d.]+) ns", log)
    wl = re.findall(r"Total wirelength: (\d+)", log)
    ok = "Circuit successfully routed" in log
    it = re.findall(r"Successfully routed after (\d+) routing iterations", log)
    return (float(cpd[-1]) if cpd and ok else None, int(wl[-1]) if wl and ok else None, ok, int(it[-1]) if it else None)


def b_table(sub):
    res = {}
    for c in CIRCUITS:
        vals = [route_info(read(f"{d}/vpr_stdout.log")) for d in seeds(f"{ROOT}/{sub}/{c}")]
        res[c] = {"cpd": ms([v[0] for v in vals]), "wl": ms([v[1] for v in vals]),
                  "routed": sum(v[2] for v in vals), "n": len(vals), "iters": ms([v[3] for v in vals]), "raw": [v[0] for v in vals]}
    return res

for sub in ["B1", "B3a", "B3b"]:
    summary[sub] = b_table(sub)
summary["B2"] = {}
for c in CIRCUITS:
    summary["B2"][c] = {}
    for d in sorted(glob.glob(f"{ROOT}/B2/{c}/*")):
        name = os.path.basename(d)
        vals = [route_info(read(f"{x}/vpr_stdout.log")) for x in seeds(d)]
        summary["B2"][c][name] = {"cpd": ms([v[0] for v in vals]), "wl": ms([v[1] for v in vals]),
                                  "routed": sum(v[2] for v in vals), "n": len(vals),
                                  "iters": ms([v[3] for v in vals]), "raw": [v[0] for v in vals]}

with open(f"{OUT}/summary.json", "w") as f:
    json.dump(summary, f, indent=1)
print("ok")
