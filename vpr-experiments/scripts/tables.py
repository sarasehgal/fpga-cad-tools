import json,sys
s=json.load(open('/home/claude/e1/results/summary.json'))
C=["or1200","clma","pdc","riscv_core_lut6"]
f=lambda m,sd: f"{m:.1f} ± {sd:.1f}"
which=sys.argv[1]
if which=="A2":
    print("| Circuit | Seed 1–5 final bb_cost | Mean ± SD | CV (SD/mean) | Range (max−min) |\n| --- | --- | --- | --- | --- |")
    for c in C:
        v=s['A2'][c]; m,sd,n=v['bb']; r=v['raw']
        print(f"| {c} | {', '.join(f'{x:.1f}' for x in r)} | {f(m,sd)} | {100*sd/m:.2f}% | {max(r)-min(r):.1f} ({100*(max(r)-min(r))/m:.1f}%) |")
if which=="A3":
    T=["1e-2","1e-3","1e-4","1e-5","1e-6"]
    print("| Circuit | default (auto) | " + " | ".join(f"init_t = {t}" for t in T) + " |\n| --- | --- | " + " | ".join("---" for _ in T)+" |")
    for c in C:
        d=s['A2'][c]['bb']
        print(f"| {c} | {f(d[0],d[1])} | " + " | ".join(f(*s['A3'][c][t]['bb'][:2]) for t in T)+" |")
if which=="A4":
    print("| Circuit | Range window ON (A.2) | Range window OFF (A.4) | Change in mean | SD ratio (off/on) |\n| --- | --- | --- | --- | --- |")
    for c in C:
        a=s['A2'][c]['bb']; b=s['A4'][c]['bb']
        print(f"| {c} | {f(a[0],a[1])} | {f(b[0],b[1])} | {100*(b[0]/a[0]-1):+.1f}% | {b[1]/a[1]:.1f}× |")
