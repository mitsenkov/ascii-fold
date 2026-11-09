import json, sys
path, chain = sys.argv[1], sys.argv[2]
ss, res = {}, {}
order = []
for l in open(path):
    if l.startswith("HELIX") and l[19] == chain:
        for r in range(int(l[21:25]), int(l[33:37]) + 1): ss[r] = "H"
    elif l.startswith("SHEET") and l[21] == chain:
        for r in range(int(l[22:26]), int(l[33:37]) + 1): ss[r] = "E"
    elif l.startswith("ATOM") and l[21] == chain and l[16] in " A":
        name = l[12:16].strip(); r = int(l[22:26])
        if name in ("CA", "O"):
            if r not in res: res[r] = {}; order.append(r)
            res[r][name] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    elif l.startswith("ENDMDL"): break
rows = [(r, res[r]["CA"], res[r].get("O")) for r in order if "CA" in res[r]]
n = len(rows)
cx = sum(c[1][0] for c in rows)/n; cy = sum(c[1][1] for c in rows)/n; cz = sum(c[1][2] for c in rows)/n
out = []
for r, ca, o in rows:
    o = o or ca
    out.append([round(ca[0]-cx,2), round(ca[1]-cy,2), round(ca[2]-cz,2), ss.get(r,"C"),
                round(o[0]-ca[0],2), round(o[1]-ca[1],2), round(o[2]-ca[2],2)])
print(json.dumps(out, separators=(",",":")))
