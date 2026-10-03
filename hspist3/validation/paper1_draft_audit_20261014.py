#!/usr/bin/env python3
"""##CHRIS 2026-10-14: audit of every Paper 1 draft number that comes from the canonical A1 v2 table, for the box-truncation
correction (methods sec. 14). For each item it prints the definition, the value recomputed from the OLD table (the dated copy
260919_A1v2_final_cs_vs_eta_pre_boxtrunc_20261014.csv), whether that reproduces the draft's printed text, the value from the
NEW (regenerated) table, and the action. Then the per-cell before/after deviations, the regeneration checks, and the numbers of
the Methods paragraph. Nothing is edited unless --apply is given; --apply edits exactly the EDIT strings in paper1_draft.tex and
inserts the Methods paragraph after the "Effective length" paragraph.

usage: python3 paper1_draft_audit_20261014.py [--apply]
"""
import csv, math, os, statistics as st, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tests_20260913 as T
import plot_speed_of_sound_edmd as sos
REPO = os.path.dirname(os.path.dirname(HERE))
TEX = os.path.join(REPO, "0000_PLAN_OVERALL", "paper1_speedofsound", "writeup", "paper1_draft.tex")
OLD = os.path.join(T.PLOTS, "260919_A1v2_final_cs_vs_eta_pre_boxtrunc_20261014.csv")
NEW = os.path.join(T.PLOTS, "260919_A1v2_final_cs_vs_eta.csv")
CENTER_X_CHECK = "5.5\\times10^{-6}"   # printed by paper1_boxtrunc_20261014.py (methods sec. 14.1): "recorded Center_X offset
                                      # vs -delta/2: max |difference| 5.5e-06 sigma" -- quoted, not recomputed (needs the traces)

def load(p):
    return [dict(r, e=float(r["eta"]), c=float(r["c_s"]), es=float(r["c_s_err_scaled"]), er=float(r["c_s_err"]),
                 sm=float(r["c_s_scatter_mass"]), chi=float(r["chi2_red"]), erec=float(r.get("eta_rec") or r["eta"]),
                 kr=float(r["KR"]) if r.get("KR") else None)
            for r in csv.DictReader(open(p))]

def ad(e):   # KR adiabatic, analytic Z' -- the Roman script's function (roman2002_remapped_20260922.py, `adiabatic`)
    a = np.array([e]); return float(sos.cs_adiabatic_2d_monatomic(sos.Z_kolafa_rottner_2006(a), sos.dZ_kolafa_rottner_2006(a), a, kbt=1, m=1)[0])

def dev(r): return 100 * (r["c"] - ad(r["e"])) / ad(r["e"])

def mean_dev_le04(t):   # roman2002_remapped_20260922.py: rows with eta <= 0.69, then eta <= 0.4
    return float(np.mean([dev(r) for r in t if r["e"] <= 0.69 and r["e"] <= 0.4]))

def max_dev_le04(t): return max(dev(r) for r in t if r["e"] <= 0.4)

def zoom_range(t):      # lowdensity_zoom_20260917.py: XMAX = 0.15, N = 100 points
    d = [dev(r) for r in t if r["e"] <= 0.15]; return min(d), max(d)

def at(t, erec):        # one cell, by its recorded eta
    return min(t, key=lambda r: abs(r["erec"] - erec))

def melt(t):            # paper1_melting_figure_20261002.py, same selection and the same nearest-cell rule
    d = sorted([r for r in t if 0.66 <= r["e"] <= 0.765], key=lambda r: r["e"])
    e = np.array([r["e"] for r in d]); imax = int(np.argmin(np.abs(e - 0.700))); imin = int(np.argmin(np.abs(e - 0.715)))
    a, b = d[imax], d[imin]; dip = a["c"] - b["c"]
    return dict(emax=a["e"], emin=b["e"], dip=dip, sg=dip / math.hypot(a["es"], b["es"]), sr=dip / math.hypot(a["er"], b["er"]))

def window(t): return [r for r in t if 0.695 <= r["e"] <= 0.720]   # the draft's 0.695 <= eta <= 0.720, on the eta column

def main():
    O, N = load(OLD), load(NEW); tex = open(TEX).read(); tl = tex.splitlines()
    items = []   # (lines, old snippet, new snippet, definition, old computed, reproduced?, new computed, action)
    def add(old_s, new_s, definition, oc, ok, nc):
        lines = [i + 1 for i, l in enumerate(tl) if old_s in l]
        items.append((lines, old_s, new_s, definition, oc, ok, nc, "UNCHANGED" if old_s == new_s else "EDIT"))
    mo, mn = mean_dev_le04(O), mean_dev_le04(N)
    add(r"$+1.04\,\%$", f"$+{mn:.2f}\\,\\%$", "mean dev from KR, eta <= 0.4 (Roman script def.; 14 cells)", f"{mo:+.2f}", f"{mo:+.2f}" == "+1.04", f"{mn:+.2f}")
    xo, xn = max_dev_le04(O), max_dev_le04(N)
    add(r"never worse than $1.68\,\%$", f"never worse than ${xn:.2f}\\,\\%$", "max dev, eta <= 0.4 (the pi/8 cell, delta = 0)", f"{xo:.2f}", f"{xo:.2f}" == "1.68", f"{xn:.2f}")
    zo, zn = zoom_range(O), zoom_range(N)
    add(r"sits $0.5$--$1.5\,\%$ above", f"sits ${zn[0]:.1f}$--${zn[1]:.1f}\\,\\%$ above", "dev range of the N = 100 points, eta <= 0.15 (zoom XMAX)",
        f"{zo[0]:.2f} to {zo[1]:.2f}", f"{zo[0]:.1f}-{zo[1]:.1f}" == "0.5-1.5", f"{zn[0]:.2f} to {zn[1]:.2f}")
    do, dn = dev(at(O, 0.1122)), dev(at(N, 0.1122))
    add(r"sits $1.0\,\%$ above Kolafa", f"sits ${dn:.1f}\\,\\%$ above Kolafa", "dev at the worked-ladder cell, eta_rec = 0.1122",
        f"{do:.2f}", f"{do:.1f}" == "1.0", f"{dn:.2f}")
    k27 = "median of c_s_scatter_mass / c_s over eta_rec <= 0.69"
    f27 = lambda t: st.median([100 * r["sm"] / r["c"] for r in t if r["erec"] <= 0.69])
    add(r"That scatter, $0.27\,\%$", f"That scatter, ${f27(N):.2f}\\,\\%$", k27, f"{f27(O):.3f}", f"{f27(O):.2f}" == "0.27", f"{f27(N):.3f}")
    a52o, a52n = dev(at(O, 0.5236)), dev(at(N, 0.5236))
    add(r"ours sits $+2.6\,\%$", f"ours sits ${a52n:+.1f}\\,\\%$", "dev at eta = 0.5236 (L_0 = 7.5, delta = 0)", f"{a52o:+.2f}", f"{a52o:+.1f}" == "+2.6", f"{a52n:+.2f}")
    m_o, m_n = melt(O), melt(N)
    add(r"local maximum & $\eta = 0.700$", f"local maximum & $\\eta = {m_n['emax']:.3f}$", "eta of the local c_s maximum (melting script)",
        f"{m_o['emax']:.4f}", f"{m_o['emax']:.3f}" == "0.700", f"{m_n['emax']:.4f}")
    add(r"local minimum & $\eta = 0.715$", f"local minimum & $\\eta = {m_n['emin']:.3f}$", "eta of the local c_s minimum (melting script)",
        f"{m_o['emin']:.4f}", f"{m_o['emin']:.3f}" == "0.715", f"{m_n['emin']:.4f}")
    add(r"$c_s(0.700) - c_s(0.715) = 2.80$", f"$c_s({m_n['emax']:.3f}) - c_s({m_n['emin']:.3f}) = {m_n['dip']:.2f}$", "depth of the dip (melting script)",
        f"{m_o['dip']:.3f}", f"{m_o['dip']:.2f}" == "2.80", f"{m_n['dip']:.3f}")
    add(r"$\mathbf{9.5\sigma}$", f"$\\mathbf{{{m_n['sg']:.1f}\\sigma}}$", "dip / plotted (scaled) errors in quadrature", f"{m_o['sg']:.2f}", f"{m_o['sg']:.1f}" == "9.5", f"{m_n['sg']:.2f}")
    add(r"$18.2\sigma$ on the propagated", f"${m_n['sr']:.1f}\\sigma$ on the propagated", "dip / propagated errors in quadrature", f"{m_o['sr']:.2f}", f"{m_o['sr']:.1f}" == "18.2", f"{m_n['sr']:.2f}")
    co, cn = [r["chi"] for r in window(O)], [r["chi"] for r in window(N)]
    fmt = lambda v: f"{v:.1f}" if v < 10 else f"{v:.0f}"
    add(r"$\chi^2_{\mathrm{red}} = 1.8$--$14$", f"$\\chi^2_{{\\mathrm{{red}}}} = {fmt(min(cn))}$--${fmt(max(cn))}$", "chi2_red range, 0.695 <= eta <= 0.720",
        f"{fmt(min(co))}-{fmt(max(co))} ({len(co)} cells)", f"{fmt(min(co))}-{fmt(max(co))}" == "1.8-14", f"{fmt(min(cn))}-{fmt(max(cn))} ({len(cn)} cells)")
    out = lambda t: [r["sm"] for r in t if 0.66 <= r["e"] <= 0.765 and not (0.695 <= r["e"] <= 0.720)]
    io, oo = st.mean([r["sm"] for r in window(O)]), st.mean(out(O)); i_n, o_n = st.mean([r["sm"] for r in window(N)]), st.mean(out(N))
    add(r"$0.40$ inside", f"${i_n:.2f}$ inside", "mean c_s_scatter_mass, 0.695 <= eta <= 0.720", f"{io:.3f}", f"{io:.2f}" == "0.40", f"{i_n:.3f}")
    add(r"against $0.08$ outside", f"against ${o_n:.2f}$ outside", "mean c_s_scatter_mass, rest of the melting-figure range 0.66-0.765",
        f"{oo:.3f}", f"{oo:.2f}" == "0.08", f"{o_n:.3f}")

    print("### Paper 1 draft audit: every number from the canonical A1 v2 table, old -> new (printed BEFORE editing)\n")
    print("| tex line(s) | old text | new text | definition | old table gives | old text reproduced? | new table gives | action |")
    print("|---|---|---|---|---|---|---|---|")
    for lines, o, n, dfn, oc, ok, nc, act in items:
        print(f"| {', '.join(map(str, lines))} | `{o}` | `{n}` | {dfn} | {oc} | {'yes' if ok else '**NO**'} | {nc} | {act} |")
    print(f"\n'a factor five' (l. {[i+1 for i, l in enumerate(tl) if 'a factor five' in l]}): inside/outside = {io/oo:.2f} before, {i_n/o_n:.2f} after -> unchanged.")
    print(f"'both within our grid spacing of 0.005': |{m_n['emax']:.4f} - 0.702| = {abs(m_n['emax']-0.702):.4f}, "
          f"|{m_n['emin']:.4f} - 0.714| = {abs(m_n['emin']-0.714):.4f} -> still true.")
    print(f"'0.73' (the nofuse sigma check) in the draft: {sum('0.73' in l for l in tl)} occurrences -> not quoted, nothing to audit.")
    print(f"abstract 'eta = 0.0065 to 0.76': corrected range {min(r['e'] for r in N):.4f} to {max(r['e'] for r in N):.4f} -> unchanged.")
    r65 = at(N, 0.65); L0t = float(r65["L0"]) - float(r65["delta_sigma"]) / 2
    print(f"l. 101 thickness factor (t/2)/(L_0 - 2r) at eta = 0.65 with L_0,true: {100*(T.WALL_T/2)/(L0t-1.0):.3f} % (text 0.50 %); "
          f"at eta = 0.0065 delta = {float(at(N, 0.0065)['delta_sigma']):.0f} -> unchanged.")
    print("Cell labels NOT edited: 'eta = 0.1122' (ladder-line and slow-mode captions, ll. 211, 163) names the cell by its recorded "
          f"eta, as do the two 261001 figures drawn from its raw traces; its corrected eta is {at(N, 0.1122)['e']:.6f}. "
          "'24 densities with eta <= 0.69' (l. 222) is the canonical KR cut, applied to the recorded eta; the 24th cell's corrected "
          f"eta is {at(N, 0.69)['e']:.4f}.")

    print("\n### Per-cell deviations before and after (KR = the table's own KR column; D in units of c_s_err_scaled)\n")
    print("| eta_rec | eta_true | c_s - KR before | c_s - KR after | dev before [%] | dev after [%] | D before | D after | change in D | abs(D) |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    moved, grew = [], []
    for o, n in zip(O, N):
        assert abs(o["erec"] - n["erec"]) < 1e-9
        if o["kr"] is None or n["kr"] is None:
            print(f"| {n['erec']:.6f} | {n['e']:.6f} | n/a | n/a | n/a | n/a | n/a | n/a | n/a | no KR in table |"); continue
        Do, Dn = (o["c"] - o["kr"]) / o["es"], (n["c"] - n["kr"]) / n["es"]
        if abs(Dn - Do) > 0.5: moved.append((n["erec"], abs(Dn) < abs(Do)))
        if abs(Dn) > abs(Do) + 5e-3: grew.append((n["erec"], Do, Dn))
        print(f"| {n['erec']:.6f} | {n['e']:.6f} | {o['c']-o['kr']:+.5f} | {n['c']-n['kr']:+.5f} | {100*(o['c']-o['kr'])/o['kr']:+.3f} | "
              f"{100*(n['c']-n['kr'])/n['kr']:+.3f} | {Do:+.2f} | {Dn:+.2f} | {Dn-Do:+.2f} | {'smaller' if abs(Dn) < abs(Do) - 5e-3 else ('LARGER' if abs(Dn) > abs(Do) + 5e-3 else 'same')} |")
    print(f"\ncells with |change in D| > 0.5: {len(moved)} at eta_rec = {', '.join(f'{e:.3f}' for e, _ in moved)}; "
          f"all toward KR: {all(t for _, t in moved)}")
    print("cells where |D| grows (by > 0.005): " + ("; ".join(f"eta_rec = {e:.3f}: {a:+.2f} -> {b:+.2f}" for e, a, b in grew) or "none"))

    print("\n### Regeneration checks\n")
    ch = [n for o, n in zip(O, N) if o["c"] != n["c"] or o["eta"] != n["eta"]]
    same = [n for o, n in zip(O, N) if o["c"] == n["c"] and o["eta"] == n["eta"]]
    print(f"rows changed: {len(ch)} of {len(N)}; unchanged (delta = 0): L_0 = {', '.join(n['L0'] for n in same)}")
    idc = max(abs(n["c"] - o["c"] * float(n["L_eff_true"]) / (float(n["L_eff_true"]) + float(n["delta_sigma"]) / 2)) for o, n in zip(O, N))
    ide = max(abs(n["e"] - o["e"] * float(n["L0"]) / (float(n["L0"]) - float(n["delta_sigma"]) / 2)) for o, n in zip(O, N))
    print(f"identity c_s,new = c_s,old * L_eff,true/L_eff,rec: max |difference| {idc:.1e} (table prints 5 decimals)")
    print(f"identity eta_new = eta_old * L_0/(L_0 - delta/2): max |difference| {ide:.1e} (table prints 6 decimals)")

    zero = sum(1 for r in N if float(r["delta_sigma"]) == 0.0)
    sh = [(r["erec"], 100 * (1 - float(r["L_eff_true"]) / (float(r["L_eff_true"]) + float(r["delta_sigma"]) / 2)), 100 * (r["e"] / r["erec"] - 1)) for r in N]
    dil = max(sh, key=lambda s: s[1] if s[0] <= 0.16 else -1); top = max(sh, key=lambda s: s[1]); tope = max(sh, key=lambda s: s[2])
    print(f"\n### Numbers for the Methods paragraph\n\ndelta = 0 at {zero} of {len(N)} densities; c_s lowered by at most {dil[1]:.2f} % for "
          f"eta <= 0.16 (at eta_rec = {dil[0]:.4f}) and by at most {top[1]:.2f} % overall (at eta_rec = {top[0]:.4f}); eta raised by at most "
          f"{tope[2]:.2f} % (at eta_rec = {tope[0]:.4f}); {len(moved)} densities changed D by more than 0.5, all toward KR: "
          f"{all(t for _, t in moved)}; Center_X check quoted from methods sec. 14.1: {CENTER_X_CHECK} sigma.")

    if "--apply" in sys.argv:
        s = tex
        for lines, o, n, dfn, oc, ok, nc, act in items:
            if act != "EDIT": continue
            assert s.count(o) == len(lines) >= 1, o
            s = s.replace(o, n)
            print(f"edited ({len(lines)}x){'' if ok else ' -- old text was NOT reproduced by the old table; new text is the new table value'}: {o} -> {n}")
        anchor = "% TODO-source: validation/paper1_canonical_20260919.py; 260919_A1v2_final_cs_vs_eta.csv\n"
        assert s.count(anchor) == 1
        para = (anchor + "\n% ##CHRIS 2026-10-14: box-truncation paragraph (methods sec. 14), numbers from validation/paper1_draft_audit_20261014.py\n"
                "\\paragraph{Integer-pixel box width} The simulation sizes its box in whole pixels of $\\sigma/24$, so when\n"
                "$2L_0$ is not a whole number of pixels the box is shorter than requested by\n"
                "$\\delta = (48L_0 - \\lfloor 48L_0 \\rfloor)/24$, less than one pixel. The divider starts $L_0$ from one wall\n"
                "and oscillates about the centre of the shortened box, so each compartment is on average $L_0 - \\delta/2$\n"
                f"long, as the recorded box centre confirms at every density to ${CENTER_X_CHECK}\\,\\sigma$. For every $N = 100$\n"
                "density we therefore use $\\eta = N_s\\pi r^2/[(L_0 - \\delta/2)H]$ and\n"
                f"$L_{{\\mathrm{{eff}}}} = L_0 - 2r - t/2 - \\delta/2$; $\\delta = 0$ at {zero} of the {len(N)} densities, and the\n"
                f"correction lowers $c_s$ by at most ${dil[1]:.2f}\\,\\%$ for $\\eta \\le 0.16$ and ${top[1]:.2f}\\,\\%$ overall and raises\n"
                f"$\\eta$ by at most ${tope[2]:.2f}\\,\\%$. The correction was applied after the measurements, under a rule fixed\n"
                "before it was evaluated: the table is regenerated if any density's deviation from Kolafa--Rottner changes\n"
                f"by more than half that density's error bar, which happened at {['no','one','two','three','four','five','six','seven','eight'][len(moved)]} densities, all of which moved\n"
                "toward the equation of state. The larger systems of \\S\\ref{sec:finitesize} carry the same shortfall and are\n"
                "not yet corrected.\n"
                "% TODO-source: methods sec. 14 (rule, 4db8c9d), 14.1 (Center_X check; 615561c), 14.2; validation/paper1_boxtrunc_20261014.py\n")
        s = s.replace(anchor, para)
        open(TEX, "w").write(s); print("\napplied: edits + Methods paragraph written to paper1_draft.tex")

def a2():
    """##CHRIS 2026-10-02 (Task F3): every draft number and statement that depends on the A2 size ladder, old -> new,
    after the methods sec. 14.3 correction. Old = the dated copies (_pre_boxtrunc_261002), new = the regenerated tables.
    Numbers for the Methods sentence come from paper1_A2_boxtrunc_261002.py's own functions. --apply edits the EDIT rows."""
    import paper1_A2_boxtrunc_261002 as A2
    tex = open(TEX).read(); tl = tex.splitlines()
    rd = lambda n: {round(float(r["eta"]), 2): r for r in csv.DictReader(open(T.plot_path(n)))}
    P, NEW = rd("260917_A2_cs_vs_N_extrapolation_pre_boxtrunc_261002.csv"), rd("260917_A2_cs_vs_N_extrapolation.csv")
    items = []
    def add(old_s, new_s, definition, oc, ok, nc, act=None):
        lines = [i + 1 for i, l in enumerate(tl) if old_s.split("\n")[0] in l]
        items.append((lines, old_s, new_s, definition, oc, ok, nc, act or ("UNCHANGED" if old_s == new_s else "EDIT")))
    for eta in (0.02, 0.05):
        p, n = P[eta], NEW[eta]
        f = lambda r: f"{float(r['c_inf']):.4f} \\pm {float(r['c_inf_err']):.4f}"
        add(f(p), f(n), f"c_inf at eta = {eta} (260917 extrapolation)", f"{p['c_inf']} ± {p['c_inf_err']}",
            f(p) in tex, f"{n['c_inf']} ± {n['c_inf_err']}")
        g = lambda r: f"${float(r['KR']):.4f}$"
        add(g(p), g(n), f"KR at eta = {eta}", p["KR"], g(p) in tex, n["KR"])
        h = lambda r: f"(${float(r['dev_KR_sigma']):+.1f}\\sigma$)"
        add(h(p), h(n), f"(c_inf - KR)/sigma at eta = {eta}", p["dev_KR_sigma"], h(p) in tex, n["dev_KR_sigma"])
    # numbers of the sec. 14.3 result, for the Methods sentence
    pts = {t: {k: {g: A2.point(v, k[0], k[1], g) for g in "BA"} for k, v in A2.cells(t).items()} for t in A2.TABLES}
    maxD = max(abs(q["A"]["D"] - q["B"]["D"]) for t in pts for q in pts[t].values())
    src = pts[A2.TABLES[0]]; etas = sorted(e for e in {e for e, _ in src} if sum(1 for (x, _) in src if x == e) >= 3)
    F = {e: {g: A2.fit(src, e, g) for g in "BA"} for e in etas}
    dc = {e: (F[e]["A"]["cinf"] - F[e]["B"]["cinf"]) / F[e]["B"]["cinf_err"] for e in etas}
    db = {e: (F[e]["A"]["b"] - F[e]["B"]["b"]) / F[e]["B"]["b_err"] for e in etas}
    big = [e for e in etas if abs(db[e]) > 1]; maxc = max(abs(v) for v in dc.values()); maxb = max(abs(db[e]) for e in big)
    verb = "lowers" if all(db[e] < 0 for e in big) else "changes"
    print(f"sec. 14.3 numbers: max |change in D| = {maxD:.2f}; max |change in c_inf| = {maxc:.2f} sigma; "
          f"|change in b| > 1 sigma at eta = {big} ({', '.join(f'{db[e]:+.2f}' for e in big)}); max {maxb:.2f}\n")
    add("The larger systems of \\S\\ref{sec:finitesize} carry the same shortfall and are\nnot yet corrected.",
        "The same rule, pre-registered separately, was applied to\nthe larger systems of \\S\\ref{sec:finitesize}: "
        f"the correction moves no point of the size ladder by more than ${maxD:.2f}$\nof its error bar and no fitted "
        f"$c_\\infty$ by more than ${maxc:.2f}$ of its error, but it {verb} the coefficient of\n$1/\\sqrt N$ by up to "
        f"${maxb:.2f}$ of its error (at $\\eta = " + "$ and $".join(f"{e:.2f}" for e in big) + "$), so the size ladder was regenerated too.",
        "Methods paragraph: the A2 status sentence", "not corrected", True, f"corrected (sec. 14.3, REGENERATE)")
    add("For every $N = 100$\ndensity we therefore use", "For every density and system\nsize we therefore use",
        "Methods paragraph: scope of the correction", "N = 100 only", True, "all A1 and A2 cells")
    add("long, as the recorded box centre confirms at every density to",
        "long, as the recorded box centre confirms at every $N = 100$ density to",
        "Methods paragraph: the Center_X check covers the A1 runs only (sec. 14.1)", "every density", True, "every N = 100 density")
    add("% TODO-source: methods sec. 14 (rule, 4db8c9d), 14.1 (Center_X check; 615561c), 14.2; validation/paper1_boxtrunc_20261014.py",
        "% TODO-source: methods sec. 14 (rule, 4db8c9d), 14.1 (Center_X check; 615561c), 14.2; validation/paper1_boxtrunc_20261014.py;\n"
        "%   sec. 14.3 (A2 rule 0ddefa9, results 14.3.1); validation/paper1_A2_boxtrunc_261002.py",
        "source comment", "-", True, "+ sec. 14.3")
    zoomD = [(e, N, q["A"]["D"]) for (e, N), q in sorted(src.items()) if e <= 0.15 and N in (900, 1600)]
    worst = max(zoomD, key=lambda z: abs(z[2]))
    add("larger systems agree with it within their errors.", "larger systems agree with it within their errors.",
        "zoom caption: D (A geometry, plotted sigma) of the N = 900/1600 points shown (eta <= 0.15): "
        + ", ".join(f"{e:.2f}/N{N}: {d:+.2f}" for e, N, d in zoomD), "(claim)", all(abs(d) <= 1 for _, _, d in zoomD),
        f"largest |D| = {abs(worst[2]):.2f} at eta = {worst[0]:.2f}, N = {worst[1]}", act="FLAG -- claim not supported at 1 sigma; not edited (wording is the plan author's call)")
    print("### Paper 1 draft audit, A2-dependent text (methods sec. 14.3), old -> new (printed BEFORE editing)\n")
    print("| tex line(s) | old text | new text | definition | old table gives | old text reproduced / claim holds? | new table gives | action |")
    print("|---|---|---|---|---|---|---|---|")
    for lines, o, n, dfn, oc, ok, nc, act in items:
        print(f"| {', '.join(map(str, lines))} | `{o.replace(chr(10), ' ')}` | `{n.replace(chr(10), ' ')}` | {dfn} | {oc} | {'yes' if ok else '**NO**'} | {nc} | {act} |")
    print("\nNot edited: 'a size ladder to N = 2500 shows the offset is finite size' (abstract) and 'The N = 100 offset closes as the box "
          f"grows' (sec. finitesize) -- after the correction c_inf sits at {NEW[0.02]['dev_KR_sigma']} and {NEW[0.05]['dev_KR_sigma']} sigma "
          "from KR at eta = 0.02 and 0.05, as before.")
    if "--apply" in sys.argv:
        s = tex
        for lines, o, n, dfn, oc, ok, nc, act in items:
            if act != "EDIT": continue
            assert s.count(o) == 1, o; s = s.replace(o, n); print(f"edited: {o.splitlines()[0][:70]} ...")
        open(TEX, "w").write(s); print("\napplied to paper1_draft.tex")

def m():
    """##CHRIS 2026-10-02 (Task M, the plan author's decisions): M1 zoom caption, M2 the worked cell's label, M3 the stale
    error-bar wording in the zoom and overlay scripts. Every text change is printed old -> new BEFORE editing; --apply
    edits the tex rows (TEX); the script rows (SCRIPT) are applied exactly as printed and checked afterwards."""
    import paper1_A2_boxtrunc_261002 as A2
    tex = open(TEX).read(); tl = tex.splitlines()
    src = {k: A2.point(v, k[0], k[1], "A") for k, v in A2.cells(A2.TABLES[0]).items()}
    zd = sorted((e, N, q["D"]) for (e, N), q in src.items() if e <= 0.15 and N in (900, 1600))
    w = max(zd, key=lambda z: abs(z[2])); rest = max(abs(d) for e, N, d in zd if (e, N) != (w[0], w[1]))
    row = [r for r in csv.DictReader(open(T.plot_path("260919_A1v2_final_cs_vs_eta.csv"))) if abs(float(r["eta_rec"]) - 0.1122) < 1e-6][0]
    e4 = f"{float(row['eta']):.4f}"
    print(f"zoom points N = 900/1600 (A geometry, plotted sigma): " + ", ".join(f"eta {e:.2f} N {N}: D = {d:+.2f}" for e, N, d in zd))
    print(f"largest |D| {abs(w[2]):.2f} at eta = {w[0]:.2f}, N = {w[1]}; all others |D| <= {rest:.2f}; worked cell eta_true = {row['eta']} -> {e4}\n")
    zoom, over = "lowdensity_zoom_20260917.py", "overlay_N100_vs_A2_20260915.py"
    new_bar = "error bars = 1σ propagated error of the through-origin slope (σ_ν = sd/√n per mass, × √χ²_red when > 1)"
    new_doc = "error bars are the 1 sigma\npropagated error of the through-origin slope (slope_with_errors: sd/sqrt(n) per mass, x sqrt(chi2_red) when > 1)"
    items = [  # (kind, where, old, new, why)
        ("TEX", "zoom caption (M1)", "larger systems agree with it within their errors.}",
         f"larger systems agree with it within their errors, except $N = {w[1]}$ at $\\eta = {w[0]:.2f}$ (${abs(w[2]):.0f}\\sigma$).}}",
         f"D = {w[2]:+.2f}; the other five points |D| <= {rest:.2f}"),
        ("TEX", "slow-mode caption (M2; first mention in the text)", "$\\eta = 0.1122$, with its five-period running mean",
         f"$\\eta = {e4}$ (recorded $0.1122$), with its five-period running mean", f"eta_true = {row['eta']}; figure title says {e4}"),
        ("TEX", "ladder-line caption (M2)", "One worked mass ladder, $\\eta = 0.1122$, nine divider masses",
         f"One worked mass ladder, $\\eta = {e4}$, nine divider masses", f"figure title says {e4}, legend gives both"),
        ("SCRIPT", f"{zoom}:107 footnote (M3)", "error bars = 1σ scatter of per-mass c_s", new_bar, "the plotted bar is slope_with_errors' err_scaled (:44)"),
        ("SCRIPT", f"{over}:99 footnote (M3)", "error bars = 1σ scatter of per-mass c_s", new_bar, "the plotted bar is slope_with_errors' err_scaled (:42)"),
        ("SCRIPT", f"{zoom}:113 printed header (M3)", "| η | N | c_s | ± scatter | KR | dev [%] |", "| η | N | c_s | ± err (propagated) | KR | dev [%] |", "same quantity"),
        ("SCRIPT", f"{zoom}:3-4 docstring (M3)", "error bars are the\n1 sigma scatter of the per-mass c_s.", new_doc + ".", "same quantity"),
        ("SCRIPT", f"{over}:4 docstring (M3)", "error bars are the 1 sigma scatter of the per-mass c_s.", new_doc + ".", "same quantity"),
    ]
    print("### Task M text changes, old -> new (printed BEFORE editing)\n")
    print("| kind | where | tex line | old | new | why |\n|---|---|---|---|---|---|")
    for k, where, o, n, why in items:
        ln = ", ".join(str(i + 1) for i, l in enumerate(tl) if o.split("\n")[0] in l) if k == "TEX" else "-"
        print(f"| {k} | {where} | {ln} | `{o.replace(chr(10), ' ')}` | `{n.replace(chr(10), ' ')}` | {why} |")
    if "--apply" in sys.argv:
        s = tex
        for k, where, o, n, why in items:
            if k != "TEX": continue
            assert s.count(o) == 1, o; s = s.replace(o, n); print(f"edited tex: {where}")
        open(TEX, "w").write(s)
    if "--check-scripts" in sys.argv:
        for k, where, o, n, why in items:
            if k != "SCRIPT": continue
            f = os.path.join(HERE, where.split(":")[0]); t = open(f).read()
            print(f"{where}: old {'ABSENT' if o not in t else 'STILL PRESENT'}, new {'present' if n in t else 'MISSING'}")

def p():
    """##CHRIS 2026-10-02 (Task P; plan author's DECISION: Paper 1 uses the Kolafa-Rottner 2006 rho_max = 0.90 fit, Eq. 7 with
    the coefficients of sec. 3.2, fitted to eta <= 0.7069; data are compared with it only for eta <= 0.69). P3: the four Z'
    values of the melting caption, where each is evaluated and whether that is inside the fit range. P2/P3/P4: every text
    change, old -> new, printed BEFORE editing. --apply edits the TEX rows; SCRIPT/DOC rows are applied exactly as printed
    and checked with --check-scripts."""
    import plot_speed_of_sound_edmd as sos
    tex = open(TEX).read(); tl = tex.splitlines()
    EMAX, ECMP = sos.KR2006_ETA_MAX, sos.KR2006_PLOT_ETA_MAX
    zan = lambda e: float(sos.dZ_kolafa_rottner_2006(np.array([e]))[0])
    zfd = lambda e, h=1e-5: float((sos.Z_kolafa_rottner_2006(np.array([e + h])) - sos.Z_kolafa_rottner_2006(np.array([e - h])))[0] / (2 * h))
    print(f"fit range: eta <= pi*0.90/4 = {EMAX:.6f} -> '{EMAX:.4f}' (KR2006_ETA_MAX); comparison cutoff KR2006_PLOT_ETA_MAX = {ECMP}\n")
    print("### P3: the Z' values quoted in the melting caption\n")
    print("| caption value | eta | Z' analytic (module) | Z' finite difference (h = 1e-5, as the melting script) | reproduced | eta <= fit range? | action |")
    print("|---|---|---|---|---|---|---|")
    keep = []
    for txt, e in (("48", 0.67), ("28", 0.69), ("-9", 0.700), ("+3386", 0.720)):
        a, d = zan(e), zfd(e); ins = e <= EMAX
        rep = (f"{a:+.0f}" if txt.startswith("+") else f"{a:.0f}") == txt
        if ins: keep.append((txt, e))
        print(f"| {txt} | {e:.3f} | {a:.4f} | {d:.4f} | {'yes' if rep else '**NO**'} | {'yes' if ins else '**NO** (beyond ' + format(EMAX, '.4f') + ')'} | "
              f"{'keep, labelled fluid-branch fit' if ins else 'REMOVE'} |")
    g = np.linspace(0.60, EMAX, 200001)
    Z, dZ = sos.Z_kolafa_rottner_2006(g), sos.dZ_kolafa_rottner_2006(g)
    c = sos.cs_adiabatic_2d_monatomic(Z, dZ, g, kbt=1, m=1); et = float(g[int(np.argmax(c))])
    print(f"\nKR's own c_s maximum on 0.60 <= eta <= {EMAX:.4f} (analytic Z'): eta = {et:.4f} -> caption says 'near eta = 0.683': "
          f"{'consistent' if f'{et:.3f}' == '0.683' else 'NOT 0.683'}\n")
    em, ec = f"{EMAX:.4f}", f"{ECMP:.2f}"
    fit = f"$\\rho_{{\\max}} = 0.90$ fit"
    vals = [(t.lstrip("+"), e) for t, e in keep]
    parts = [f"$Z' = {vals[0][0]}$ at $\\eta = {vals[0][1]:.2f}$"] + [f"${t}$ at ${e:.3f}$" if e == 0.700 else f"${t}$ at ${e:.2f}$" for t, e in vals[1:]]
    zt = (", ".join(parts[:-1]) + " and " + parts[-1]) if len(parts) > 1 else parts[0]
    items = [  # (kind, where, old, new, why)
        ("TEX", "Methods (sec. Model and apparatus): new paragraph after 'Integer-pixel box width' (P2: fit named once)",
         "%   sec. 14.3 (A2 rule 0ddefa9, results 14.3.1); validation/paper1_A2_boxtrunc_261002.py\n",
         "%   sec. 14.3 (A2 rule 0ddefa9, results 14.3.1); validation/paper1_A2_boxtrunc_261002.py\n\n"
         "% ##CHRIS 2026-10-02: the reference equation of state, named and bounded (plan author's decision, methods sec. 15.1)\n"
         "\\paragraph{Reference equation of state} We compare with the hard-disk equation of state of Kolafa and\n"
         f"Rottner~\\cite{{kolafa2006}}: the {fit} of their Eq.~(7), with the coefficients of their \\S3.2, which is\n"
         f"fitted to $\\eta \\le {em}$. $Z'$ is taken from the same fit, and data are compared with it only for\n"
         f"$\\eta \\le {ec}$.\n"
         "% TODO-source: hspist3/plot_speed_of_sound_edmd.py:620-640 (KR2006_COEFFICIENTS, KR2006_ETA_MAX); validation/paper1_kr_sanity_261002.py (PASS)\n",
         f"fit range {em} = pi*0.90/4; comparison cutoff {ec} = KR2006_PLOT_ETA_MAX"),
        ("TEX", "Fig. csvseta caption (P2)",
         "beyond $\\eta \\approx 0.69$ lie in that region, where no fluid reference exists; they are shown but\n"
         "are not counted as a deviation from Kolafa--Rottner.}",
         f"beyond $\\eta = {ec}$ lie in or near that region; Kolafa--Rottner ({fit}, fitted to\n"
         f"$\\eta \\le {em}$) is compared with data only for $\\eta \\le {ec}$, so those points are shown but are not\n"
         "counted as a deviation from it.}", "names the range exactly; 'no fluid reference' was imprecise below 0.7069"),
        ("TEX", "sec:ordering, first sentence (P2)",
         "Kolafa--Rottner is a \\emph{fluid} equation of state fitted below the transition, so a deviation\n"
         "from it above $\\eta \\approx 0.69$ measures the extrapolation, not the gas.",
         f"Kolafa--Rottner is a \\emph{{fluid}} equation of state; the {fit} used here is fitted to\n"
         f"$\\eta \\le {em}$ and compared with data only for $\\eta \\le {ec}$, so above that a deviation from it measures\n"
         "the fluid-branch fit, not the gas.", f"the fit reaches {em}, into the coexistence interval; 'fitted below the transition' was imprecise"),
        ("TEX", "melting caption (P3)",
         "Kolafa--Rottner is solid inside its fitted range and\n"
         "dashed beyond it, and is cut off at $\\eta = 0.705$ deliberately: its fitted $Z'$ runs $48$ at\n"
         "$\\eta = 0.67$, $28$ at $0.69$, $-9$ at $0.700$ and $+3386$ at $0.720$, so a $c_s$ built from it\n"
         "there is meaningless rather than merely uncertain. The negative $Z'$ already inside the stated\n"
         "fitted range is also why KR's own $c_s$ turns over near $\\eta = 0.683$.}",
         f"Kolafa--Rottner ({fit}, fitted to\n"
         f"$\\eta \\le {em}$) is solid where it is compared with data, $\\eta \\le {ec}$, and dashed from there to\n"
         f"$\\eta = 0.705$, inside its fit range; the fluid-branch fit gives {zt}.\n"
         "Beyond its fit range the fluid-branch fit is an extrapolation and its derivative carries no\n"
         "information. The negative $Z'$ already inside the fit range is also why KR's own $c_s$ turns over\n"
         f"near $\\eta = {et:.3f}$.}}", f"+3386 is at 0.720 > {em}: removed; the three values inside the range stay"),
        ("TEX", "zoom caption (P4)", "larger systems agree with it within their errors, except",
         "larger systems agree with it within $2\\sigma$, except", "plan author's decision 2"),
        ("UNCHANGED", "sec:ordering title, l. 342", "\\subsection{Above $\\eta \\approx 0.69$: the right comparison is not Kolafa--Rottner}",
         "(unchanged)", "states where the comparison stops, consistent with the decision"),
        ("UNCHANGED", "l. 398", "$c_s(\\eta)$ above $0.69$ is the melting transition", "(unchanged)", "about our data, not the KR range"),
        ("UNCHANGED", "l. 404 (Outlook)", "Above $\\eta \\approx 0.69$ the comparison is with", "(unchanged)", "where the comparison stops, consistent"),
        ("SCRIPT", "paper1_canonical_20260919.py:74 legend, main figure (first KR legend: fit named)",
         'label="Kolafa-Rottner 2006 (valid to eta = 0.69)")',
         f'label="Kolafa-Rottner 2006, rho_max = 0.90 fit (Eq. 7), fitted to eta <= {em};\\ncompared with data for eta <= {ec}")',
         "the first figure that draws the KR equation of state"),
        ("SCRIPT", "overlay_N100_vs_A2_20260915.py:68 legend", 'label="Kolafa–Rottner 2006 (valid to η ≈ 0.69)")',
         f'label="Kolafa–Rottner 2006 (fitted to η ≤ {em}; compared for η ≤ {ec})")', "range bounded"),
        ("SCRIPT", "paper1_melting_figure_20261002.py:68 legend", 'label="Kolafa–Rottner 2006 (fitted range, $\\\\eta \\\\leq 0.69$)")',
         f'label="Kolafa–Rottner 2006, compared with data ($\\\\eta \\\\leq {ec}$)")', f"0.69 is the comparison cutoff, not the fit range ({em})"),
        ("SCRIPT", "paper1_melting_figure_20261002.py:70 legend", 'label="Kolafa–Rottner, EXTRAPOLATED (no fluid branch here)")',
         f'label="same fit, not compared ($0.69 < \\\\eta \\\\leq 0.705$; fit range $\\\\eta \\\\leq {em}$)")', f"0.69-0.705 is inside the fit range: not an extrapolation of the fit"),
        ("SCRIPT", "paper1_melting_figure_20261002.py annotation",
         '"KR\'s fitted $Z\'$ changes sign near $\\\\eta = 0.70$\\nand diverges by $0.72$: $c_s$ from it is\\n"\n                "meaningless here, not merely uncertain"',
         f'"KR\'s fluid-branch fit: $Z\'$ changes sign near $\\\\eta = 0.70$,\\ninside its fit range; beyond $\\\\eta = {em}$ the fit is\\n"\n                "an extrapolation and its derivative carries no information"',
         "no statement about the fit beyond its range"),
        ("DOC", "writeup/papers/README.md, KolafaRottner2006 row",
         "| 10.1080/00268970600880574 | the reference EOS, valid to η ≈ 0.69 |",
         f"| 10.1080/00268970600967963 | the reference EOS: the ρ_max = 0.90 fit (Eq. 7, coefficients of §3.2), fitted to η ≤ {em}; compared with data for η ≤ {ec} |",
         "range named; DOI corrected from the paper's own title page [SOURCE: PDF in ZZZ_PAPER/]"),
    ]
    print("### P2 / P3 / P4: every text change, old -> new (printed BEFORE editing)\n")
    print("| kind | where | tex line | old | new | why |\n|---|---|---|---|---|---|")
    for k, where, o, n, why in items:
        ln = ", ".join(str(i + 1) for i, l in enumerate(tl) if o.split("\n")[0] in l) if k in ("TEX", "UNCHANGED") else "-"
        print(f"| {k} | {where} | {ln} | `{o.strip().replace(chr(10), ' ')}` | `{n.strip().replace(chr(10), ' ')}` | {why} |")
    if "--apply" in sys.argv:
        s = tex
        for k, where, o, n, why in items:
            if k != "TEX": continue
            assert s.count(o) == 1, o; s = s.replace(o, n); print(f"edited tex: {where}")
        open(TEX, "w").write(s)
    path = lambda k, where: (os.path.join(REPO, "0000_PLAN_OVERALL/paper1_speedofsound/writeup/papers/README.md") if k == "DOC"
                             else os.path.join(HERE, where.split(":")[0].split(" ")[0]))
    if "--apply-scripts" in sys.argv:       # the SCRIPT/DOC rows, applied from the same strings that were printed
        for k, where, o, n, why in items:
            if k not in ("SCRIPT", "DOC"): continue
            f = path(k, where); t = open(f).read(); assert t.count(o) == 1, (where, o)
            open(f, "w").write(t.replace(o, n)); print(f"edited {k.lower()}: {where}")
    if "--check-scripts" in sys.argv:
        for k, where, o, n, why in items:
            if k not in ("SCRIPT", "DOC"): continue
            t = open(path(k, where)).read()
            print(f"{where}: old {'ABSENT' if o not in t else 'STILL PRESENT'}, new {'present' if n in t else 'MISSING'}")

if __name__ == "__main__":
    a2() if "--a2" in sys.argv else (m() if "--m" in sys.argv else (p() if "--p" in sys.argv else main()))
