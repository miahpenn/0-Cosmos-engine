"""One ruler run after PR #565 (pre-run entry in the coordination log 2026-11-23): R 12, xAR 5.9 at dr 0.05,
otherwise identical to research/frozen_s_threshold_followup_v1.py (run() unchanged, T 10000, estimator B).
python research/frozen_s_threshold_followup_5p9_v1.py [outdir]"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import frozen_s_threshold_followup_v1 as F

od = sys.argv[1] if len(sys.argv) > 1 else "results/frozen-s-threshold-followup-2026-11-23"
_, o = F.job(dict(R=12, xAR=5.9, dr=0.05))
c = F.cls(o)
txt = (f"# 5.9 ruler check (R 12, dr 0.05, T 10000): lambda_B {o['lam']:+.4f} [A {o['lamA']:+.4f}] class {c} "
       f"span {o['span']:.2f} floor {o['floor']:.2e} t_end {o['tend']:.0f}\n"
       f"# primary dr 0.1: +0.0035 G. Reading: " + ("growth -> (5.7, 5.9) is a measured cross-grid bracket" if c == "G"
       else "no growth -> the 5.9 edge stays provisional") + "; R 12 stays MIXED.\n")
open(os.path.join(od, "ruler_5p9.txt"), "w").write(txt); print(txt)
