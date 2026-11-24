"""The coded terms of the production A_a rate, one tag per term (design section 1).

Production A_a rate = vacuum_bssn.primary_l2_rhs (L2) + v55_pirk_adapter.primary_l3_with_matter (L3), where
  L2: Aa = -(DrDr_alpha - lap_alpha/3) + alpha*(Rr_r - R/3), with DrDr_alpha, lap_alpha, Rr_r = -I*radial_bracket
      and R = -I*scalar_bracket coded in vacuum_bssn.geometry_terms (I = inverse_physical_radial = X**2/a);
  L3: Aa = beta*d(Aa) + alpha*K*Aa + KO(Aa) - (16 pi/3)*alpha*(pr - pt).
Rule, fixed before any split run: the named variables of the code (a, b, X, alpha, K, Aa, beta, Lambda, r, ap, bp,
app, bpp, chip, chipp, Lambda_p, metric_lambda, inverse_physical_radial, alphap, alphapp, pr, pt and the KO
filter) are kept whole; every parenthesised sum written inline in those expressions is distributed, so each
additive term of the code becomes one tag carrying its outer factors. Nothing is combined. A tag is "origin" when
its term carries an explicit 1/r (written /r, /(r*b), / r**2 through metric_lambda, or metric_lambda * r), else
"bulk". The lapse piece A_a*(d alpha) is its own tag. Atoms are evaluated in binary64 exactly as the code does;
the tag products are formed in extended precision.
"""
import numpy as np

WIDE = np.longdouble

ORIGIN = frozenset({
    "lap.alphap_r", "Rr.ap_rb", "Rr.metric_lambda", "Rr.metric_lambda_rbp", "Rr.chip_r",
    "R.bp_rb_3", "R.bp_rb_ab", "R.metric_lambda", "R.chip_r",
})
L2_NAMES = (
    "DrDr.alphapp", "DrDr.alphap_ap", "DrDr.alphap_chip",
    "lap.alphapp", "lap.alphap_ap", "lap.alphap_bp", "lap.alphap_chip", "lap.alphap_r",
    "Rr.app", "Rr.Lambda_p", "Rr.ap2", "Rr.bp2", "Rr.Lambda_ap", "Rr.ap_rb", "Rr.metric_lambda",
    "Rr.metric_lambda_rbp", "Rr.chipp", "Rr.chip_ap", "Rr.chip_bp", "Rr.chip_r",
    "R.app", "R.bpp", "R.Lambda_p", "R.ap2", "R.bp2", "R.bp_rb_3", "R.bp_rb_ab", "R.metric_lambda",
    "R.chipp", "R.chip2", "R.chip_ap", "R.chip_bp", "R.chip_r",
)
L3_NAMES = ("L3.beta_dAa", "L3.alphaKAa", "L3.ko", "matter.pr", "matter.pt")
NAMES = L2_NAMES + L3_NAMES + ("lapse",)

# Secondary reading, declared with the primary one: the rate functions' own top-level terms.
COARSE = {
    "DrDr_alpha": tuple(n for n in L2_NAMES if n.startswith("DrDr.")),
    "lap_alpha/3": tuple(n for n in L2_NAMES if n.startswith("lap.")),
    "alpha*Rr_r": tuple(n for n in L2_NAMES if n.startswith("Rr.")),
    "alpha*R/3": tuple(n for n in L2_NAMES if n.startswith("R.")),
    "beta*dAa": ("L3.beta_dAa",), "alpha*K*Aa": ("L3.alphaKAa",), "KO": ("L3.ko",),
    "matter": ("matter.pr", "matter.pt"), "lapse": ("lapse",),
}
COARSE_ORIGIN = frozenset({"lap_alpha/3", "alpha*Rr_r", "alpha*R/3"})


def l2_terms(vacuum, grid, s):
    d, d2 = vacuum._d, vacuum._d2
    r = grid.centers
    a, b, X = s.a, s.b, s.X
    ap, bp, app, bpp = d(grid, a, 1), d(grid, b, 1), d2(grid, a, 1), d2(grid, b, 1)
    chip = -0.5 * d(grid, X, 1) / X
    chipp = -0.5 * (d2(grid, X, 1) / X - (d(grid, X, 1) / X) ** 2)
    Lambda_p = d(grid, s.Lambda, -1)
    metric_lambda = (1.0 - a / b) / r**2
    inv = X**2 / a
    alphap, alphapp = d(grid, s.alpha, 1), d2(grid, s.alpha, 1)
    (r, a, b, ap, bp, app, bpp, chip, chipp, Lambda_p, ml, I, alpha, alphap, alphapp, Lam) = (
        np.asarray(v, dtype=WIDE) for v in
        (r, a, b, ap, bp, app, bpp, chip, chipp, Lambda_p, metric_lambda, inv, s.alpha, alphap, alphapp, s.Lambda))
    aI, third = alpha*I, WIDE(1)/3
    t = {
        # -DrDr_alpha = -I*(alphapp - alphap*(0.5*ap/a + 2*chip))
        "DrDr.alphapp": -I*alphapp,
        "DrDr.alphap_ap": I*alphap*(WIDE(.5)*ap/a),
        "DrDr.alphap_chip": I*alphap*(2*chip),
        # +lap_alpha/3 = (I/3)*(alphapp - alphap*(0.5*ap/a - bp/b - 2*chip - 2/r))
        "lap.alphapp": I*third*alphapp,
        "lap.alphap_ap": -I*third*alphap*(WIDE(.5)*ap/a),
        "lap.alphap_bp": I*third*alphap*(bp/b),
        "lap.alphap_chip": I*third*alphap*(2*chip),
        "lap.alphap_r": I*third*alphap*(2/r),
        # alpha*Rr_r = -alpha*I*radial_bracket
        "Rr.app": -aI*(WIDE(.5)*app/a),
        "Rr.Lambda_p": aI*(a*Lambda_p),
        "Rr.ap2": aI*(WIDE(.75)*(ap/a)**2),
        "Rr.bp2": -aI*(WIDE(.5)*(bp/b)**2),
        "Rr.Lambda_ap": aI*(WIDE(.5)*Lam*ap),
        "Rr.ap_rb": -aI*(ap/(r*b)),
        "Rr.metric_lambda": -aI*(2*ml),
        "Rr.metric_lambda_rbp": -aI*(2*ml*(r*bp/b)),
        "Rr.chipp": -aI*(4*chipp),
        "Rr.chip_ap": aI*(2*chip*(ap/a)),
        "Rr.chip_bp": -aI*(2*chip*(bp/b)),
        "Rr.chip_r": -aI*(2*chip*(2/r)),
        # -alpha*R/3 = (alpha*I/3)*scalar_bracket
        "R.app": aI*third*(WIDE(.5)*app/a),
        "R.bpp": aI*third*(bpp/b),
        "R.Lambda_p": -aI*third*(a*Lambda_p),
        "R.ap2": -aI*third*((ap/a)**2),
        "R.bp2": aI*third*(WIDE(.5)*(bp/b)**2),
        "R.bp_rb_3": aI*third*(2*3*bp/(r*b)),
        "R.bp_rb_ab": -aI*third*(2*(a/b)*bp/(r*b)),
        "R.metric_lambda": aI*third*(4*ml),
        "R.chipp": aI*third*(8*chipp),
        "R.chip2": aI*third*(8*chip**2),
        "R.chip_ap": -aI*third*(8*chip*(WIDE(.5)*ap/a)),
        "R.chip_bp": aI*third*(8*chip*(bp/b)),
        "R.chip_r": aI*third*(8*chip*(2/r)),
    }
    assert tuple(t) == L2_NAMES
    zero = np.zeros_like(alpha)
    t.update({k: zero for k in L3_NAMES})
    return t


def l3_terms(vacuum, vm, grid, s, scalars, matter):
    total = vm.total_matter_projection(grid, s, scalars, matter)
    dAa = vacuum._d(grid, s.Aa, 1)
    ko = vacuum.ko_dissipation(grid, s.Aa, 1)
    beta, dAa, alpha, K, Aa, ko, pr, pt = (np.asarray(v, dtype=WIDE) for v in
                                           (s.beta, dAa, s.alpha, s.K, s.Aa, ko, total["pr"], total["pt"]))
    c = WIDE(16)*WIDE(np.pi)/3
    t = {k: np.zeros_like(alpha) for k in L2_NAMES}
    t.update({
        "L3.beta_dAa": beta*dAa,
        "L3.alphaKAa": alpha*K*Aa,
        "L3.ko": ko,
        "matter.pr": -c*alpha*pr,
        "matter.pt": c*alpha*pt,
    })
    return t
