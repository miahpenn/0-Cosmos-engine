"""Extended-B: diagnostic-only accepted-metric recovery tests."""
from __future__ import annotations
import math
from types import SimpleNamespace
import numpy as np
import pytest
from engine.matter_system import BSSNMetricSlice, ConservedSpecies, Species, _metric_arrays, primitives, geometry_metric_derivatives

def _slice(n=4, alpha=0.9, a_scale=1.0):
    r=1.0+np.arange(n,dtype=float); ones=np.ones(n); zeros=np.zeros(n)
    return BSSNMetricSlice(r=r,a=a_scale*ones.copy(),b=ones.copy(),X=ones.copy(),
        alpha=alpha*ones.copy(),beta=zeros.copy(),Aa=zeros.copy(),K=zeros.copy(),Lambda=zeros.copy(),B=zeros.copy())

def _radiation_state(metrics,E=1e-6,S_abs=0.2e-6):
    energy=np.empty(len(metrics)); momentum=np.empty(len(metrics))
    for i,m in enumerate(metrics):
        Sr=S_abs/math.sqrt(m.gamma_rr_inv)
        energy[i]=m.sqrt_gamma*E; momentum[i]=m.sqrt_gamma*Sr
    return ConservedSpecies(np.zeros(len(metrics)),energy,momentum)

def _ratio(U,metrics,i=0):
    E=U.energy_t[i]/metrics[i].sqrt_gamma
    Sr=U.momentum_r[i]/metrics[i].sqrt_gamma
    return math.sqrt(max(metrics[i].gamma_rr_inv*Sr*Sr,0.0))/max(abs(E),1e-300)

def test_recovery_metrics_default_is_stage_metric():
    arrays=_metric_arrays(_slice()); U=_radiation_state(arrays)
    p0=primitives(arrays,U,Species.RADIATION)
    p1=primitives(arrays,U,Species.RADIATION,recovery_metrics=None)
    p2=primitives(arrays,U,Species.RADIATION,recovery_metrics=arrays)
    assert [p.rho for p in p0] == pytest.approx([p.rho for p in p1])
    assert [p.rho for p in p0] == pytest.approx([p.rho for p in p2])

def test_stage_cone_failure_is_recoverable_on_accepted_metric():
    accepted=_metric_arrays(_slice(a_scale=1.0)); stage=_metric_arrays(_slice(a_scale=0.25))
    U=_radiation_state(accepted,E=1e-6,S_abs=0.95e-6)
    assert _ratio(U,accepted)<1.0 and _ratio(U,stage)>1.0
    with pytest.raises(ValueError,match=r"violates E>=\|S\|"):
        primitives(stage,U,Species.RADIATION)
    assert all(p.rho>0 for p in primitives(stage,U,Species.RADIATION,recovery_metrics=accepted))

def test_evolve_species_keeps_stage_metric_for_faces_and_sources(monkeypatch):
    import engine.matter_system as ms
    stage=_slice(alpha=0.9); accepted=_slice(alpha=0.6)
    stage_metrics=_metric_arrays(stage); accepted_metrics=_metric_arrays(accepted)
    U=_radiation_state(accepted_metrics,E=1e-6,S_abs=0.1e-6)
    md=geometry_metric_derivatives(stage,lambda v,p:np.gradient(v))
    seen={"recovery_alpha":[],"face_alpha":[],"source_alpha":[]}
    real_primitives=ms.primitives
    def prim_spy(metrics,state,species,recovery_metrics=None):
        if recovery_metrics is not None: seen["recovery_alpha"]=[m.alpha for m in recovery_metrics]
        if recovery_metrics is None: return real_primitives(metrics,state,species)
        return real_primitives(metrics,state,species,recovery_metrics=recovery_metrics)
    def face_spy(metric,ql,qr,species):
        seen["face_alpha"].append(metric.alpha); return np.zeros(3)
    def boundary_spy(metric,prim):
        seen["face_alpha"].append(metric.alpha); return np.zeros(3)
    def source_spy(metric,prim,*args):
        seen["source_alpha"].append(metric.alpha); return 0.0,0.0
    monkeypatch.setattr(ms,"primitives",prim_spy)
    monkeypatch.setattr(ms,"_hll_flux",face_spy)
    monkeypatch.setattr(ms,"_radiation_flux",boundary_spy)
    monkeypatch.setattr(ms,"_radiation_source",source_spy)
    out=ms.evolve_species(stage,md,U,Species.RADIATION,dt=1e-6,validate_physical_state=False,recovery_metric=accepted)
    assert out is not None
    assert seen["recovery_alpha"] == pytest.approx([0.6]*len(stage_metrics))
    assert seen["face_alpha"] and seen["source_alpha"]
    assert seen["face_alpha"] == pytest.approx([0.9]*len(seen["face_alpha"]))
    assert seen["source_alpha"] == pytest.approx([0.9]*len(seen["source_alpha"]))

def test_validation_override_is_radiation_only(monkeypatch):
    from engine import true_cmc_pirk_kernel as tmod
    from engine.true_cmc_pirk_kernel import V55TrueCMCPIRKKernel
    import engine.matter_system as ms
    accepted=_slice(alpha=0.6); stage=_slice(alpha=0.9,a_scale=0.25)
    Urad=_radiation_state(_metric_arrays(accepted),E=1e-6,S_abs=0.95e-6)
    class Matter:
        dark_matter=ConservedSpecies(np.ones(4),np.ones(4),np.zeros(4))
        baryons=ConservedSpecies(np.ones(4),np.ones(4),np.zeros(4))
        radiation=Urad
    monkeypatch.setattr(tmod,"metric_slice_from_q",lambda grid,geometry:geometry)
    calls=[]; real=ms.primitives
    def spy(metrics,U,species,recovery_metrics=None):
        calls.append((species,recovery_metrics))
        if recovery_metrics is None: return real(metrics,U,species)
        return real(metrics,U,species,recovery_metrics=recovery_metrics)
    monkeypatch.setattr(ms,"primitives",spy)
    with pytest.raises(ValueError,match=r"violates E>=\|S\|"):
        V55TrueCMCPIRKKernel._validate_matter_state(object(),stage,Matter())
    calls.clear()
    V55TrueCMCPIRKKernel._validate_matter_state(object(),stage,Matter(),recovery_geometry=accepted)
    assert len(calls)==3 and calls[0][1] is None and calls[1][1] is None
    assert calls[2][0] is Species.RADIATION and calls[2][1] is not None
    assert calls[2][1][0].alpha == pytest.approx(0.6)

def test_diagnostic_off_rhs_delegates_to_unchanged_production_rhs(monkeypatch):
    from engine.production_kernel import V55ProductionKernel
    from engine.true_cmc_pirk_kernel import V55TrueCMCPIRKKernel
    sentinel=object(); seen=[]
    def rhs(self,state): seen.append((self,state)); return sentinel
    monkeypatch.setattr(V55ProductionKernel,"_rhs",rhs)
    kernel=V55TrueCMCPIRKKernel(); state=object()
    assert kernel.use_accepted_metric_for_predictor_radiation_recovery is False
    assert kernel._rhs(state) is sentinel and seen==[(kernel,state)]

def test_predictor_rhs_passes_explicit_recovery_metric(monkeypatch):
    from engine import true_cmc_pirk_kernel as tmod
    from engine import matter_rhs as mr
    from engine.production_kernel import ProductionState
    from engine.scalar_system import ScalarFields
    from engine.true_cmc_pirk_kernel import V55TrueCMCPIRKKernel
    import engine.scalar_system as scalar_mod
    import engine.v55_matter as v55m
    z=np.zeros(4); o=np.ones(4)
    stage=_slice(alpha=0.9,a_scale=0.25); accepted=_slice(alpha=0.6)
    scalars=ScalarFields(o,z,z,z,z,z); rad=_radiation_state(_metric_arrays(_slice()),E=1e-6,S_abs=0.1e-6)
    matter=SimpleNamespace(dark_matter=ConservedSpecies(o,o,z),baryons=ConservedSpecies(o,o,z),radiation=rad)
    state=ProductionState(grid=SimpleNamespace(cell_derivative_fourth=lambda vals,parity=1:np.gradient(vals)),
        geometry=stage,scalars=scalars,matter=matter,t=0.0,tau=0.0,e_folds=0.0)
    monkeypatch.setattr(tmod,"metric_slice_from_q",lambda grid,geometry:geometry)
    monkeypatch.setattr(v55m,"dm_density",lambda metric,matter:z)
    class Rhs:
        S=PS=D=PD=phi=Pi=z
    monkeypatch.setattr(scalar_mod,"scalar_rhs_arrays",lambda *a,**k:Rhs())
    seen=[]
    def species_spy(metric,md,U,species,**kwargs):
        if species is Species.RADIATION: seen.append(kwargs.get("recovery_metric"))
        return ConservedSpecies(np.zeros_like(U.rest),np.zeros_like(U.energy_t),np.zeros_like(U.momentum_r))
    monkeypatch.setattr(mr,"species_rhs",species_spy)
    V55TrueCMCPIRKKernel()._rhs(state,radiation_recovery_metric=accepted)
    assert seen==[accepted]


def test_total_fluid_projection_passes_recovery_metric_only_to_radiation(monkeypatch):
    import engine.v55_matter as vm
    stage = _slice(alpha=0.9)
    accepted = _slice(alpha=0.6)
    matter = SimpleNamespace(dark_matter=object(), baryons=object(), radiation=object())
    seen = []
    def project(metric, state, species, recovery_metric=None):
        seen.append((metric, species, recovery_metric))
        return {k: np.zeros(4) for k in ("rho", "pr", "pt", "j")}
    monkeypatch.setattr(vm, "project_species", project)
    vm.total_fluid_projection(stage, matter, radiation_recovery_metric=accepted)
    assert len(seen) == 3
    assert all(row[0] is stage for row in seen)
    assert seen[0][1] is Species.DARK_MATTER and seen[0][2] is None
    assert seen[1][1] is Species.BARYON and seen[1][2] is None
    assert seen[2][1] is Species.RADIATION and seen[2][2] is accepted


def test_production_kernel_lapse_solver_forwards_recovery_metric(monkeypatch):
    import engine.production_kernel as pk
    from engine.production_kernel import V55ProductionKernel
    accepted = _slice(alpha=0.6)
    calls = []
    def fake_solve(*args, **kwargs):
        calls.append(kwargs)
        return np.ones(4), 0.0
    monkeypatch.setattr(pk, "solve_cmc_lapse", fake_solve)
    V55ProductionKernel._solve_lapse(object(), object(), object(), object(),
                                     radiation_recovery_metric=accepted)
    assert calls[-1].get("radiation_recovery_metric") is accepted
    V55ProductionKernel._solve_lapse(object(), object(), object(), object())
    assert "radiation_recovery_metric" not in calls[-1]


def test_production_kernel_rhs_passes_recovery_metric(monkeypatch):
    import engine.production_kernel as pk
    from engine.production_kernel import ProductionState, V55ProductionKernel
    from engine.scalar_system import ScalarFields
    z, o = np.zeros(4), np.ones(4)
    stage = _slice(alpha=0.9, a_scale=0.25)
    accepted = _slice(alpha=0.6)
    state = ProductionState(
        grid=SimpleNamespace(cell_derivative_fourth=lambda vals, parity=1: np.gradient(vals)),
        geometry=stage, scalars=ScalarFields(o, z, z, z, z, z),
        matter=SimpleNamespace(
            dark_matter=ConservedSpecies(o,o,z),
            baryons=ConservedSpecies(o,o,z),
            radiation=_radiation_state(_metric_arrays(_slice())),
        ),
        t=0.0, tau=0.0, e_folds=0.0,
    )
    monkeypatch.setattr(pk, "metric_slice_from_q", lambda grid, geometry: geometry)
    monkeypatch.setattr(pk, "dm_density", lambda metric, matter: z)
    class Rhs:
        S = PS = D = PD = phi = Pi = z
    monkeypatch.setattr(pk, "scalar_rhs_arrays", lambda *a, **k: Rhs())
    seen = []
    def species_spy(metric, md, U, species, **kwargs):
        if species is Species.RADIATION:
            seen.append(kwargs.get("recovery_metric"))
        return ConservedSpecies(np.zeros_like(U.rest), np.zeros_like(U.energy_t), np.zeros_like(U.momentum_r))
    monkeypatch.setattr(pk, "species_rhs", species_spy)
    V55ProductionKernel()._rhs(state, radiation_recovery_metric=accepted)
    assert seen == [accepted]


def test_cmc_lapse_threads_recovery_metric_to_stress_and_target(monkeypatch):
    import engine.cmc_gauge as cg
    geometry = _slice(n=8)
    grid = SimpleNamespace(
        centers=geometry.r.copy(), n=8, dr=1.0,
        cell_derivative_fourth=lambda values, parity=1: np.gradient(values),
    )
    accepted = _slice(n=8, alpha=0.6)
    calls = []
    def total_spy(grid_arg, geometry_arg, scalars_arg, matter_arg, radiation_recovery_metric=None):
        calls.append(("total", radiation_recovery_metric))
        z = np.zeros(8)
        return {"rho": z, "pr": z, "pt": z, "j": z}
    def target_spy(grid_arg, geometry_arg, scalars_arg, matter_arg, outer_frac=0.2, radiation_recovery_metric=None):
        calls.append(("target", radiation_recovery_metric))
        return 0.0
    monkeypatch.setattr(cg, "total_matter_projection", total_spy)
    monkeypatch.setattr(cg, "target_kdot", target_spy)
    monkeypatch.setattr(cg, "solve_banded", lambda *a, **k: np.ones(8))
    cg.solve_cmc_lapse(
        grid, geometry, object(), object(),
        radiation_recovery_metric=accepted,
    )
    assert calls == [("total", accepted), ("target", accepted)]


def test_target_kdot_passes_recovery_metric_into_matter_rhs(monkeypatch):
    import engine.cmc_gauge as cg
    geometry = _slice(n=8)
    grid = SimpleNamespace(centers=geometry.r.copy(), n=8, dr=1.0)
    accepted = _slice(n=8, alpha=0.6)
    class Vacuum:
        @staticmethod
        def primary_l2_rhs(grid_arg, geometry_arg):
            return {"K": np.zeros(8)}
    monkeypatch.setattr(cg.adapter, "vendor_modules", lambda: (None, Vacuum(), None))
    calls = []
    def l3_spy(grid_arg, geometry_arg, scalars_arg, matter_arg, radiation_recovery_metric=None):
        calls.append(radiation_recovery_metric)
        return {"K": np.zeros(8)}
    monkeypatch.setattr(cg.adapter, "primary_l3_with_matter", l3_spy)
    out = cg.target_kdot(grid, geometry, object(), object(), radiation_recovery_metric=accepted)
    assert out == pytest.approx(0.0)
    assert calls == [accepted]
