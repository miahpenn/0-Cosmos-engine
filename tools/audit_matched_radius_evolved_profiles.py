#!/usr/bin/env python3
"""Reproduce matched-physical-radius comparison from accepted-state artifacts.

No evolution is run. The script verifies archived report/profile hashes and
source/vendor pins, interpolates existing t=0.60 profiles to common radii, and
writes a machine-readable comparison. Research diagnostic only.
"""
from __future__ import annotations
import hashlib, json, os
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE="4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR="d6052d605673ce9d82cdc99b0df79855fa2ea215"
JOINT=Path(os.environ.get("JOINT_ARTIFACT_ROOT","reference/joint"))
FULL=Path(os.environ.get("FULL_ARTIFACT_ROOT","reference/full"))
OUT=Path(os.environ.get("MATCHED_RADIUS_OUT","results/matched-radius-full-profile-evolved-2026-10-25"))
NS=(60,120,240)
FIELDS=("profile_H","profile_M","profile_C_lambda","profile_Aa")
RADIUS=np.linspace(0.5,10.0,1000)
RADIUS_2_10=np.linspace(2.0,10.0,800)

def sha(path:Path)->str:
 h=hashlib.sha256()
 with path.open("rb") as f:
  for chunk in iter(lambda:f.read(1<<20),b""): h.update(chunk)
 return h.hexdigest()

def load_arm(root:Path,n:int,arm:str)->dict:
 d=root/f"N{n}"
 report_path=d/"report.json"; profiles_path=d/"profiles.npz"; sums_path=d/"SHA256SUMS.txt"
 if not all(p.is_file() for p in (report_path,profiles_path,sums_path)):
  raise RuntimeError(f"missing accepted-state artifact files under {d}")
 sums={}
 for line in sums_path.read_text().splitlines():
  parts=line.split()
  if len(parts)==2:sums[parts[1]]=parts[0]
 for p in (report_path,profiles_path):
  if sums.get(p.name)!=sha(p): raise RuntimeError(f"artifact checksum mismatch: {p}")
 report=json.loads(report_path.read_text())
 pin=report["source_pins"]["source_pin_check"]
 if pin.get("source_commit")!=SOURCE or pin.get("vendor_commit")!=VENDOR or pin.get("vendor_checkout")!=VENDOR:
  raise RuntimeError(f"source/vendor pin mismatch in {report_path}")
 if pin.get("protected_engine_vendor_dependency_diff_empty") is not True:
  raise RuntimeError(f"protected source diff not empty in {report_path}")
 if abs(float(report["configuration"]["target_time"])-0.60)>1e-12:
  raise RuntimeError(f"wrong accepted time in {report_path}")
 z=np.load(profiles_path,allow_pickle=False)
 if arm=="full":
  prefix=f"N{n}_full_continuous_metric_profile_final_"
 else:
  prefix=f"N{n}_joint_integral_and_cmc_row_final_"
 out={}
 for name in ("r",)+FIELDS:
  key=prefix+name
  if key not in z.files: raise RuntimeError(f"missing {key} in {profiles_path}")
  out[name]=np.asarray(z[key],dtype=np.float64)
 if not all(np.all(np.isfinite(v)) for v in out.values()):
  raise RuntimeError(f"nonfinite profile in {profiles_path}")
 if np.any(np.diff(out["r"])<=0): raise RuntimeError(f"nonmonotone radii in {profiles_path}")
 return out

def metrics(profiles:dict)->dict:
 out={}
 for f in FIELDS:
  y=np.interp(RADIUS,profiles["r"],profiles[f])
  y2=np.interp(RADIUS_2_10,profiles["r"],profiles[f])
  out[f]={"max_abs_r_0p5_to_10":float(np.max(np.abs(y))),
          "rms_r_0p5_to_10":float(np.sqrt(np.mean(y*y))),
          "rms_r_2_to_10":float(np.sqrt(np.mean(y2*y2)))}
 return out

def main()->int:
 arms={a:{n:load_arm(root,n,a) for n in NS} for a,root in (("partial_joint",JOINT),("full_continuous_profile",FULL))}
 result={
  "schema":"matched-radius-evolved-profile-crosswalk/v1",
  "status":"completed_offline_diagnostic",
  "source_commit":SOURCE,"vendor_commit":VENDOR,"accepted_time":0.60,
  "method":{"common_radius_window":[0.5,10.0],"common_samples":len(RADIUS),
            "rms_inner_window":[2.0,10.0],"rms_inner_samples":len(RADIUS_2_10),
            "interpolation":"linear in physical radius","comparison_by_cell_index":False},
  "artifacts_verified":True,
  "full_continuous_profile":{str(n):metrics(arms["full_continuous_profile"][n]) for n in NS},
  "partial_joint":{str(n):metrics(arms["partial_joint"][n]) for n in NS},
  "fixed_radius_H_full_profile":{},
  "acceptance":{"all_input_hashes_verified":True,"source_vendor_pins_match":True,
    "protected_source_diff_empty":True,"all_profiles_finite":True,
    "production_files_changed":False,"full_matter_coupled_propagation_identity_closed":False,
    "total_stress_energy_conservation_closed":False,"scientific_admission":False}
 }
 for r in (1.0,2.0,3.0,5.0,8.0):
  result["fixed_radius_H_full_profile"][str(int(r))]={str(n):float(np.interp(r,arms["full_continuous_profile"][n]["r"],arms["full_continuous_profile"][n]["profile_H"])) for n in NS}
 OUT.mkdir(parents=True,exist_ok=True)
 path=OUT/"summary.json";path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
 print(json.dumps({"status":result["status"],"verified":result["artifacts_verified"],
  "full_profile_H_max_abs":{n:result["full_continuous_profile"][str(n)]["profile_H"]["max_abs_r_0p5_to_10"] for n in NS},
  "partial_H_max_abs":{n:result["partial_joint"][str(n)]["profile_H"]["max_abs_r_0p5_to_10"] for n in NS},
  "scientific_admission":False},indent=2))
 return 0
if __name__=="__main__":raise SystemExit(main())
