"""
dev/raytheory/reference/make_taup_reference.py

TauP (ObsPy implementation) travel times used to validate dev/raytheory.
To compare the ray code and TauP on EXACTLY the same velocities, the
flexOPT DSM model is first written as a TauP .nd file by
    julia dev/raytheory/reference/writeND.jl dataInput/model/earth/PREM.psv.model premDSM.nd
(cubic zones sampled every 10 km, linear in between: < 1e-3 s effect), then
    python dev/raytheory/reference/make_taup_reference.py premDSM.nd taup_premDSM.csv
The CSV is committed, so TauP/ObsPy is not needed to run the dev scripts.
"""
import sys, csv
from obspy.taup.taup_create import build_taup_model
from obspy.taup import TauPyModel
import os, tempfile

nd, out = sys.argv[1], sys.argv[2]
tmp = tempfile.mkdtemp()
build_taup_model(nd, output_folder=tmp, verbose=False)
model = TauPyModel(model=os.path.join(tmp, os.path.splitext(os.path.basename(nd))[0] + ".npz"))

# our phase name -> TauP phase name (TauP needs the exact discontinuity depth)
phases = {"P": "P", "S": "S", "PP": "PP", "SS": "SS", "PS": "PS", "PcP": "PcP", "ScS": "ScS",
          "PKP": "PKP", "PKIKP": "PKIKP", "PKiKP": "PKiKP", "SKS": "SKS", "pP": "pP", "sP": "sP",
          "sS": "sS", "P410s": "P400s", "P660s": "P670s", "Pms": "Pms"}
depths = [0.0, 100.0, 500.0]
dists = [5.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 150.0, 160.0, 170.0]
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["phase", "taupName", "depth_km", "distance_deg", "time_s", "p_s_per_rad", "takeoff_deg"])
    for ours, tp in phases.items():
        for z in depths:
            for d in dists:
                for a in model.get_travel_times(source_depth_in_km=z, distance_in_degree=d, phase_list=[tp]):
                    if a.name != tp:     # TauP may return e.g. Pdiff for P
                        continue
                    w.writerow([ours, tp, z, d, "%.5f" % a.time, "%.6f" % a.ray_param, "%.4f" % a.takeoff_angle])
print("written", out)
