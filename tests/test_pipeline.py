import csv
import hashlib
import json
import numpy as np
import pytest
from tbg_et.cli import main
from tbg_et.data import load_dataset, load_csv
from tbg_et.demo import make_demo
from tbg_et.pipeline import analyze


def test_end_to_end_and_roundtrip(tmp_path):
    data = make_demo(grid_size=6, energy_points=801)
    result = analyze(data, tmp_path / "run", "AB-center", "AA-center", plots=False)
    assert result["input_kind"] == "synthetic"
    assert abs(result["decomposition"]["closure_error"]) < 1e-12
    loaded = load_dataset(tmp_path / "run/input.npz", tmp_path / "run/metadata.json")
    np.testing.assert_array_equal(loaded.dos, data.dos)
    manifest = json.loads((tmp_path / "run/manifest.json").read_text())
    for name, digest in manifest["sha256"].items():
        assert hashlib.sha256((tmp_path / "run" / name).read_bytes()).hexdigest() == digest
    with pytest.raises(ValueError, match="not empty"):
        analyze(data, tmp_path / "run", "AB-center", "AA-center", plots=False)


@pytest.mark.parametrize("problem", ["units", "negative", "duplicate", "missing_provenance"])
def test_bad_input_rejected(problem):
    data = make_demo(grid_size=6, energy_points=401)
    if problem == "units":
        data.metadata["dos_units"] = "arbitrary"
    elif problem == "negative":
        data.dos[0, 0] = -1
    elif problem == "duplicate":
        data.site_ids[1] = data.site_ids[0]
    else:
        data.metadata.pop("coupling_model")
    with pytest.raises(ValueError):
        data.validate()


def test_cli_demo_and_external_analysis(tmp_path):
    main(["demo", "--angles", "1.1", "--grid-size", "6", "--energy-points", "401", "--no-plots",
          "--output", str(tmp_path / "demo")])
    root = tmp_path / "demo/theta_1.1"
    main(["analyze", "--input", str(root / "input.npz"), "--metadata", str(root / "metadata.json"),
          "--reference", "AB-center", "--target", "SP-center", "--output", str(tmp_path / "analysis"), "--no-plots"])
    assert (tmp_path / "analysis/rates.csv").exists()


def test_csv_join_uses_ids_and_rejects_duplicate_energy(tmp_path):
    data = make_demo(grid_size=6, energy_points=401)
    (tmp_path / "metadata.json").write_text(json.dumps(data.metadata))
    (tmp_path / "sites.csv").write_text("site_id,registry,x_nm,y_nm,lambda_ev,coupling_ev\na,AB,0,0,1,0.001\nb,AA,1,0,0.8,0.002\n")
    ldos = "site_id,energy_ev,dos_states_per_ev_nm2\nb,1,4\na,0,1\nb,0,3\na,1,2\na,-1,2\nb,-1,4\n"
    (tmp_path / "ldos.csv").write_text(ldos)
    paths = [tmp_path / "sites.csv", tmp_path / "ldos.csv", tmp_path / "metadata.json"]
    loaded = load_csv(*paths)
    np.testing.assert_array_equal(loaded.dos, [[2, 1, 2], [4, 3, 4]])
    (tmp_path / "ldos.csv").write_text(ldos + "a,0,1\n")
    with pytest.raises(ValueError, match="duplicate energy"):
        load_csv(*paths)


def test_figures_render(tmp_path):
    analyze(make_demo(grid_size=6, energy_points=401), tmp_path / "fig", "AB-center", "AA-center")
    assert (tmp_path / "fig/overview.png").read_bytes().startswith(b"\x89PNG")
    assert (tmp_path / "fig/overview.pdf").read_bytes().startswith(b"%PDF")
