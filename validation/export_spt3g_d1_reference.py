#!/usr/bin/env python3
"""Export reproducible SPT-3G D1 full-likelihood references from candl."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import yaml



def _relative_to_repo(path: Path) -> str:
    """Return ``path`` relative to its enclosing ``spt_candl_data`` checkout.

    Reference metadata is committed, so it must not record an absolute path
    specific to the machine that produced it.
    """
    parts = Path(path).resolve().parts
    if "spt_candl_data" in parts:
        index = len(parts) - 1 - parts[::-1].index("spt_candl_data")
        return str(Path(*parts[index:]))
    return Path(path).name


CASES: dict[str, dict[str, float]] = {
    "foregrounds": {
        "TT_Poisson_90x90": 6.4,
        "TT_Poisson_90x150": 9.1,
        "TT_Poisson_90x220": 19.5,
        "TT_Poisson_150x150": 14.7,
        "TT_Poisson_150x220": 37.2,
        "TT_Poisson_220x220": 81.3,
        "TT_CIB_150x150": 4.8,
        "TT_CIB_150x220": 12.6,
        "TT_CIB_220x220": 42.5,
        "TT_CIBClustering_Alpha": 0.61,
        "TT_tSZ_Amp": 2.3,
        "TT_kSZ_Amp": 2.8,
        "TT_GalCirrus_Amp": 2.7,
        "TT_GalCirrus_Alpha": -2.46,
        "TT_GalCirrus_Beta": 1.53,
        "TE_PolGalDust_Amp": 0.15,
        "TE_PolGalDust_Alpha": -2.34,
        "TE_PolGalDust_Beta": 1.47,
        "EE_PolGalDust_Amp": 0.07,
        "EE_PolGalDust_Alpha": -2.37,
        "EE_PolGalDust_Beta": 1.55,
    },
    "instrument": {
        "Kappa": 2.1e-4,
        "T2P2_90": -0.0081,
        "T2P2_150": -0.0107,
        "T2P2_220": -0.0251,
        "beta_1": 0.32,
        "beta_2": -0.41,
        "beta_3": 0.23,
        "beta_4": -0.18,
        "beta_5": 0.11,
        "beta_6": -0.29,
        "beta_7": 0.07,
        "beta_8": -0.15,
        "beta_9": 0.21,
        "beta_pol_90": 0.34,
        "beta_pol_150": 0.58,
        "beta_pol_220": 0.76,
        "Tcal_ext150": 1.0031,
        "Tcal_rel90": 0.9962,
        "Tcal_rel220": 1.0084,
        "Ecal_ext150": 0.9927,
        "Ecal_rel90": 1.0053,
        "Ecal_rel220": 0.9894,
    },
}
CASES["combined"] = CASES["foregrounds"] | CASES["instrument"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candl-source", required=True, type=Path)
    parser.add_argument("--data-source", required=True, type=Path)
    parser.add_argument("--test-yaml", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git_revision(path: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def load_test_case(path: Path) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    record = yaml.safe_load(path.read_text())
    spectrum = np.loadtxt(path.parent / record["test_spectrum"])
    params = copy.deepcopy(record["param_values"])
    params["Dl"] = {
        spectrum_type: spectrum[:, column]
        for column, spectrum_type in enumerate(("ell", "TT", "TE", "EE", "BB", "pp", "kk"))
    }
    return params, params["Dl"]


def flattened_input(like: Any, params: dict[str, Any]) -> np.ndarray:
    return np.concatenate([np.asarray(params["Dl"][spectrum_type]) for spectrum_type in like.spec_types])


def transformation_key(index: int, transformation: Any) -> str:
    descriptor = getattr(transformation, "descriptor", type(transformation).__name__)
    slug = "".join(character.lower() if character.isalnum() else "_" for character in descriptor)
    return f"{index:02d}_{slug.strip('_')}"


def run_case(like: Any, params: dict[str, Any], *, export_stages: bool) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    theory = flattened_input(like, params)
    stages: dict[str, np.ndarray] = {"00_input": theory}

    for index, transformation in enumerate(like.data_model, start=1):
        theory = np.asarray(transformation.transform(theory, params))
        if export_stages:
            stages[transformation_key(index, transformation)] = theory

    binned = np.asarray(like.bin_model_specs(theory))
    residual = np.asarray(like.data_bandpowers) - binned
    data_chi2 = float(like.chi_square(params))
    prior_penalty = float(like.prior_logl(params))
    loglike = float(like.log_like(params))

    stages["final_unbinned"] = theory
    stages["binned_model"] = binned
    stages["residual"] = residual
    return stages, {
        "data_chi2": data_chi2,
        "prior_penalty": prior_penalty,
        "loglike": loglike,
    }


def source_hashes(data_source: Path) -> dict[str, str]:
    root = data_source / "spt_candl_data"
    d1 = root / "SPT3G_D1_TnE_v0"
    paths = [
        d1 / "SPT3G_D1_TnE.yaml",
        d1 / "SPT3G_D1_TnE_bdp.txt",
        d1 / "SPT3G_D1_TnE_cov.dat",
        d1 / "effective_frequencies.yaml",
        root / "foreground_templates" / "dl_tsz_from_agora.txt",
        root / "foreground_templates" / "dl_ksz_from_agora.txt",
        d1 / "beams_templates" / "cov_eigenmodes_300_4100.npz",
        d1 / "beams_templates" / "B_ell_300_4100_main_rc4.npz",
        d1 / "beams_templates" / "beam_modes_prior_covariance.txt",
    ]
    paths.extend(sorted((d1 / "windows").glob("*.txt")))
    return {str(path.relative_to(data_source)): sha256(path) for path in paths}


def main() -> None:
    args = parse_args()
    import sys

    sys.path[:0] = [str(args.candl_source), str(args.data_source)]
    import candl
    import spt_candl_data

    args.output.mkdir(parents=True, exist_ok=True)
    test_record = yaml.safe_load(args.test_yaml.read_text())
    base_params, cmb_spectra = load_test_case(args.test_yaml)
    likelihood_path = Path(spt_candl_data.SPT3G_D1_TnE_multifreq)
    like = candl.Like(str(likelihood_path), feedback=False)

    baseline_stages, baseline_metrics = run_case(like, base_params, export_stages=True)
    ells = cmb_spectra["ell"].astype(int)
    n_ell = len(ells)

    # 1. Input theory and templates
    np.savetxt(
        args.output / "cmb_spectra.txt",
        np.column_stack([cmb_spectra["ell"], cmb_spectra["TT"], cmb_spectra["TE"], cmb_spectra["EE"]]),
        fmt=["%d", "%.12e", "%.12e", "%.12e"],
        header="ell TT TE EE",
    )
    with open(args.output / "baseline_parameters.txt", "w") as f:
        for name, value in sorted(test_record["param_values"].items()):
            f.write(f"{name}\t{value}\n")

    data_root = args.data_source / "spt_candl_data"
    tsz_template = np.loadtxt(data_root / "foreground_templates" / "dl_tsz_from_agora.txt")[:, 1][ells]
    ksz_template = np.loadtxt(data_root / "foreground_templates" / "dl_ksz_from_agora.txt")[:, 1][ells]
    np.savetxt(
        args.output / "foreground_templates.txt",
        np.column_stack([tsz_template, ksz_template]),
        fmt="%.12e",
        header="tsz_template ksz_template",
    )

    beam_dir = data_root / "SPT3G_D1_TnE_v0" / "beams_templates"
    with np.load(beam_dir / "cov_eigenmodes_300_4100.npz") as archive:
        mask = np.isin(archive["ell"], ells)
        beam_modes = archive["modes"][np.tile(mask, 3), :9].reshape(3, n_ell, 9)
    with np.load(beam_dir / "B_ell_300_4100_main_rc4.npz") as archive:
        mask = np.isin(archive["ell"], ells)
        main_beams = np.stack([archive[frequency][mask] for frequency in ("90", "150", "220")])

    np.savetxt(
        args.output / "beam_main.txt",
        main_beams.T,
        fmt="%.12e",
        header="main_90 main_150 main_220",
    )
    # beam_modes: reshape (3, n_ell, 9) -> (n_ell, 27)
    beam_modes_flat = beam_modes.transpose(1, 0, 2).reshape(n_ell, 27)
    np.savetxt(
        args.output / "beam_modes.txt",
        beam_modes_flat,
        fmt="%.12e",
        header="beam_modes_27_cols (9 for 90, 9 for 150, 9 for 220)",
    )

    # 2. Intermediate stages and unbinned vectors (shape: (n_ell, 21))
    def _save_21(filename: str, array_1d: np.ndarray) -> None:
        reshaped = array_1d.reshape(21, n_ell).T
        np.savetxt(args.output / filename, reshaped, fmt="%.12e")

    _save_21("stage_preinstrument.txt", baseline_stages["10_ee_polarised_galactic_dust"])
    _save_21("stage_leakage.txt", baseline_stages["11_t2p_leakage"])
    _save_21("stage_beam.txt", baseline_stages["12_polarized_beam_model__rc3"])
    _save_21("baseline_final_unbinned.txt", baseline_stages["final_unbinned"])

    # Checkpoint ells for earlier foreground stages
    chk_indices = np.array([0, 50, 100, 250, 500, 1000, 1500, 2000, 2500, 3000, 3500, n_ell - 1])
    chk_stages = [
        "00_input", "01_super_sample_lensing", "02_aberration", "03_tt_poisson",
        "04_tt_cib_clustering_multi_amp", "05_tsz", "06_tt_tsz_cib_correlation",
        "07_ksz", "08_cirrus", "09_te_polarised_galactic_dust", "10_ee_polarised_galactic_dust",
    ]
    chk_data = []
    for s_idx, stage_name in enumerate(chk_stages):
        stage_arr = baseline_stages[stage_name].reshape(21, n_ell)
        for c_idx in chk_indices:
            ell_val = float(ells[c_idx])
            spec_vals = [stage_arr[p, c_idx] for p in range(21)]
            chk_data.append([float(s_idx), ell_val] + spec_vals)
    np.savetxt(
        args.output / "foreground_checkpoints.txt",
        np.array(chk_data),
        fmt=["%d", "%d"] + ["%.12e"] * 21,
        header="stage_index ell 21_spectrum_values",
    )

    # 3. Multipoint evaluations
    multipoint_parameters: dict[str, float] = {}
    metrics: dict[str, dict[str, float]] = {"baseline": baseline_metrics}
    binned_models = {"baseline": baseline_stages["binned_model"]}
    residuals = {"baseline": baseline_stages["residual"]}

    for name, updates in CASES.items():
        params = copy.deepcopy(base_params)
        params.update(updates)
        for parameter, value in params.items():
            if parameter != "Dl":
                multipoint_parameters[f"{name}__{parameter}"] = float(value)
        stages, case_metrics = run_case(like, params, export_stages=False)
        _save_21(f"{name}_final_unbinned.txt", stages["final_unbinned"])
        binned_models[name] = stages["binned_model"]
        residuals[name] = stages["residual"]
        metrics[name] = case_metrics

    with open(args.output / "multipoint_parameters.txt", "w") as f:
        for name, value in sorted(multipoint_parameters.items()):
            f.write(f"{name}\t{value}\n")

    case_names = ["baseline", "foregrounds", "instrument", "combined"]
    np.savetxt(
        args.output / "binned_models.txt",
        np.column_stack([binned_models[c] for c in case_names]),
        fmt="%.12e",
        header=" ".join(case_names),
    )
    np.savetxt(
        args.output / "residuals.txt",
        np.column_stack([residuals[c] for c in case_names]),
        fmt="%.12e",
        header=" ".join(case_names),
    )

    metadata = {
        "likelihood": "SPT3G_D1_TnE_v0 full multifrequency",
        "candl_revision": git_revision(args.candl_source),
        "spt_candl_data_revision": git_revision(args.data_source),
        # Record a repository-relative descriptor path, never this workstation's
        # absolute checkout location.
        "test_yaml": _relative_to_repo(args.test_yaml),
        "test_yaml_sha256": sha256(args.test_yaml),
        "declared_test_loglike": float(test_record["test_logl"]),
        "source_hashes": source_hashes(args.data_source),
        "ell_min": int(like.ells[0]),
        "ell_max": int(like.ells[-1]),
        "number_of_theory_multipoles": int(len(like.ells)),
        "spectrum_order": list(like.spec_order),
        "spectrum_types": list(like.spec_types),
        "bins_per_spectrum": [int(value) for value in like.N_bins],
        "number_of_bandpowers": int(len(like.data_bandpowers)),
        "base_parameters": test_record["param_values"],
        "multipoint_updates": CASES,
        "metrics": metrics,
        "baseline_stage_keys": list(baseline_stages),
    }
    (args.output / "reference.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
