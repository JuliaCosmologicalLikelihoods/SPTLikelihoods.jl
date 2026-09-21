#!/usr/bin/env python3
"""Convert the official SPT-3G D1 release into numeric Julia-artifact inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tarfile
from pathlib import Path

import numpy as np
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--archive", type=Path)
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


def source_paths(data_source: Path, descriptor: dict[str, object]) -> list[Path]:
    root = data_source / "spt_candl_data"
    d1 = root / "SPT3G_D1_TnE_v0"
    paths = [
        d1 / "SPT3G_D1_TnE.yaml",
        d1 / str(descriptor["band_power_file"]),
        d1 / str(descriptor["covariance_file"]),
        d1 / "effective_frequencies.yaml",
        root / "foreground_templates" / "dl_tsz_from_agora.txt",
        root / "foreground_templates" / "dl_ksz_from_agora.txt",
        d1 / "beams_templates" / "cov_eigenmodes_300_4100.npz",
        d1 / "beams_templates" / "B_ell_300_4100_main_rc4.npz",
        d1 / "beams_templates" / "beam_modes_prior_covariance.txt",
    ]
    paths.extend(sorted((d1 / str(descriptor["window_functions_folder"])).glob("*.txt")))
    return paths


def write_array(path: Path, array: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, array)


def create_archive(source: Path, archive: Path) -> None:
    archive.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "w:xz", preset=9) as tar:
        for path in sorted(source.rglob("*")):
            if path.is_file():
                tar.add(path, arcname=path.relative_to(source))


def main() -> None:
    args = parse_args()
    root = args.data_source / "spt_candl_data"
    d1 = root / "SPT3G_D1_TnE_v0"
    descriptor_path = d1 / "SPT3G_D1_TnE.yaml"
    descriptor = yaml.safe_load(descriptor_path.read_text())

    if args.output.exists():
        if any(args.output.iterdir()):
            raise FileExistsError(f"output directory must be empty: {args.output}")
    else:
        args.output.mkdir(parents=True)

    spectra_info = descriptor["spectra_info"]
    spectrum_order = [next(iter(entry)) for entry in spectra_info]
    bins_per_spectrum = [int(next(iter(entry.values()))) for entry in spectra_info]
    ells = np.arange(2, 4096, dtype=np.int64)

    data = np.loadtxt(d1 / descriptor["band_power_file"], dtype=np.float64)
    covariance = np.loadtxt(d1 / descriptor["covariance_file"], dtype=np.float64)
    if data.shape != (sum(bins_per_spectrum),):
        raise ValueError(f"unexpected data-vector shape: {data.shape}")
    if covariance.shape != (len(data), len(data)):
        raise ValueError(f"unexpected covariance shape: {covariance.shape}")
    if not np.allclose(covariance, covariance.T, rtol=0.0, atol=1e-12):
        raise ValueError("released covariance is not symmetric")

    write_array(args.output / "data_vector.npy", data)
    write_array(args.output / "covariance.npy", covariance)
    write_array(args.output / "ells.npy", ells)

    windows_dir = d1 / descriptor["window_functions_folder"]
    for index, (spectrum, n_bins) in enumerate(zip(spectrum_order, bins_per_spectrum), start=1):
        path = windows_dir / f"{spectrum.replace(' ', '_')}_window_functions.txt"
        raw = np.loadtxt(path, dtype=np.float64)
        if raw.shape != (len(ells), n_bins + 1):
            raise ValueError(f"unexpected window shape for {spectrum}: {raw.shape}")
        if not np.array_equal(raw[:, 0].astype(np.int64), ells):
            raise ValueError(f"unexpected ell grid for {spectrum}")
        write_array(args.output / "windows" / f"{index:02d}.npy", raw[:, 1:])

    for name in ("tsz", "ksz"):
        source = root / "foreground_templates" / f"dl_{name}_from_agora.txt"
        raw = np.loadtxt(source, dtype=np.float64)
        if not np.array_equal(raw[:, 0].astype(np.int64), np.arange(len(raw))):
            raise ValueError(f"unexpected template ell grid: {source}")
        write_array(args.output / "templates" / f"{name}.npy", raw[ells, 1])

    with np.load(d1 / "beams_templates" / "cov_eigenmodes_300_4100.npz") as archive:
        beam_ells = archive["ell"]
        raw_modes = archive["modes"]
    beam_mask = np.isin(beam_ells, ells)
    if not np.array_equal(beam_ells[beam_mask], ells):
        raise ValueError("beam eigenmodes do not cover the D1 ell grid")
    selected_modes = raw_modes[np.tile(beam_mask, 3), :9]
    write_array(args.output / "beams" / "eigenmodes.npy", selected_modes.reshape(3, len(ells), 9))

    with np.load(d1 / "beams_templates" / "B_ell_300_4100_main_rc4.npz") as archive:
        main_ells = archive["ell"]
        main_mask = np.isin(main_ells, ells)
        if not np.array_equal(main_ells[main_mask], ells):
            raise ValueError("main beams do not cover the D1 ell grid")
        main_beams = np.stack([archive[frequency][main_mask] for frequency in ("90", "150", "220")])
    write_array(args.output / "beams" / "main_temperature.npy", main_beams)
    write_array(
        args.output / "beams" / "prior_covariance.npy",
        np.loadtxt(d1 / "beams_templates" / "beam_modes_prior_covariance.txt", dtype=np.float64),
    )

    effective_frequencies = yaml.safe_load((d1 / "effective_frequencies.yaml").read_text())
    metadata = {
        "artifact": "SPT3G_D1_TnE_v0",
        "source_revision": git_revision(args.data_source),
        "descriptor_sha256": sha256(descriptor_path),
        "source_sha256": {
            str(path.relative_to(args.data_source)): sha256(path)
            for path in source_paths(args.data_source, descriptor)
        },
        "ell_min": int(ells[0]),
        "ell_max": int(ells[-1]),
        "spectrum_order": spectrum_order,
        "bins_per_spectrum": bins_per_spectrum,
        "effective_frequencies_ghz": effective_frequencies,
        "arrays": {
            "data_vector": [len(data)],
            "covariance": list(covariance.shape),
            "windows": [len(spectrum_order), len(ells), bins_per_spectrum],
            "eigenmodes": [3, len(ells), 9],
            "main_temperature_beams": list(main_beams.shape),
        },
    }
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    if args.archive is not None:
        create_archive(args.output, args.archive)
        print(f"archive_sha256={sha256(args.archive)}")


if __name__ == "__main__":
    main()
