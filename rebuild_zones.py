"""Rebuild the statistical Köppen–Geiger zones from the aligned five-class grid."""

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import Resampling, reproject
from scipy.ndimage import label


ROOT = Path(__file__).resolve().parent
BASELINE = ROOT / "kg_1991_2020_china_0p1_unmerged.tif"
REFERENCE = ROOT / "kg_1991_2020_china_0p1_zones.tif"
AREA_RADIUS_KM = 6371.0088
EDGE_RADIUS_KM = 6371.0
THRESHOLD_KM2 = 25_000
CLASSES = ("Tropical", "Arid", "Temperate", "Cold", "Polar")


def grid_areas(shape, transform):
    height, width = shape
    lat_edges = transform.f + transform.e * np.arange(height + 1)
    lon_edges = transform.c + transform.a * np.arange(width + 1)
    lat_centres = np.round((lat_edges[:-1] + lat_edges[1:]) / 2, 8)
    meridional = abs(np.diff(np.sin(np.deg2rad(lat_edges))))
    zonal = abs(np.diff(np.deg2rad(lon_edges)))
    return AREA_RADIUS_KM**2 * meridional[:, None] * zonal[None, :], lat_centres


def components(grid, areas):
    ids = np.zeros(grid.shape, dtype=np.int32)
    classes = [0]
    for code in range(1, 6):
        patches, count = label(grid == code, structure=np.ones((3, 3), dtype=bool))
        cells = patches > 0
        ids[cells] = patches[cells] + len(classes) - 1
        classes.extend([code] * count)
    sizes = np.bincount(ids.ravel(), weights=areas.ravel())
    sizes[0] = 0
    return ids, np.asarray(classes), sizes


def patch_borders(ids, latitudes):
    nominal = EDGE_RADIUS_KM * np.deg2rad(0.1)
    parallel_latitudes = np.deg2rad((latitudes[:-1] + latitudes[1:]) / 2)
    borders = defaultdict(lambda: defaultdict(float))
    for left, right, lengths in (
        (ids[:, :-1], ids[:, 1:], nominal),
        (ids[:-1], ids[1:], nominal * np.cos(parallel_latitudes[:, None])),
    ):
        touching = (left > 0) & (right > 0) & (left != right)
        weights = np.broadcast_to(lengths, left.shape)[touching]
        for first, second, length in zip(left[touching], right[touching], weights):
            borders[first][second] += length
            borders[second][first] += length
    return borders


def merge_patches(grid, areas, latitudes, threshold=THRESHOLD_KM2):
    ids, classes, sizes = components(grid, areas)
    adjacency = patch_borders(ids, latitudes)
    largest = [np.argmax(np.where(classes == code, sizes, -1)) for code in range(1, 6)]
    roots = np.where(sizes >= threshold, np.arange(len(sizes)), 0)
    roots[largest] = largest
    while True:
        pending = {}
        for patch in np.flatnonzero(roots[1:] == 0) + 1:
            shared = defaultdict(float)
            for neighbor, length in adjacency.get(patch, {}).items():
                if roots[neighbor]:
                    shared[roots[neighbor]] += length
            if shared:
                pending[patch] = max(
                    shared, key=lambda core: (round(shared[core], 8), sizes[core], -core)
                )
        if not pending:
            break
        for patch, core in pending.items():
            roots[patch] = core
    unassigned = (roots == 0) & (classes > 0)
    roots[unassigned] = np.flatnonzero(unassigned)
    return classes[roots[ids]].astype(np.uint8)


def align_global(source_path, profile, valid):
    fine = np.zeros(valid.shape, dtype=np.uint8)
    with rasterio.open(source_path) as source:
        if source.crs.to_epsg() != 4326 or not np.allclose(source.res, (0.1, 0.1)):
            raise ValueError("Supply the 0.1° WGS 84 global class raster for 1991–2020.")
        reproject(
            rasterio.band(source, 1), fine, src_nodata=0,
            dst_transform=profile["transform"], dst_crs=profile["crs"],
            dst_nodata=0, resampling=Resampling.nearest,
        )
    lookup = np.repeat(np.arange(6, dtype=np.uint8), (1, 3, 4, 9, 12, 2))
    return np.where(valid, lookup[fine], 0)


def summarize(before, after, areas):
    _, _, original_patches = components(before, areas)
    _, _, revised_patches = components(after, areas)
    valid = before > 0
    changed = before != after
    total_area = float(areas[valid].sum())
    changed_area = float(areas[changed].sum())
    return {
        "area_threshold_km2": THRESHOLD_KM2,
        "valid_cells": int(valid.sum()),
        "valid_area_km2": total_area,
        "original_patches": len(original_patches) - 1,
        "final_patches": len(revised_patches) - 1,
        "final_patches_below_threshold": int((revised_patches[1:] < THRESHOLD_KM2).sum()),
        "reassigned_cells": int(changed.sum()),
        "reassigned_area_km2": changed_area,
        "reassigned_area_percent": 100 * changed_area / total_area,
        "class_summary": {
            name: {
                "cells_before": int((before == code).sum()),
                "cells_after": int((after == code).sum()),
                "area_before_km2": float(areas[before == code].sum()),
                "area_after_km2": float(areas[after == code].sum()),
            }
            for code, name in enumerate(CLASSES, start=1)
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, help="Optional original 0.1° global GeoTIFF.")
    parser.add_argument("--output", type=Path, default=Path("rebuilt_zones.tif"))
    parser.add_argument("--verify", action="store_true", help="Compare with the published zones.")
    args = parser.parse_args()
    if args.output.resolve() in (BASELINE, REFERENCE):
        parser.error("Choose a separate output path for the rebuilt GeoTIFF.")

    with rasterio.open(BASELINE) as dataset:
        before, profile = dataset.read(1), dataset.profile
    if args.source:
        aligned = align_global(args.source, profile, before > 0)
        np.testing.assert_array_equal(aligned, before, err_msg="Baseline grid does not match.")
        before = aligned
    areas, latitudes = grid_areas(before.shape, profile["transform"])
    after = merge_patches(before, areas, latitudes)
    np.testing.assert_array_equal(after > 0, before > 0)
    if args.verify:
        with rasterio.open(REFERENCE) as expected:
            np.testing.assert_array_equal(after, expected.read(1), err_msg="Zone values differ.")
            if expected.transform != profile["transform"] or expected.crs != profile["crs"]:
                raise ValueError("Reference GeoTIFF has different georeferencing.")

    with rasterio.open(args.output, "w", **profile) as dataset:
        dataset.write(after, 1)
        dataset.update_tags(period="1991–2020", area_threshold_km2=THRESHOLD_KM2)
    print(json.dumps(summarize(before, after, areas), indent=2))


if __name__ == "__main__":
    main()
