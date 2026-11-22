# Five Köppen–Geiger climate zones for regional drought analysis in China

[中文处理说明](METHODS.zh-CN.md)

This dataset provides a fixed five-zone classification on a 0.1° ERA5-Land grid for a soil-moisture drought analysis covering 1981–2025. The climate zones represent the 1991–2020 climatology. The processed GeoTIFF is the layer used for the regional statistics in the study.

## Files

| File | Contents |
| --- | --- |
| [kg_1991_2020_china_0p1_zones.tif](kg_1991_2020_china_0p1_zones.tif) | Final five-zone raster used in the study |
| [kg_1991_2020_china_0p1_unmerged.tif](kg_1991_2020_china_0p1_unmerged.tif) | Five-class raster on the same grid and analysis mask, before patch consolidation |
| [class_codes.csv](class_codes.csv) | Raster values, Köppen–Geiger group letters, names in English and Chinese, and source code ranges |
| [rebuild_zones.py](rebuild_zones.py) and [requirements.txt](requirements.txt) | Python code and dependencies for rebuilding and verifying the final raster |
| [validation.json](validation.json) | Grid and area checks, counts by class, and raster checksums |

Both GeoTIFFs have one unsigned 8-bit band and use WGS 84 (EPSG:4326). They contain 617 columns and 355 rows, with 0.1° pixels. Outer bounds are 73.45°E to 135.15°E and 18.15°N to 53.65°N. Valid values are 1 (A, Tropical), 2 (B, Arid), 3 (C, Temperate), 4 (D, Cold) and 5 (E, Polar). Value 0 is NoData. In a GIS, display the band as unique values and use nearest-neighbour resampling for these categorical data.

## Source and attribution

Beck, H. E., McVicar, T. R., Vergopolan, N., Berg, A., Lutsko, N. J., Dufour, A., Zeng, Z., Jiang, X., van Dijk, A. I. J. M., & Miralles, D. G. (2023). High-resolution (1 km) Köppen–Geiger maps for 1901–2099 based on constrained CMIP6 projections. *Scientific Data*, 10, 724. [Article](https://doi.org/10.1038/s41597-023-02549-6)

The climate input is `1991_2020/koppen_geiger_0p1.tif` in `koppen_geiger_tif.zip`, Figshare record 21789074, archive version 2 (14 January 2026). The archive's MD5 is `7fc2f5a15d4f5fe0ce59c9a9b502aa09`. Its 1991–2020 0.1° raster matches the study's local source in class values and georeferencing. The two climate rasters here are adaptations of the Beck et al. dataset, shared with attribution under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Please cite Beck et al. and link to this repository when reusing the processed layer. [Source archive](https://doi.org/10.6084/m9.figshare.21789074) · [GloH2O data and licence](https://www.gloh2o.org/koppen/)

ERA5-Land monthly soil moisture provided the target grid and contributed to the definition of the fixed analysis mask. The two GeoTIFFs encode climate class values. [Copernicus Climate Data Store](https://doi.org/10.24381/cds.68d2bb30) · [Muñoz-Sabater et al. (2021)](https://doi.org/10.5194/essd-13-4349-2021)

## Processing

1. **Align and mask.** Align the global 0.1° classification to the study's 0.1° ERA5-Land grid using nearest-neighbour resampling because the grid centres differ by approximately 0.05°. Retain the intersection of the study's fixed soil-moisture coverage and valid climate classes. This coverage includes mainland China, Hainan and valid cells in Taiwan, Hong Kong and Macao. Remote South China Sea islands are outside the analysis mask. For exact reproduction, valid cells are encoded by positive values in the supplied unmerged raster.
2. **Recode.** Map the source codes into their five Köppen–Geiger groups: 1–3 to A, 4–7 to B, 8–16 to C, 17–28 to D and 29–30 to E. Store the groups as integer values 1–5, respectively. Sum individual spherical grid-cell areas, using an Earth radius of 6,371.0088 km, to determine actual patch areas.
3. **Identify patches and retain cores.** Label eight-neighbour connected patches within each class. Retain every original patch of at least 25,000 km² and the original largest patch of each class regardless of area.
4. **Assign smaller patches.** In each round, assign an unassigned patch to the retained core with which it shares the longest combined boundary. The boundary includes shared edges with patches already assigned to that core. Evaluate all assignments using the previous round's state and apply them simultaneously. Edge lengths along latitude parallels vary with latitude. Use shared sides to measure contact, so a corner contributes zero length. Shared lengths rounded to eight decimal places in km are compared first. Ties favour a larger original core area and then the smaller deterministic patch index. Assignments stop when no further eligible shared edges exist, and unassigned patches keep their original classes.
5. **Check the outcome.** The 10,000, 25,000 and 50,000 km² trials yielded 32, 23 and 25 patches, respectively, from 559 original patches. The intermediate threshold had the fewest patches among these trials and reassigned less area than the 50,000 km² trial. It reassigned 4,899 cells (488,837.88 km², or 5.17% of the analysis area). Eight final patches below 25,000 km² are isolated by the shared-edge criterion. The valid footprint stayed at 95,859 cells and 9,456,590.82 km². Class areas changed, including a 9.34% decrease for Cold and a 10.76% increase for Polar. The consolidated layer is an adaptation for regional aggregation, with altered local class boundaries.

## Rebuild and verify

Use Python 3.11. Run the following from the downloaded directory. The supplied unmerged GeoTIFF contains the aligned five-class grid and the fixed valid mask.

```sh
python -m pip install -r requirements.txt
python rebuild_zones.py --verify
```

To independently verify the global-source alignment and five-class recoding on that mask, extract the 1991–2020 0.1° raster from the source archive and supply it as well.

```sh
python rebuild_zones.py --source 1991_2020/koppen_geiger_0p1.tif --verify
```

Both commands generate `rebuilt_zones.tif` and report the processing checks. Verification compares the aligned class values and final class values as applicable, along with the grid coordinates and CRS of the final raster. GeoTIFF binary checksums can differ when metadata or software versions differ.
