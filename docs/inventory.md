# Academic folder inventory

Root: `E:\Projects Academia Ongoing\Ebola VHH (Hidde Lab)`

All of these stay on E:\. `config/paths.json` points at them. Original MATLAB scripts also mention `G:\LLSM_VSV_VHH`, `F:\LLSM_VSV_VHH`, and Dropbox copies of the same project.

## Paper and figure composites

| Path | What |
| --- | --- |
| `Ebolavirus VHHs v2 230104.docx` | Current draft. Figure 5 is the LLSM/cell-binding figure. |
| `FIGURE FILES/FINAL FIGURES 220317/` | 17 Mar 2022 freeze: Illustrator/PDF composites, per-panel `.fig`/EPS/TIFF, `Main&Suppl_Figs_Text.docx`. |
| `FIGURE FILES/FINAL FIGURES 220317/EboV_VHH_Main_Fig.ai` | Original main-figure composite (this is paper Fig 5). |
| `FIGURE FILES/FINAL FIGURES 220317/EboV_VHH_Suppl_Figures.ai` | Original supplementary composite (S3 and related). |

## Semi-processed inputs for Figure 5

| Path | What |
| --- | --- |
| `TEMPORARYDATA/CS3_G10_Hidde1/` | G10 LLSM: `M488.mat`, `M560.mat`, `data488master.mat`, Ex10–Ex14 (Ex11 used in 5D). |
| `TEMPORARYDATA/CS_G84/` | G84 LLSM: same tables, Ex05–Ex09 (Ex07 used in 5C). |
| `FIGURE FILES/CFDs/` | Source CFD `.fig` files before the 2022 restyling. |
| `FIGURE FILES/MIPs_Cells/` | MIPs, cell-shape `.fig` files, AP2 intensity distribution. |
| `FIGURE FILES/CELLS/` | Stacked histograms for precipitation and AP2 association (S3C). |
| `FIGURE FILES/CFDs.mat`, `NumbersFromInVitroFigs_220206.mat` | Extracted numeric dumps. |

## Semi-processed inputs for S3 (in vitro / SDCM)

| Path | What |
| --- | --- |
| `FIGURE FILES/InVitro FIGURES/` | `SglVir_GFP.fig`, `VHHpVir_AF647.fig`, `Dist_RFP_clusters.fig`, `FIGURES_EboV.m`. |
| `VSV Ebola (Jason)/SDCM/` | Spinning-disk experiments (44 nM / 440 nM / SM calibration). |
| `VSV Ebola (Jason)/FIGURES/` | Per-condition in-vitro figure folders. |

Gap: `VSV_binding.mat` is **not** in this E:\ tree. `FIGURES_EboV.m` originally loaded it from Dropbox (`...\VSV Ebola (Jason)\SDCM\VSV_binding.mat`). Current S3A/B work can start from the `.fig` files; regenerating from detections needs that `.mat` (or a rebuild from the SDCM `detection_v2.mat` files).

## Processing scripts (copied into `src/matlab/legacy/`)

| Original | Role |
| --- | --- |
| `2021 Processing In Vivo/DataProc_211102_VHH_EboV_data_in_cells.m` | LLSM detection concat → `M488` / cell association. |
| `2021 Processing In Vitro/.../VSV_Ebola_InVitro_DataProcessing211002.m` | SDCM pooling / bleed-through / occupancy. |
| `FIGURE FILES/FINAL FIGURES 220317/FinalFigs_EboV_VHH_CellSHP_and_ecdfs.m` | Interactive cell-shape styling + CFD restyle used for the 2022 freeze. |
| `FIGURE FILES/InVitro FIGURES/FIGURES_EboV.m` | In-vitro paper panels. |
| `FIGURE FILES/Figures_Cells_VHH_EboV.m` | Cellular overview figures. |
| `VSV Ebola (Jason)/SDCM/VSV_VHH_BindingDist_1903xx.m` | Binding-distribution plots. |

## Raw / near-raw LLSM (not needed to restyle current panels)

| Path | What |
| --- | --- |
| `Rasmus (from TKLAB)/LLSM_VSV_VHH/` | Copy of the lattice dataset. |
| `20190415_p5_p55_sCMOS_RH_SM_MeGFP_AF647/` | Another LLSM/SM copy at the project root. |
| `SUM51 for hidde summary/`, `SUM52/` | Per-cell condition folders; G10 analysis in SUM51 was noted as damaged, which is why `TEMPORARYDATA/CS3_G10_Hidde1` was used instead. |

## Representative cells used in the 2022 freeze

- 5A: Ex01, no MOI
- 5B: Ex15, no VHH
- 5C: Ex07, G84
- 5D: Ex11, G10
