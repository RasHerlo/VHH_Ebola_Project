# Panel numbering and claims

The 17 Mar 2022 figure dump (`FIGURE FILES/FINAL FIGURES 220317`) labels this figure as **Fig 1**. The 4 Jan 2023 draft (`Ebolavirus VHHs v2 230104.docx`) labels it **Figure 5**. This repo uses paper numbering.

## Figure 5 — G10 blocks cell binding and alters localization

| Paper | 2022 folder | MATLAB call | Claim |
| --- | --- | --- | --- |
| 5A | Fig 1A Ex01_MOI | `ebola_vhh.fig5.A.mip` / `.cell_shape` / `.ap2_dist` | AP2-RFP clusters define the cell surface; intensities span ~an order of magnitude. |
| 5B | Fig 1B Ex15_NoVHH | `ebola_vhh.fig5.B.mip` / `.cell_shape` | Without VHH, MeGFP-VSV binds the reconstructed surface. |
| 5C | Fig 1C Ex7 G84 | `ebola_vhh.fig5.C.full_view` / `.side_view` / `.int_zoom` | G84 still allows cell binding; zooms: internalized examples. Captions TBD. |
| 5D | Fig 1D Ex11 G10 | `ebola_vhh.fig5.D.full_view` / `.side_view` / `.precip_zoom` | G10 yields more off-cell particles; zooms: precipitation. Captions TBD. |
| 5E | Fig 1E ecdf G10 CellBind | `ebola_vhh.fig5.E.cfd` | More G10 per virion → less cell association. |
| 5F | Fig 1F ecdf G84 CellBind | `ebola_vhh.fig5.F.cfd` | G84 occupancy does not reduce binding the way G10 does. |
| 5G | Fig 1G ecdf G10 AP2ass | `ebola_vhh.fig5.G.cfd` | G10 occupancy increases AP2 association; not different from G84 (not a pathway shift). |
| 5H | Fig 1H ecdf G84 APass | `ebola_vhh.fig5.H.cfd` | G84 occupancy increases AP2 association similarly to G10. |
| 5I | Fig 1I ecdf G10 Int | `ebola_vhh.fig5.I.cfd` | High G10 occupancy: virions less likely to be inside the cell. |
| 5J | Fig 1J ecdf G84 Int | `ebola_vhh.fig5.J.cfd` | High G84 occupancy: virions still bind and are more likely in internalized clusters. |

Open items carried over from `Main&Suppl_Figs_Text.docx`: Kolmogorov–Smirnov tests on the CFDs were never added; several zoom captions are still TBD.

Design discussion (19 Aug 2026): [docs/meetings/260819.md](meetings/260819.md). Intended argument for the redraw: (1) G10 occupancy reduces binding; (2) G84 occupancy enriches internalized clusters; (3) AP2-association does not differ between the two, so this is not a pathway shift. That I–J read is sharper than the 2023 draft caption.

## Figure S3 — calibration and localization summaries

The Results text calls occupancy counts S3A/S3B. The legend file uses A = GFP intensity, B = VHHs per virion, C = localization summary. This repo follows the legend file.

| Paper | Source | MATLAB call | Claim |
| --- | --- | --- | --- |
| S3A | `InVitro FIGURES/SglVir_GFP.fig` | `ebola_vhh.s3.A.vsv_intensity` | MeGFP-VSV on glass is a reasonably homogeneous intensity population. |
| S3B | `InVitro FIGURES/VHHpVir_AF647.fig` | `ebola_vhh.s3.B.vhh_per_virion` | At 44 nM: ~45 G10, ~13 G84, ~1 G68 fragments per virion. |
| S3C | `CELLS/StkHist_Prec.fig`, `StkHist_AP2ass.fig` | `ebola_vhh.s3.C.precip_hist` / `.ap2_hist` | G10/G84 increase off-cell precipitation; cell-bound virus shows higher AP2 association, especially inside. |

S3 is not embedded in the Word draft; the Illustrator/PDF composite is `EboV_VHH_Suppl_Figures.ai/.pdf`.

## Candidate (not in the draft)

`ebola_vhh.candidates.ap2_cluster_dist` — AP2 cluster distribution from the in-vitro figure folder.
