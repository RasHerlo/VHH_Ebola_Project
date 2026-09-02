# VHH Ebola figures

Workspace for redrawing **Figure 5** and its associated supplementary panels (S3) from *Anti-ebolavirus single domain antibody fragments block infection at distinct steps during virus infection*.

Data stay on disk at `E:\Projects Academia Ongoing\Ebola VHH (Hidde Lab)`. This repo holds code, a path map, and exported figures. The 2022 figure-file folders are labeled Fig 1A–J; the paper draft numbers that same figure as Figure 5. Code and docs here use the paper numbering.

## MATLAB session

From the repository root:

```matlab
setup
ebola_vhh.list_panels          % story + every panel/subpart
ebola_vhh.fig5.E.cfd           % one subpart
ebola_vhh.run_panel('fig5','C') % all subparts of one panel, then compose
ebola_vhh.fig5.assemble        % preview composite from exports already on disk
ebola_vhh.s3.assemble('RunPanels', true)
```

Edit `config/paths.json` if the academic folder moves. The committed template is `config/paths.example.json`.

## Layout

```
config/                 paths + figure_manifest (claims, layout, calls)
docs/                   inventory, panel map, meeting notes (see docs/meetings/)
src/matlab/+ebola_vhh/  current figure code (one function per subpart)
src/matlab/legacy/      unmodified 2021–2022 scripts
src/python/ebola_vhh/   path/manifest loaders for a later rewrite
output/                 generated PDF/PNG (gitignored except this README)
```

Work on a single subpart by opening its function, for example `src/matlab/+ebola_vhh/+fig5/+E/cfd.m`. Those files currently re-export the 2022 `.fig` / TIFF sources. Replace the `ebola_vhh.render(...)` call with new plotting when that panel is being redrawn. `ebola_vhh.style` is the shared color/font set.

`ebola_vhh.fig5.assemble` and `ebola_vhh.s3.assemble` build a labeled preview PDF/PNG from the panel exports so the full story can be checked without Illustrator. The original Illustrator composites remain the reference on E:\ (`EboV_VHH_Main_Fig.ai`, `EboV_VHH_Suppl_Figures.ai`).

## Figure 5 story (short)

G10 and G84 both neutralize VSV-EBOV, but not at the same step. The 19 Aug 2026 discussion ([docs/meetings/260819.md](docs/meetings/260819.md)) frames the redraw around three points: **G10 occupancy reduces binding**; **G84 occupancy enriches internalized clusters**; **AP2-association does not differ**, so the effects are not a pathway shift. Panels **A–D** are the representative cells; **E–J** are the occupancy CFDs.

S3 is the glass-coverslip calibration (MeGFP intensity, VHHs per virion at 44 nM) plus a localization summary (precipitation / AP2 association by condition).
