# build_colour_index — build checklist

- [X] Rename `build_index.nf`/`.config` to `build_colour_index.nf`/`.config`
- [X] Write `schema.json`
- [X] Copy `color_mapping.py` into `build_colour_index/bin/`
- [X] Write `modules/colour_mapping.nf`
- [X] Write `modules/ggcat.nf`
- [X] Write `modules/sbwt_build.nf`
- [X] Write `modules/themisto2.nf` (Themisto2 build + export)
- [X] Write `build_colour_index.nf` orchestrator
- [X] Write `build_colour_index.config`
- [ ] Write `README.md`
- [ ] Wire into lsmd `main.nf` and test on v_tarriae
- [ ] Run `nextflow run main.nf -profile test` and see what actually happens
