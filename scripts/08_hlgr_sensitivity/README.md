# HLGR sensitivity module

This downstream-only module compares the frozen HLGR genomic proxy calls between
the reproductive and bacteraemia groups, then repeats the comparison after
removing all 13 ST6 genomes. It reads the canonical metadata and integrated
MLST/AMR proxy table; it does not rerun or rewrite either upstream pipeline.

From the repository root:

```bash
python3 -m pip install -r scripts/downstream_requirements.txt
python3 scripts/08_hlgr_sensitivity/hlgr_sensitivity.py
```

The full module run creates both an ST × proxy prevalence dot plot and a grouped
ST bar chart (reproductive versus bacteraemia prevalence, with source-specific
genome denominators). To regenerate only those charts from the existing
downstream ST table, run:

```bash
python3 scripts/08_hlgr_sensitivity/plot_st_proxy_dotplot.py --overwrite
python3 scripts/08_hlgr_sensitivity/plot_st_proxy_grouped_bar.py --overwrite
```

Use `--overwrite` to replace this module's deterministic outputs. Tables are
written to `results/tables/hlgr_sensitivity/` and plots to
`results/figures/hlgr_sensitivity/`. The odds ratio is reproductive versus
bacteraemia, the prevalence difference is reproductive minus bacteraemia, and
the reported two-sided test is Fisher's exact test. The 95% interval is a Woolf
log-odds interval (with a documented 0.5 correction only if a cell is zero).
