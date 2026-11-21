# Climatic limits to the introduction of hops, pistachio, hazelnut and almond in Brazil

Code and result tables for the article:

> Bronzeri, J.V., Nassif, D.S.P., Gazzola, J. Climatic limits to the introduction of hops, pistachio, hazelnut and almond in Brazil: validated climate envelopes, limiting factors and CMIP6 projections. Submitted to *European Journal of Agronomy*.

The pipeline calibrates climate envelopes for four crops from the climate of their global harvested area (CROPGRIDS v1.08), validates them outside Brazil, calibrates the estimate of annual chilling hours against 231 INMET automatic weather stations, and applies the envelopes to Brazil under the reference climate (1991-2020) and under a six-model CMIP6 ensemble (SSP2-4.5, SSP3-7.0 and SSP5-8.5; 2041-2060 and 2071-2090).

## Repository structure

```
scripts/    13 scripts, run in numerical order
results/    result tables (CSV/JSON) used in the article tables and text
results/figures/   article figures (PNG)
```

Code comments and printed messages are in Portuguese; variable names follow the article (tmean, tmin, tmax, precip, frio = chilling hours).

## How to run

The scripts were written and run as Google Colab cells, reading and writing to a Google Drive folder (`/content/drive/MyDrive/era5_brasil`). To reproduce:

1. Open a Colab notebook, mount Google Drive (`from google.colab import drive; drive.mount("/content/drive")`) and install the dependencies (`pip install -r requirements.txt`).
2. Create an account at the Copernicus Climate Data Store (CDS), accept the ERA5-Land and CMIP6 licences, and replace `SEU-TOKEN-AQUI` with your personal token in scripts 01, 05 and 11.
3. Download the MapBiomas Collection 9 land cover map for 2023 (GeoTIFF, 30 m) from https://brasil.mapbiomas.org and set its path in script 04.
4. Paste and run each script, in order, as a Colab cell.

| Step | Script | What it does | Main outputs |
|---|---|---|---|
| 1 | `01_era5_land_download.py` | ERA5-Land monthly means 1991-2020 (2 m temperature, precipitation, potential evaporation) | `brasil_*_1991_2020.nc` |
| 2 | `02_worldclim_brazil.py` | WorldClim v2.1 (10 arc-min) monthly Tmin and Tmax, clipped to Brazil | `worldclim_*_brasil.nc` |
| 3 | `03_baseline_climatology.py` | Reference climatology on the 0.1° ERA5-Land grid | `climatologia_completa_brasil.nc` |
| 4 | `04_mapbiomas_mask.py` | Majority MapBiomas class per 0.1° cell and agricultural-use mask | `mapbiomas_mask_era5.nc` |
| 5 | `05_miroc6_download.py` | MIROC6 tas and pr (historical, SSP2-4.5, SSP5-8.5) for the single-model analysis | `cmip6/*.nc` |
| 6 | `06_inmet_chill_validation.py` | Observed chilling hours at INMET stations (2008-2020) and comparison with monthly-data methods | `frio_inmet_anual_v4.csv`, `validacao_frio_inmet_v4.csv`, `validacao_frio_resumo_v4.csv` |
| 7 | `07_chill_sigma_calibration.py` | Calibration of day-to-day variability (σ) in the sinusoidal chill model, with 5-fold cross-validation | `frio_sigma_v4.json`, `frio_sigma_curva_v4.csv`, `frio_sigma_resumo_v4.csv`, Figure 1 |
| 8 | `08_global_data.py` | Global WorldClim variables (with calibrated chilling hours), CROPGRIDS maps, municipality coordinates | `worldclim_vars.npz`, `cropgrids/` |
| 9 | `09_calibration_global_validation.py` | Envelope calibration (area-weighted percentiles), global validation, leave-one-region-out, independent producing regions, sensitivity | `parametros_v4.csv/json`, `percentis_v4.json`, `validacao_global_v4.csv`, `validacao_pontos_v4.csv`, `sensibilidade_global_v4.csv` |
| 10 | `10_brazil_suitability.py` | Application to Brazil, limiting factors, MIROC6 scenarios, sensitivity, hop-producing municipalities | `brasil_v4.csv`, `fatores_limitantes_v4.nc`, `sensibilidade_brasil_v4.csv`, `lupulo_municipios_v4.csv` |
| 11 | `11_cmip6_ensemble_download.py` | CMIP6 ensemble (MIROC6, GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-LR, MRI-ESM2-0, UKESM1-0-LL; tas, tasmin, tasmax, pr). The download queue resumes if the session drops. | `cmip6/ensemble/{model}/{var}_{experiment}.nc` |
| 12 | `12_cmip6_ensemble.py` | Monthly delta method per variable, model bias, warming, thermal suitability by model, scenario and period, model agreement | `ensemble_modelos_v4.csv`, `ensemble_resumo_v4.csv`, `ensemble_vies_v4.csv`, `ensemble_concordancia_v4.nc` |
| 13 | `13_figures.py` | Article figures (`IDIOMA = "en"`) | `figuras/fig*.png` |

Steps 6 and 7 must run before steps 8, 10 and 12, which read the calibrated σ from `frio_sigma_v4.json`.

## Result files (`results/`)

| File | Content | Used in |
|---|---|---|
| `parametros_v4.csv`, `parametros_v4.json` | Envelope parameters [a, b, c, d] by crop and variable | Table 3 |
| `percentis_v4.json` | Area-weighted percentiles used for calibration and sensitivity | Section 2.5 |
| `validacao_global_v4.csv` | Harvested area suitable, global land suitable, AUC, leave-one-region-out | Table 4 |
| `validacao_pontos_v4.csv` | Climate and suitability at 55 independent producing regions | Table 4, Section 3.2 |
| `sensibilidade_global_v4.csv`, `sensibilidade_brasil_v4.csv` | Sensitivity to percentiles, threshold and aggregation rule | Table 8 |
| `brasil_v4.csv` | Suitable and thermally suitable area in Brazil, share of territory excluded by each criterion | Table 5 |
| `fatores_limitantes_v4.nc` | Maps of criteria met and limiting factors in Brazil | Figs. 3 and 4 |
| `lupulo_municipios_v4.csv` | Climate and criteria met at 90 hop-producing municipalities | Section 3.3 |
| `ensemble_modelos_v4.csv`, `ensemble_resumo_v4.csv`, `ensemble_vies_v4.csv` | Results by CMIP6 model, ensemble summary, model bias | Tables 6 and 7, Fig. 6 |
| `ensemble_concordancia_v4.nc` | Number of models placing each cell within the thermal envelope | Fig. 5 |
| `frio_inmet_anual_v4.csv` | Observed annual chilling hours by INMET station and year | Section 2.3 |
| `validacao_frio_inmet_v4.csv`, `validacao_frio_resumo_v4.csv` | Observed and estimated chilling hours by station; error statistics | Table 2 |
| `frio_sigma_v4.json`, `frio_sigma_curva_v4.csv`, `frio_sigma_resumo_v4.csv` | Calibrated σ, error curve and cross-validation | Table 2, Fig. 1 |

The two NetCDF files were saved with zlib compression after the analysis; their values are identical to the files written by the scripts.

## Correspondence between figures and files

| Article | File |
|---|---|
| Fig. 1 | `fig6_validacao_frio.png` |
| Fig. 2 | `fig5_indice_global.png` |
| Fig. 3 | `fig1_criterios_atendidos.png` |
| Fig. 4 | `fig2_fatores_limitantes.png` |
| Fig. 5 | `fig3_concordancia_modelos.png` |
| Fig. 6 | `fig4_area_por_cenario.png` |

## Notes

- Scripts 01 to 05 consolidate the preprocessing cells of the original exploratory notebook into single scripts; the data requests and processing steps are the same. Script 03 was checked against the original `climatologia_completa_brasil.nc` (identical temperature and precipitation fields).
- Monthly precipitation totals were obtained from ERA5-Land mean daily rates multiplied by 30 days, as in the original processing.
- Tmin and Tmax come from WorldClim v2.1 (1970-2000), whereas mean temperature and precipitation come from ERA5-Land (1991-2020). This mismatch is discussed in the article.
- The INMET hourly data (step 6) are downloaded from https://portal.inmet.gov.br/dadoshistoricos.

## Data sources

All input data are public:

- ERA5-Land monthly means (Muñoz-Sabater et al., 2021) and CMIP6 projections, Copernicus Climate Data Store
- WorldClim v2.1 (Fick and Hijmans, 2017)
- CROPGRIDS v1.08 (Tang et al., 2024), https://doi.org/10.6084/m9.figshare.22491997
- MapBiomas Collection 9 (Souza et al., 2020)
- INMET automatic weather stations
- Natural Earth admin-1 boundaries
- Municipality coordinates: https://github.com/kelvins/municipios-brasileiros

## Licence

Code: MIT (see `LICENSE`). Input datasets keep their original licences.

## How to cite

See `CITATION.cff` or the Zenodo record of this repository.
