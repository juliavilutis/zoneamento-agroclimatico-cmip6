# ═══ CÉLULA 5: ensemble CMIP6 (6 modelos x 3 cenários x 2 períodos) aplicado ao Brasil ═══
import os, glob, json, zipfile, requests, numpy as np, pandas as pd, xarray as xr, geopandas as gpd
from rasterio import features
from rasterio.transform import from_origin

DIR = "/content/drive/MyDrive/era5_brasil"
V4  = f"{DIR}/v4"
ENS = f"{DIR}/cmip6/ensemble"
LIMIAR = 0.5
VARS = ["tmean", "tmin", "tmax", "precip", "frio"]
TERM = [0, 1, 2, 4]
MODELOS  = ["miroc6", "gfdl_esm4", "ipsl_cm6a_lr", "mpi_esm1_2_lr", "mri_esm2_0", "ukesm1_0_ll"]
CENARIOS = {"ssp2_4_5": "SSP2-4.5", "ssp3_7_0": "SSP3-7.0", "ssp5_8_5": "SSP5-8.5"}
PERIODOS = {"2041-2060": (2041, 2060), "2071-2090": (2071, 2090)}
PAR = json.load(open(f"{V4}/parametros_v4.json"))

# ── Clima atual do Brasil, máscara, área e uso do solo (igual à célula 3) ───────
ds = xr.open_dataset(f"{DIR}/climatologia_completa_brasil.nc", engine="netcdf4", decode_timedelta=False)
lat, lon = ds.lat.values, ds.lon.values
tmean_m, tmin_m, tmax_m = ds.tmean.values, ds.tmin.values, ds.tmax.values
precip_anual = ds.precip.sum("month").values

if not glob.glob("/tmp/estados_shp/*.shp"):
    r = requests.get("https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_1_states_provinces.zip", stream=True, timeout=120)
    with open("/tmp/estados.zip", "wb") as f:
        for parte in r.iter_content(1 << 20): f.write(parte)
    with zipfile.ZipFile("/tmp/estados.zip") as z: z.extractall("/tmp/estados_shp/")
estados = gpd.read_file(glob.glob("/tmp/estados_shp/*.shp")[0]); estados = estados[estados["iso_a2"] == "BR"]
br = features.rasterize([(g, 1) for g in estados.geometry], out_shape=(lat.size, lon.size),
                        transform=from_origin(lon[0] - 0.05, lat[0] + 0.05, 0.1, 0.1)).astype(bool)
br &= ~np.isnan(tmean_m[0]) & ~np.isnan(tmin_m[0])
area_ha = np.repeat(((111.32 * 0.1) ** 2 * np.cos(np.deg2rad(lat)) * 100)[:, None], lon.size, axis=1)
agri = xr.open_dataset(f"{DIR}/mapbiomas_mask_era5.nc", engine="netcdf4")["mascara"].values == 1
mha = lambda cond, uso=True: area_ha[br & cond & (agri if uso else True)].sum() / 1e6

H    = np.cos(2 * np.pi * (np.arange(24) + 0.5) / 24)
DIAS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
import json
from scipy import stats as sp
SIGMA = json.load(open(f"{V4}/frio_sigma_v4.json"))["sigma"]      # variabilidade entre dias, calibrada na célula 8
def horas_frio(tmin, tmax, sigma=SIGMA, sup=7.2, inf=0.0):
    med, amp = (tmin + tmax) / 2, (tmax - tmin) / 2
    total = np.zeros(tmin.shape[1:])
    for m in range(12):
        for h in H:
            temp = med[m] + amp[m] * h
            total += DIAS[m] * (sp.norm.cdf((sup - temp) / sigma) - sp.norm.cdf((inf - temp) / sigma))
    return total

def trap(v, a, b, c, d):
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(v <= a, 0., np.where(v < b, (v - a) / (b - a),
               np.where(v <= c, 1., np.where(v < d, (d - v) / (d - c), 0.)))).astype(float)

def avaliar(par, tm, tn, tx, pr):
    X = {"tmean": tm.mean(0), "tmin": tn.min(0), "tmax": tx.max(0), "precip": pr, "frio": horas_frio(tn, tx)}
    s = [trap(X[k], *par[k]) for k in VARS]; media = sum(s) / 5
    return (np.where(np.minimum.reduce(s) > 0, media, 0.) > LIMIAR,          # apta: todos os critérios
            np.minimum.reduce([s[q] for q in TERM]) > 0)                     # térmica: regime térmico no envelope

# ── Leitura de um arquivo do ensemble: climatologia mensal na grade do Brasil ───
def clim(modelo, var, exp, anos):
    d = xr.open_dataset(f"{ENS}/{modelo}/{var}_{exp}.nc")[var]
    d = d.rename({k: v for k, v in {"latitude": "lat", "longitude": "lon"}.items() if k in d.dims})
    if float(d.lon.max()) > 180: d = d.assign_coords(lon=((d.lon + 180) % 360) - 180).sortby("lon")
    d = d.sortby("lat")
    ano = d.time.dt.year
    d = d.sel(time=(ano >= anos[0]) & (ano <= anos[1])).groupby("time.month").mean("time")
    return d.interp(lat=lat, lon=lon, method="linear", kwargs={"fill_value": "extrapolate"}).values

# ── Aptidão atual (referência) ──────────────────────────────────────────────────
ATUAL = {c: avaliar(p, tmean_m, tmin_m, tmax_m, precip_anual) for c, p in PAR.items()}

# ── Laço principal: modelo x cenário x período ──────────────────────────────────
linhas, vies, conc = [], [], {}
for modelo in MODELOS:
    hist = {v: clim(modelo, v, "historical", (1991, 2014)) for v in ("tas", "tasmin", "tasmax", "pr")}
    vies.append({"modelo": modelo, "vies_tmean_Brasil_C": np.nanmean((hist["tas"].mean(0) - 273.15 - tmean_m.mean(0))[br]),
                 "vies_tmean_Sul_C": np.nanmean((hist["tas"].mean(0) - 273.15 - tmean_m.mean(0))[br & (lat[:, None] < -20)])})
    for exp, nome_exp in CENARIOS.items():
        for nome_per, anos in PERIODOS.items():
            d = {v: clim(modelo, v, exp, anos) - hist[v] for v in ("tas", "tasmin", "tasmax")}
            razao = np.clip(clim(modelo, "pr", exp, anos).sum(0) / np.maximum(hist["pr"].sum(0), 1e-12), 0.3, 1.5)
            for cultura, par in PAR.items():
                apta, term = avaliar(par, tmean_m + d["tas"], tmin_m + d["tasmin"], tmax_m + d["tasmax"], precip_anual * razao)
                linhas.append({"cultura": cultura, "cenario": nome_exp, "periodo": nome_per, "modelo": modelo,
                               "dT_C": np.nanmean(d["tas"].mean(0)[br]), "dTmin_C": np.nanmean(d["tasmin"].mean(0)[br]),
                               "dTmax_C": np.nanmean(d["tasmax"].mean(0)[br]), "chuva_%": np.nanmean(razao[br]) * 100 - 100,
                               "apta_MapBiomas_Mha": mha(apta), "termica_bruta_Mha": mha(term, uso=False),
                               "termica_MapBiomas_Mha": mha(term)})
                k = f"{cultura}|{nome_exp}|{nome_per}"
                conc[k] = conc.get(k, np.zeros(br.shape)) + term
    print("processado:", modelo)

tab = pd.DataFrame(linhas).round(2); tab.to_csv(f"{V4}/ensemble_modelos_v4.csv", index=False)
tab_vies = pd.DataFrame(vies).round(2); tab_vies.to_csv(f"{V4}/ensemble_vies_v4.csv", index=False)

# ── Resumo do ensemble e mapas de concordância entre modelos ────────────────────
g = tab.groupby(["cultura", "cenario", "periodo"], sort=False)
resumo = g.agg(dT_medio=("dT_C", "mean"), dT_min=("dT_C", "min"), dT_max=("dT_C", "max"),
               apta_max=("apta_MapBiomas_Mha", "max"),
               termica_media=("termica_MapBiomas_Mha", "mean"), termica_min=("termica_MapBiomas_Mha", "min"),
               termica_max=("termica_MapBiomas_Mha", "max"),
               modelos_com_area=("termica_MapBiomas_Mha", lambda s: int((s > 0.05).sum()))).reset_index()
resumo.insert(3, "termica_atual", resumo.cultura.map({c: mha(ATUAL[c][1]) for c in PAR}))
resumo["area_4_de_6_modelos"] = [mha(conc[f"{r.cultura}|{r.cenario}|{r.periodo}"] >= 4) for r in resumo.itertuples()]
resumo = resumo.round(2); resumo.to_csv(f"{V4}/ensemble_resumo_v4.csv", index=False)

mapas = {}
for c in PAR:
    a = ATUAL[c][1].astype(float); a[~br] = np.nan; mapas[f"{c}|atual"] = (("lat", "lon"), a)
for k, v in conc.items():
    v = v.astype(float); v[~br] = np.nan; mapas[k] = (("lat", "lon"), v)
xr.Dataset(mapas, coords={"lat": lat, "lon": lon}).to_netcdf(f"{V4}/ensemble_concordancia_v4.nc", mode="w")

pd.set_option("display.width", 250)
print("\nVIÉS DE CADA MODELO (histórico 1991-2014 menos ERA5-Land), em °C\n" + tab_vies.to_string(index=False))
print("\nAQUECIMENTO MÉDIO NO BRASIL POR MODELO (°C)\n" +
      tab[tab.cultura == list(PAR)[0]].pivot_table(index="modelo", columns=["cenario", "periodo"], values="dT_C", sort=False).round(2).to_string())
print("\nRESUMO DO ENSEMBLE (áreas em Mha, com filtro MapBiomas)\n" + resumo.to_string(index=False))
print("\nÁREA TERMICAMENTE APTA POR MODELO (Mha, com MapBiomas)")
for c in PAR:
    print(f"\n{c}\n" + tab[tab.cultura == c].pivot_table(index="modelo", columns=["cenario", "periodo"], values="termica_MapBiomas_Mha", sort=False).to_string())
print("\nCélula 5 concluída. Arquivos em", V4)
