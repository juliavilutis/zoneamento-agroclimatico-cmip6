# ═══ CÉLULA 3 de 3: aplicar ao Brasil (atual e cenários), fatores limitantes e sensibilidade ═══
import os, glob, json, zipfile, requests, unicodedata, numpy as np, pandas as pd, xarray as xr, geopandas as gpd
from rasterio import features
from rasterio.transform import from_origin

DIR = "/content/drive/MyDrive/era5_brasil"
V4  = f"{DIR}/v4"
CMIP6 = f"{DIR}/cmip6"
LIMIAR = 0.5
VARS = ["tmean", "tmin", "tmax", "precip", "frio"]
TERM = [0, 1, 2, 4]                               # critérios térmicos (todos menos a chuva)
PAR  = json.load(open(f"{V4}/parametros_v4.json"))
PERC = json.load(open(f"{V4}/percentis_v4.json"))

# ── Clima do Brasil, máscara, área por pixel e uso do solo ──────────────────────
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
BRASIL = mha(np.ones_like(br), uso=False)
print(f"Brasil: {BRASIL:.0f} Mha | uso agrícola compatível: {mha(np.ones_like(br)) / BRASIL * 100:.1f}%")

# ── Funções (as mesmas da célula 2) ─────────────────────────────────────────────
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

def clima(tm, tn, tx, pr):
    return {"tmean": tm.mean(0), "tmin": tn.min(0), "tmax": tx.max(0), "precip": pr, "frio": horas_frio(tn, tx)}

def avaliar(par, X):
    s = [trap(X[k], *par[k]) for k in VARS]
    media = sum(s) / 5
    return {"notas": s, "media": media,
            "apta": np.where(np.minimum.reduce(s) > 0, media, 0.) > LIMIAR,       # todos os critérios
            "termica": np.minimum.reduce([s[q] for q in TERM]) > 0,               # regime térmico dentro do envelope
            "n_ok": sum((x > 0).astype(int) for x in s)}

# ── Cenários do MIROC6: delta mensal em relação ao histórico do próprio modelo ──
def cmip(nome, var):
    d = xr.open_dataset(f"{CMIP6}/{nome}.nc", engine="netcdf4")[var]
    d = d.assign_coords(lon=d.lon.values - 360).groupby("time.month").mean("time")
    return d.interp(lat=lat, lon=lon, method="linear", kwargs={"fill_value": "extrapolate"}).values

tas_hist = cmip("tas_historical_1991_2014", "tas"); pr_hist = cmip("pr_historical_1991_2014", "pr").sum(0)
CLIMAS = {"atual": clima(tmean_m, tmin_m, tmax_m, precip_anual)}
for key in [f"{s}_{p}" for s in ("ssp245", "ssp585") for p in ("2041_2060", "2071_2090")]:
    delta = cmip(f"tas_{key}", "tas") - tas_hist
    razao = np.clip(cmip(f"pr_{key}", "pr").sum(0) / np.maximum(pr_hist, 1e-12), 0.3, 1.5)
    CLIMAS[key] = clima(tmean_m + delta, tmin_m + delta, tmax_m + delta, precip_anual * razao)
    print(f"{key}: ΔT = {np.nanmean(delta.mean(0)[br]):.2f} °C | chuva = {np.nanmean(razao[br]) * 100 - 100:+.1f}%")

# ── Tabela principal: área apta, área termicamente apta e critérios que excluem ─
linhas, mapas = [], {}
for cultura, par in PAR.items():
    for cen, X in CLIMAS.items():
        r = avaliar(par, X)
        linha = {"cultura": cultura, "cenario": cen,
                 "apta_bruta_Mha": mha(r["apta"], uso=False), "apta_MapBiomas_Mha": mha(r["apta"]),
                 "termica_bruta_Mha": mha(r["termica"], uso=False), "termica_MapBiomas_Mha": mha(r["termica"]),
                 "4_de_5_criterios_Mha": mha(r["n_ok"] == 4, uso=False)}
        for q, k in enumerate(VARS):
            linha[f"zerado_{k}_%"] = mha(r["notas"][q] == 0, uso=False) / BRASIL * 100
        linhas.append(linha)
        if cen == "atual":
            n = r["n_ok"].astype(float); n[~br] = np.nan
            mapas[f"{cultura}_criterios_atendidos"] = (("lat", "lon"), n)
            for q, k in enumerate(VARS):
                z = (r["notas"][q] == 0).astype(float); z[~br] = np.nan
                mapas[f"{cultura}_zerado_{k}"] = (("lat", "lon"), z)
tab = pd.DataFrame(linhas).round(2); tab.to_csv(f"{V4}/brasil_v4.csv", index=False)
xr.Dataset(mapas, coords={"lat": lat, "lon": lon}).to_netcdf(f"{V4}/fatores_limitantes_v4.nc", mode="w")
pd.set_option("display.width", 250)
print("\nBRASIL COM PARÂMETROS CALIBRADOS\n" + tab.to_string(index=False))

# ── Sensibilidade: percentis do envelope e limiar ───────────────────────────────
QS = PERC["percentis"]; linhas = []
for P in [(0.5, 5, 95, 99.5), (1, 5, 95, 99), (2.5, 10, 90, 97.5), (5, 15, 85, 95)]:
    for cultura in PAR:
        par = {k: [PERC["valores"][cultura][k][QS.index(p)] for p in P] for k in VARS}
        par["frio"][2:] = [1e9, 2e9]
        s = [trap(CLIMAS["atual"][k], *par[k]) for k in VARS]; media = sum(s) / 5
        for lim in (0.4, 0.5, 0.6):
            linhas.append({"percentis": str(P), "limiar": lim, "cultura": cultura,
                           "limitante_Mha": mha((np.minimum.reduce(s) > 0) & (media > lim), uso=False),
                           "so_termico_exclui_Mha": mha((np.minimum.reduce([s[q] for q in TERM]) > 0) & (media > lim), uso=False),
                           "media_simples_Mha": mha(media > lim, uso=False)})
tab_s = pd.DataFrame(linhas).round(2); tab_s.to_csv(f"{V4}/sensibilidade_brasil_v4.csv", index=False)
print("\nSENSIBILIDADE (área bruta no Brasil, Mha, limiar 0,5)\n" + tab_s[tab_s.limiar == 0.5].drop(columns="limiar").to_string(index=False))

# ── Lúpulo: os 90 municípios produtores (Aprolúpulo, 2024) e o fotoperíodo ──────
LUPULO = {"AL":["União dos Palmares"],"BA":["Vitória da Conquista"],"DF":["Brasília"],"ES":["Alfredo Chaves","Laranja da Terra","Viana"],
 "GO":["Hidrolândia","Santa Cruz de Goiás","Trindade"],"MT":["Tesouro"],
 "MG":["Alvorada de Minas","Andrelândia","Araguari","Carmo da Cachoeira","Carmo do Cajuru","Coqueiral","Dom Silvério","Gonçalves","Igarapé","Itajubá","Madre de Deus de Minas","Mateus Leme","Monte Santo de Minas","Poços de Caldas","Ponte Nova","Prata","Santa Rita do Sapucaí","Uberaba","Viçosa"],
 "PR":["Araucária","Carambeí","Guarapuava","Imbituva","Iporã","Nova Olímpia","Palotina","Pato Branco","Santa Izabel do Oeste","Umuarama"],
 "PI":["Curimatá"],"RJ":["Petrópolis","Sapucaia","Teresópolis"],
 "RS":["Bento Gonçalves","Bom Jesus","Capão do Leão","Forquetinha","Garibaldi","Maratá","Marau","Roca Sales","São Francisco de Paula","Venâncio Aires"],
 "SC":["Araranguá","Bocaina do Sul","Bom Retiro","Campo Alegre","Itá","Jacinto Machado","Joinville","Paulo Lopes","Pouso Redondo","São Joaquim","Urubici","Videira"],
 "SP":["Aguaí","Analândia","Anhembi","Araçoiaba da Serra","Araraquara","Botucatu","Cajuru","Cedral","Cruzeiro","Cunha","Fartura","Iguape","Itapetininga","Itaporanga","Jacareí","Mococa","Mogi das Cruzes","Pariquera-Açu","Ribeirão Bonito","Saltinho","Santa Rosa de Viterbo","São José do Barreiro","São Roque","Sorocaba","Taguaí"]}
sem_acento = lambda t: unicodedata.normalize("NFKD", t).encode("ascii", "ignore").decode().lower().strip()
mun = pd.read_csv(f"{V4}/municipios.csv").merge(pd.read_csv(f"{V4}/estados.csv", encoding="utf-8-sig")[["codigo_uf", "uf"]], on="codigo_uf")
mun["chave"] = mun.nome.map(sem_acento) + "|" + mun.uf
def horas_luz(lat_graus):
    x = np.clip(np.tan(np.deg2rad(np.abs(lat_graus))) * np.tan(np.deg2rad(23.44)), -1, 1)
    return 24 / np.pi * np.arccos(-x)
r = avaliar(PAR["Lúpulo"], CLIMAS["atual"]); validos = np.where(br); linhas = []
for uf, nomes in LUPULO.items():
    for nome in nomes:
        m = mun[mun.chave == sem_acento(nome) + "|" + uf]
        if len(m) != 1: print("Município não encontrado:", nome, uf); continue
        la, lo = float(m.latitude.iloc[0]), float(m.longitude.iloc[0])
        q = ((lat[validos[0]] - la) ** 2 + (lon[validos[1]] - lo) ** 2).argmin(); i, j = validos[0][q], validos[1][q]
        linhas.append({"municipio": nome, "uf": uf, **{k: round(float(CLIMAS["atual"][k][i, j]), 1) for k in VARS},
                       "criterios_atendidos": int(r["n_ok"][i, j]), "horas_luz_solsticio": round(float(horas_luz(la)), 1),
                       "criterios_zerados": ", ".join(VARS[x] for x in range(5) if r["notas"][x][i, j] == 0)})
tab_l = pd.DataFrame(linhas); tab_l.to_csv(f"{V4}/lupulo_municipios_v4.csv", index=False)
luz_br = np.repeat(horas_luz(lat)[:, None], lon.size, 1); p25 = PERC["fotoperiodo_lupulo_P1_P2.5_P5_P10_P50"][1]
print(f"\nLÚPULO: {len(tab_l)} municípios produtores | critérios atendidos (de 5):", tab_l.criterios_atendidos.value_counts().sort_index().to_dict())
print(f"Horas de luz no solstício: municípios de {tab_l.horas_luz_solsticio.min()} a {tab_l.horas_luz_solsticio.max()} h | Brasil, máximo {luz_br[br].max():.1f} h | "
      f"áreas de lúpulo no mundo, P2,5 = {p25:.1f} h | área do Brasil acima disso: {mha(luz_br >= p25, uso=False):.2f} Mha")
print("\nCélula 3 concluída. Arquivos em", V4)
