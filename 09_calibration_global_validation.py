# ═══ CÉLULA 2 de 3: calibrar parâmetros com a área colhida global e validar ═══
import os, glob, json, zipfile, requests, numpy as np, pandas as pd, xarray as xr, geopandas as gpd
from rasterio import features
from rasterio.transform import Affine

DIR = "/content/drive/MyDrive/era5_brasil"
V4  = f"{DIR}/v4"
LIMIAR = 0.5
PERCENTIS = (2.5, 10, 90, 97.5)                 # a, b, c, d do trapézio
VARS = ["tmean", "tmin", "tmax", "precip", "frio"]
NOMES = {"hop": "Lúpulo", "pistachio": "Pistache", "hazelnut": "Avelã", "almond": "Amêndoa"}

WC = dict(np.load(f"{V4}/worldclim_vars.npz"))  # grade de 10 min, latitude decrescente

# ── Grade do CROPGRIDS (0,05 grau, latitude crescente) e máscara do Brasil ──────
clat = -89.975 + 0.05 * np.arange(3600); clon = -179.975 + 0.05 * np.arange(7200)
LAT = np.repeat(clat[:, None], 7200, 1); LON = np.repeat(clon[None, :], 3600, 0)
iw = np.clip(((90 - clat) * 6).astype(int), 0, 1079); jw = np.clip(((clon + 180) * 6).astype(int), 0, 2159)

if not glob.glob("/tmp/estados_shp/*.shp"):
    r = requests.get("https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_1_states_provinces.zip", stream=True, timeout=120)
    with open("/tmp/estados.zip", "wb") as f:
        for parte in r.iter_content(1 << 20): f.write(parte)
    with zipfile.ZipFile("/tmp/estados.zip") as z: z.extractall("/tmp/estados_shp/")
estados = gpd.read_file(glob.glob("/tmp/estados_shp/*.shp")[0]); estados = estados[estados["iso_a2"] == "BR"]
brasil_cg = features.rasterize([(g, 1) for g in estados.geometry], out_shape=(3600, 7200),
                               transform=Affine(0.05, 0, -180, 0, 0.05, -90)).astype(bool)

# ── Funções ─────────────────────────────────────────────────────────────────────
def percentil_pond(x, w, q):
    o = np.argsort(x); x, w = x[o], w[o]
    return np.interp(np.array(q) / 100, np.cumsum(w) / w.sum(), x)

def trap(v, a, b, c, d):
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(v <= a, 0., np.where(v < b, (v - a) / (b - a),
               np.where(v <= c, 1., np.where(v < d, (d - v) / (d - c), 0.)))).astype(float)

def calibrar(X, w, P=PERCENTIS):
    par = {}
    for k in VARS:
        a, b, c, d = percentil_pond(X[k], w, P)
        if k == "frio": c, d = 1e9, 2e9          # frio: só limite inferior
        par[k] = [float(a), float(b), float(c), float(d)]
    return par

def indices(par, X):
    s = [trap(X[k], *par[k]) for k in VARS]
    media = sum(s) / 5
    limitante = np.where(np.minimum.reduce(s) > 0, media, 0.)
    return s, media, limitante

def amostra(c):
    """Clima de cada pixel com área colhida (fora do Brasil), e a área como peso."""
    h = xr.open_dataset(f"{V4}/cropgrids/CROPGRIDSv1.08_{c}.nc")["harvarea"].values
    sel = (h > 0) & ~brasil_cg
    if c == "hop":                               # na Etiópia, o "hop" da FAO é gesho (Rhamnus prinoides)
        sel &= ~((LAT >= 3) & (LAT < 15) & (LON >= 33) & (LON < 48))
    yi, xi = np.where(sel)
    X = {k: WC[k][iw[yi], jw[xi]] for k in VARS}
    ok = ~np.isnan(X["tmean"])
    return {k: X[k][ok] for k in VARS}, h[sel][ok], LON[sel][ok], LAT[sel][ok]

# ── Calibração ──────────────────────────────────────────────────────────────────
AM  = {c: amostra(c) for c in NOMES}
PAR = {c: calibrar(*AM[c][:2]) for c in NOMES}
json.dump({NOMES[c]: PAR[c] for c in NOMES}, open(f"{V4}/parametros_v4.json", "w"), ensure_ascii=False, indent=1)
tab_par = pd.DataFrame([{"cultura": NOMES[c], "variavel": k, "a": PAR[c][k][0], "b": PAR[c][k][1],
                         "c": None if k == "frio" else PAR[c][k][2], "d": None if k == "frio" else PAR[c][k][3]}
                        for c in NOMES for k in VARS]).round(1)
tab_par.to_csv(f"{V4}/parametros_v4.csv", index=False)
print("PARÂMETROS CALIBRADOS (percentis", PERCENTIS, "ponderados pela área colhida)\n" + tab_par.to_string(index=False))

# ── Validação global: sensibilidade, seletividade, AUC e teste por região ───────
glat = 90 - (np.arange(1080) + 0.5) / 6
peso_terra = np.cos(np.deg2rad(glat))[:, None] * np.ones((1, 2160))
terra = ~np.isnan(WC["tmean"]) & (glat[:, None] > -60)
REGIOES = {"Américas": lambda lo: lo < -30, "Europa/África/Or. Médio": lambda lo: (lo >= -30) & (lo < 60), "Ásia/Oceania": lambda lo: lo >= 60}
rng = np.random.default_rng(1); linhas = []
for c in NOMES:
    X, w, lo, la = AM[c]
    _, media, lim = indices(PAR[c], X)
    _, _, lim_global = indices(PAR[c], WC)
    fundo = rng.choice(lim_global[terra], 5000, p=peso_terra[terra] / peso_terra[terra].sum())
    pres  = rng.choice(lim, 20000, p=w / w.sum())
    auc = (pres[:, None] > fundo[None, :]).mean() + 0.5 * (pres[:, None] == fundo[None, :]).mean()
    linha = {"cultura": NOMES[c], "area_colhida_mil_ha": w.sum() / 1e3, "pixels": len(w),
             "apta_limitante_%": np.average(lim > LIMIAR, weights=w) * 100,
             "apta_media_%": np.average(media > LIMIAR, weights=w) * 100,
             "terra_global_apta_%": (peso_terra[terra] * (lim_global[terra] > LIMIAR)).sum() / peso_terra[terra].sum() * 100,
             "AUC": auc}
    for nome, f in REGIOES.items():
        teste = f(lo)
        if w[teste].sum() / w.sum() < 0.02: linha[f"fora: {nome}"] = np.nan; continue
        par_sem = calibrar({k: X[k][~teste] for k in VARS}, w[~teste])
        _, _, l2 = indices(par_sem, {k: X[k][teste] for k in VARS})
        linha[f"fora: {nome}"] = np.average(l2 > LIMIAR, weights=w[teste]) * 100
    linhas.append(linha)
tab_val = pd.DataFrame(linhas).round(2); tab_val.to_csv(f"{V4}/validacao_global_v4.csv", index=False)
print("\nVALIDAÇÃO GLOBAL ('fora: região' = calibrado sem a região e testado nela)\n" + tab_val.to_string(index=False))

# ── Validação independente: 55 regiões produtoras conhecidas ────────────────────
PONTOS = {
 "hop": [("Hallertau (Alemanha)",48.60,11.63),("Tettnang (Alemanha)",47.67,9.59),("Zatec (Rep. Tcheca)",50.33,13.55),("Zalec (Eslovênia)",46.25,15.16),("Lublin (Polônia)",51.25,22.57),("Kent (Reino Unido)",51.20,0.75),("Yakima (EUA)",46.38,-120.31),("Willamette (EUA)",45.14,-122.86),("Treasure Valley (EUA)",43.66,-116.69),("Motueka (Nova Zelândia)",-41.12,173.01),("Bushy Park (Austrália)",-42.71,146.90),("George (África do Sul)",-33.96,22.46),("El Bolsón (Argentina)",-41.96,-71.53)],
 "pistachio": [("Rafsanjan (Irã)",30.41,55.99),("Kerman (Irã)",30.28,57.08),("Yazd (Irã)",31.90,54.37),("Madera (EUA)",36.96,-120.06),("Bakersfield (EUA)",35.37,-119.02),("Gaziantep (Turquia)",37.07,37.38),("Sanliurfa (Turquia)",37.16,38.79),("Siirt (Turquia)",37.93,41.94),("Aleppo (Síria)",36.20,37.16),("Bronte (Itália)",37.79,14.83),("Egina (Grécia)",37.75,23.43),("Ciudad Real (Espanha)",38.99,-3.93),("Sfax (Tunísia)",34.74,10.76)],
 "hazelnut": [("Ordu (Turquia)",40.98,37.88),("Giresun (Turquia)",40.91,38.39),("Trabzon (Turquia)",41.00,39.72),("Sakarya (Turquia)",40.78,30.40),("Viterbo (Itália)",42.42,12.10),("Alba (Itália)",44.70,8.03),("Avellino (Itália)",40.91,14.79),("McMinnville (EUA)",45.21,-123.19),("Reus (Espanha)",41.15,1.11),("Zaqatala (Azerbaijão)",41.63,46.64),("Zugdidi (Geórgia)",42.51,41.87),("Talca (Chile)",-35.43,-71.66),("Temuco (Chile)",-38.74,-72.59)],
 "almond": [("Fresno (EUA)",36.75,-119.77),("Modesto (EUA)",37.64,-120.99),("Chico (EUA)",39.73,-121.84),("Bakersfield (EUA)",35.37,-119.02),("Albacete (Espanha)",38.99,-1.86),("Murcia (Espanha)",37.99,-1.13),("Lleida (Espanha)",41.62,0.62),("Granada (Espanha)",37.18,-3.60),("Mildura (Austrália)",-34.21,142.14),("Renmark (Austrália)",-34.17,140.75),("Beja (Portugal)",38.01,-7.86),("Mirandela (Portugal)",41.49,-7.18),("Meknes (Marrocos)",33.89,-5.55),("Bari (Itália)",41.11,16.87),("Agrigento (Itália)",37.31,13.58),("Sfax (Tunísia)",34.74,10.76)],
}
linhas = []
for c, pts in PONTOS.items():
    for nome, la, lo in pts:
        i, j = int((90 - la) * 6), int((lo + 180) * 6)
        x = {k: np.array([WC[k][i, j]]) for k in VARS}
        s, media, lim = indices(PAR[c], x)
        linhas.append({"cultura": NOMES[c], "regiao": nome, **{k: round(float(x[k][0]), 1) for k in VARS},
                       "indice": round(float(lim[0]), 2), "apta": bool(lim[0] > LIMIAR),
                       "criterios_zerados": ", ".join(VARS[q] for q in range(5) if s[q][0] == 0)})
tab_pts = pd.DataFrame(linhas); tab_pts.to_csv(f"{V4}/validacao_pontos_v4.csv", index=False)
print("\nREGIÕES PRODUTORAS INDEPENDENTES: % classificadas aptas")
print((tab_pts.groupby("cultura")["apta"].mean() * 100).round(0).to_string())
print("\nNão classificadas como aptas:\n" + tab_pts[~tab_pts.apta][["cultura", "regiao", "indice", "criterios_zerados"]].to_string(index=False))
# ── Percentis para a análise de sensibilidade e fotoperíodo do lúpulo ───────────
QS = [0.5, 1, 2.5, 5, 10, 15, 85, 90, 95, 97.5, 99, 99.5]
perc = {NOMES[c]: {k: [float(v) for v in percentil_pond(AM[c][0][k], AM[c][1], QS)] for k in VARS} for c in NOMES}
def horas_luz(lat_graus):                        # duração do dia no solstício de verão
    x = np.clip(np.tan(np.deg2rad(np.abs(lat_graus))) * np.tan(np.deg2rad(23.44)), -1, 1)
    return 24 / np.pi * np.arccos(-x)
luz = percentil_pond(horas_luz(AM["hop"][3]), AM["hop"][1], [1, 2.5, 5, 10, 50])
json.dump({"percentis": QS, "valores": perc, "fotoperiodo_lupulo_P1_P2.5_P5_P10_P50": [float(v) for v in luz]},
          open(f"{V4}/percentis_v4.json", "w"), ensure_ascii=False, indent=1)
print("\nFotoperíodo nas áreas de lúpulo (h de luz no solstício): P1=%.1f  P2,5=%.1f  P5=%.1f  P10=%.1f  mediana=%.1f" % tuple(luz))

linhas = []
for P in [(0.5, 5, 95, 99.5), (1, 5, 95, 99), (2.5, 10, 90, 97.5), (5, 15, 85, 95)]:
    for c in NOMES:
        X, w = AM[c][:2]; _, media, lim = indices(calibrar(X, w, P), X)
        linhas.append({"percentis": str(P), "cultura": NOMES[c], "apta_limitante_%": np.average(lim > LIMIAR, weights=w) * 100,
                       "apta_media_%": np.average(media > LIMIAR, weights=w) * 100})
tab_sg = pd.DataFrame(linhas).round(1); tab_sg.to_csv(f"{V4}/sensibilidade_global_v4.csv", index=False)
print("\nSENSIBILIDADE GLOBAL: % da área colhida classificada apta\n" + tab_sg.pivot(index="percentis", columns="cultura", values="apta_limitante_%").to_string())
print("\nCélula 2 concluída. Arquivos em", V4)
