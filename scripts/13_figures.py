# ═══ CÉLULA 6: figuras do artigo (a partir dos arquivos em v4) ═══
import os, glob, json, zipfile, requests, unicodedata, numpy as np, pandas as pd, xarray as xr, geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

DIR = "/content/drive/MyDrive/era5_brasil"
V4  = f"{DIR}/v4"
FIG = f"{V4}/figuras"; os.makedirs(FIG, exist_ok=True)
IDIOMA = "en"                                    # "en" para a revista, "pt" para a versão em português

T = {"en": {"Lúpulo": "Hops", "Pistache": "Pistachio", "Avelã": "Hazelnut", "Amêndoa": "Almond",
            "tmean": "Mean temperature", "tmin": "Coldest-month Tmin", "tmax": "Warmest-month Tmax", "precip": "Annual precipitation", "frio": "Chilling hours",
            "n_crit": "Climatic criteria met (of 5)", "excl": "Excluded by this criterion", "nao_excl": "Within the envelope",
            "n_mod": "Models placing the pixel within the thermal envelope (of 6)", "atual": "Current (1991-2020)",
            "area": "Thermally suitable area (Mha)", "media": "Ensemble mean", "modelo": "Individual model", "faixa": "Model range",
            "indice": "Climatic suitability index", "mun": "Hop-producing municipality", "nota": "Map lines delineate study areas and do not necessarily depict accepted national boundaries."},
     "pt": {"Lúpulo": "Lúpulo", "Pistache": "Pistache", "Avelã": "Avelã", "Amêndoa": "Amêndoa",
            "tmean": "Temperatura média", "tmin": "Tmin do mês mais frio", "tmax": "Tmax do mês mais quente", "precip": "Precipitação anual", "frio": "Horas de frio",
            "n_crit": "Critérios climáticos atendidos (de 5)", "excl": "Excluído por este critério", "nao_excl": "Dentro do envelope",
            "n_mod": "Modelos que colocam o pixel no envelope térmico (de 6)", "atual": "Atual (1991-2020)",
            "area": "Área termicamente apta (Mha)", "media": "Média do ensemble", "modelo": "Modelo individual", "faixa": "Faixa entre modelos",
            "indice": "Índice de aptidão climática", "mun": "Município produtor de lúpulo", "nota": "As linhas do mapa delimitam a área de estudo e não representam necessariamente fronteiras oficiais."}}[IDIOMA]
CULTURAS = ["Lúpulo", "Pistache", "Avelã", "Amêndoa"]
VARS = ["tmean", "tmin", "tmax", "precip", "frio"]
TINTA, CINZA, FUNDO = "#0b0b0b", "#52514e", "#f0efec"
AZUIS = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]   # sequencial, um matiz
SERIES = {"SSP2-4.5": "#2a78d6", "SSP3-7.0": "#eb6834", "SSP5-8.5": "#1baf7a"}          # categórica, validada para daltonismo
plt.rcParams.update({"font.size": 13, "axes.titlesize": 15, "axes.labelsize": 13, "axes.edgecolor": CINZA, "text.color": TINTA,
                     "axes.labelcolor": TINTA, "xtick.color": CINZA, "ytick.color": CINZA, "figure.facecolor": "white", "savefig.facecolor": "white"})

def natural_earth(pasta, url):
    if not glob.glob(f"/tmp/{pasta}/*.shp"):
        r = requests.get(url, stream=True, timeout=120)
        with open(f"/tmp/{pasta}.zip", "wb") as f:
            for parte in r.iter_content(1 << 20): f.write(parte)
        with zipfile.ZipFile(f"/tmp/{pasta}.zip") as z: z.extractall(f"/tmp/{pasta}/")
    return gpd.read_file(glob.glob(f"/tmp/{pasta}/*.shp")[0])
estados = natural_earth("estados_shp", "https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_1_states_provinces.zip")
estados = estados[estados["iso_a2"] == "BR"]
paises  = natural_earth("paises_shp", "https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip")

def mapa_br(ax, dados, lat, lon, cmap, norm, titulo, zoom=None):
    ax.pcolormesh(lon, lat, np.ma.masked_invalid(dados), cmap=cmap, norm=norm, shading="nearest", rasterized=True)
    estados.boundary.plot(ax=ax, linewidth=0.5, color=CINZA)
    ax.set_xlim(*(zoom[:2] if zoom else (-74.5, -34))); ax.set_ylim(*(zoom[2:] if zoom else (-34.5, 6)))
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(""); ax.set_ylabel(""); ax.set_title(titulo, fontweight="bold")
    for s in ax.spines.values(): s.set_visible(False)

def discreto(cores): return ListedColormap(cores), BoundaryNorm(np.arange(-0.5, len(cores)), len(cores))
def legenda_discreta(fig, cores, rotulos, titulo, y=0.02, ncol=None):
    fig.legend(handles=[Patch(facecolor=c, edgecolor=CINZA, linewidth=0.5, label=r) for c, r in zip(cores, rotulos)],
               loc="lower center", bbox_to_anchor=(0.5, y), ncol=ncol or len(cores), frameon=False, title=titulo, fontsize=13, title_fontsize=14)
def salvar(fig, nome):
    fig.savefig(f"{FIG}/{nome}.png", dpi=300, bbox_inches="tight"); plt.show(); print("salvo:", f"{FIG}/{nome}.png")

fl = xr.open_dataset(f"{V4}/fatores_limitantes_v4.nc"); lat, lon = fl.lat.values, fl.lon.values

# ── Figura 1: número de critérios climáticos atendidos, clima atual ─────────────
sem_acento = lambda t: unicodedata.normalize("NFKD", t).encode("ascii", "ignore").decode().lower().strip()
mun = pd.read_csv(f"{V4}/municipios.csv").merge(pd.read_csv(f"{V4}/estados.csv", encoding="utf-8-sig")[["codigo_uf", "uf"]], on="codigo_uf")
mun["chave"] = mun.nome.map(sem_acento) + "|" + mun.uf
lup = pd.read_csv(f"{V4}/lupulo_municipios_v4.csv"); lup["chave"] = lup.municipio.map(sem_acento) + "|" + lup.uf
lup = lup.merge(mun[["chave", "latitude", "longitude"]], on="chave")
cores6 = [FUNDO] + AZUIS[1:6]; cmap, norm = discreto(cores6)
fig, axes = plt.subplots(2, 2, figsize=(13, 13.5))
for ax, c in zip(axes.flat, CULTURAS):
    mapa_br(ax, fl[f"{c}_criterios_atendidos"].values, lat, lon, cmap, norm, T[c])
    if c == "Lúpulo":
        ax.scatter(lup.longitude, lup.latitude, s=22, facecolor="white", edgecolor=TINTA, linewidth=0.9, zorder=5)
legenda_discreta(fig, cores6, [str(i) for i in range(6)], T["n_crit"], y=0.045)
fig.legend(handles=[Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor=TINTA, markersize=7, label=T["mun"])],
           loc="lower center", bbox_to_anchor=(0.5, 0.012), frameon=False, fontsize=13)
fig.text(0.5, 0.0, T["nota"], ha="center", fontsize=10, color=CINZA)
fig.subplots_adjust(bottom=0.13, wspace=0.03, hspace=0.08); salvar(fig, "fig1_criterios_atendidos")

# ── Figura 2: qual critério exclui cada região (culturas x critérios) ───────────
cmap2, norm2 = discreto([FUNDO, AZUIS[4]])
fig, axes = plt.subplots(4, 5, figsize=(17, 14.5))
for i, c in enumerate(CULTURAS):
    for j, k in enumerate(VARS):
        mapa_br(axes[i, j], fl[f"{c}_zerado_{k}"].values, lat, lon, cmap2, norm2, T[k] if i == 0 else "")
    axes[i, 0].text(-0.06, 0.5, T[c], transform=axes[i, 0].transAxes, rotation=90, va="center", ha="right", fontsize=16, fontweight="bold")
legenda_discreta(fig, [AZUIS[4], FUNDO], [T["excl"], T["nao_excl"]], None, y=0.03)
fig.text(0.5, 0.015, T["nota"], ha="center", fontsize=10, color=CINZA)
fig.subplots_adjust(bottom=0.08, wspace=0.02, hspace=0.03); salvar(fig, "fig2_fatores_limitantes")

# ── Figura 3: concordância entre modelos, Sul do Brasil ─────────────────────────
co = xr.open_dataset(f"{V4}/ensemble_concordancia_v4.nc")
cores7 = [FUNDO] + AZUIS[1:]; cmap3, norm3 = discreto(cores7); SUL = (-60, -38, -34.5, -13)
colunas = [("atual", T["atual"]), ("SSP2-4.5|2041-2060", "SSP2-4.5, 2041-2060"), ("SSP2-4.5|2071-2090", "SSP2-4.5, 2071-2090"), ("SSP5-8.5|2071-2090", "SSP5-8.5, 2071-2090")]
fig, axes = plt.subplots(2, 4, figsize=(17, 9.5))
for i, c in enumerate(["Pistache", "Amêndoa"]):
    for j, (chave, titulo) in enumerate(colunas):
        d = co[f"{c}|{chave}"].values * (6 if chave == "atual" else 1)       # o mapa atual é sim/não: usa a cor do máximo
        mapa_br(axes[i, j], d, lat, lon, cmap3, norm3, titulo if i == 0 else "", zoom=SUL)
    axes[i, 0].text(-0.05, 0.5, T[c], transform=axes[i, 0].transAxes, rotation=90, va="center", ha="right", fontsize=16, fontweight="bold")
legenda_discreta(fig, cores7, [str(i) for i in range(7)], T["n_mod"], y=0.035)
fig.text(0.5, 0.01, T["nota"], ha="center", fontsize=10, color=CINZA)
fig.subplots_adjust(bottom=0.15, wspace=0.03, hspace=0.04); salvar(fig, "fig3_concordancia_modelos")

# ── Figura 4: área termicamente apta por cenário, com a dispersão entre modelos ─
em = pd.read_csv(f"{V4}/ensemble_modelos_v4.csv"); rs = pd.read_csv(f"{V4}/ensemble_resumo_v4.csv")
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
for ax, c in zip(axes, ["Pistache", "Amêndoa"]):
    atual = float(rs[rs.cultura == c].termica_atual.iloc[0])
    for d, (cen, cor) in zip((-0.07, 0, 0.07), SERIES.items()):
        sub = em[(em.cultura == c) & (em.cenario == cen)]
        medias = [atual] + [sub[sub.periodo == p].termica_MapBiomas_Mha.mean() for p in ("2041-2060", "2071-2090")]
        ax.plot([0, 1 + d, 2 + d], medias, color=cor, linewidth=2, marker="o", markersize=8, zorder=3)
        for x, p in ((1 + d, "2041-2060"), (2 + d, "2071-2090")):
            v = sub[sub.periodo == p].termica_MapBiomas_Mha
            ax.plot([x, x], [v.min(), v.max()], color=cor, linewidth=1.2, alpha=0.6, zorder=2)
            ax.scatter([x] * len(v), v, s=14, color=cor, alpha=0.55, zorder=2)
        ax.annotate(cen, (2.07, medias[-1]), xytext=(12, {"SSP2-4.5": 14, "SSP3-7.0": 4, "SSP5-8.5": -6}[cen]), textcoords="offset points", va="center", fontsize=12, color=TINTA)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels([T["atual"], "2041-2060", "2071-2090"]); ax.set_xlim(-0.25, 2.65); ax.set_ylim(0, None)
    ax.set_title(T[c], fontweight="bold"); ax.grid(axis="y", color=FUNDO, linewidth=1); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[0].set_ylabel(T["area"])
fig.legend(handles=[Line2D([], [], color=cor, linewidth=2, marker="o", markersize=8, label=cen) for cen, cor in SERIES.items()] +
                   [Line2D([], [], color=CINZA, linewidth=2, marker="o", markersize=8, label=T["media"]),
                    Line2D([], [], color=CINZA, linestyle="", marker="o", markersize=4, alpha=0.55, label=T["modelo"])],
           loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.06))
fig.tight_layout(); salvar(fig, "fig4_area_por_cenario")

# ── Figura 5: índice de aptidão calibrado no mundo ──────────────────────────────
WC = dict(np.load(f"{V4}/worldclim_vars.npz")); PAR = json.load(open(f"{V4}/parametros_v4.json"))
def trap(v, a, b, c, d):
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(v <= a, 0., np.where(v < b, (v - a) / (b - a), np.where(v <= c, 1., np.where(v < d, (d - v) / (d - c), 0.)))).astype(float)
classes = [0.001, 0.25, 0.5, 0.75, 1.0001]; cores5 = [FUNDO, AZUIS[1], AZUIS[3], AZUIS[4], AZUIS[6]]
cmap5 = ListedColormap(cores5); norm5 = BoundaryNorm([-1] + classes, 5)
fig, axes = plt.subplots(2, 2, figsize=(17, 9.5))
for ax, c in zip(axes.flat, CULTURAS):
    s = [trap(WC[k], *PAR[c][k]) for k in VARS]; ind = np.where(np.minimum.reduce(s) > 0, sum(s) / 5, 0.); ind[np.isnan(WC["tmean"])] = np.nan
    ax.imshow(np.ma.masked_invalid(ind), extent=(-180, 180, -90, 90), cmap=cmap5, norm=norm5, interpolation="nearest", rasterized=True)
    paises.boundary.plot(ax=ax, linewidth=0.35, color=CINZA)
    ax.set_xlim(-170, 180); ax.set_ylim(-58, 75); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(""); ax.set_ylabel(""); ax.set_title(T[c], fontweight="bold")
    for sp in ax.spines.values(): sp.set_visible(False)
legenda_discreta(fig, cores5, ["0", "0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1"], T["indice"], y=0.03)
fig.text(0.5, 0.01, T["nota"], ha="center", fontsize=10, color=CINZA)
fig.subplots_adjust(bottom=0.13, wspace=0.02, hspace=0.08); salvar(fig, "fig5_indice_global")
print("\nCélula 6 concluída. Figuras em", FIG)
