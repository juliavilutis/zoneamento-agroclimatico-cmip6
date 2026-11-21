# ═══ CÉLULA 8: calibrar a variabilidade entre dias no cálculo das horas de frio ═══
import json, numpy as np, pandas as pd, xarray as xr
import matplotlib.pyplot as plt
from scipy import stats as sp

DIR = "/content/drive/MyDrive/era5_brasil"
V4  = f"{DIR}/v4"
est = pd.read_csv(f"{V4}/validacao_frio_inmet_v4.csv")

# ── Clima mensal no pixel de cada estação ───────────────────────────────────────
ds = xr.open_dataset(f"{DIR}/climatologia_completa_brasil.nc", engine="netcdf4", decode_timedelta=False)
lat, lon = ds.lat.values, ds.lon.values; tn, tx = ds.tmin.values, ds.tmax.values
valido = np.where(~np.isnan(tn[0]) & ~np.isnan(ds.tmean.values[0]))
pix = [((lat[valido[0]] - a) ** 2 + (lon[valido[1]] - b) ** 2).argmin() for a, b in zip(est.lat, est.lon)]
TN = np.array([tn[:, valido[0][q], valido[1][q]] for q in pix]).T        # (12, estações)
TX = np.array([tx[:, valido[0][q], valido[1][q]] for q in pix]).T

# ── Horas de frio: ciclo diário senoidal + anomalia do dia ~ Normal(0, sigma) ───
H    = np.cos(2 * np.pi * (np.arange(24) + 0.5) / 24)
DIAS = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
def horas_frio(tmin, tmax, sigma=0.0, sup=7.2, inf=0.0):
    """tmin, tmax: (12, ...) em °C. sigma: desvio-padrão da temperatura entre dias do mesmo mês (°C)."""
    med, amp = (tmin + tmax) / 2, (tmax - tmin) / 2
    total = np.zeros(tmin.shape[1:])
    for m in range(12):
        for h in H:
            temp = med[m] + amp[m] * h
            if sigma > 0: p = sp.norm.cdf((sup - temp) / sigma) - sp.norm.cdf((inf - temp) / sigma)
            else:         p = ((temp <= sup) & (temp > inf)).astype(float)
            total += DIAS[m] * p
    return total

obs = est.frio_obs.values
def metricas(pred, o):
    e = pred - o
    return {"vies_h": e.mean(), "erro_abs_medio_h": np.abs(e).mean(), "RMSE_h": np.sqrt((e ** 2).mean()), "r": np.corrcoef(pred, o)[0, 1]}

SIGMAS = np.arange(0, 6.01, 0.25)
PRED = {s: horas_frio(TN, TX, s) for s in SIGMAS}
curva = pd.DataFrame([{"sigma": s, **metricas(PRED[s], obs)} for s in SIGMAS]).round(2)
curva.to_csv(f"{V4}/frio_sigma_curva_v4.csv", index=False)
melhor = float(curva.loc[curva.RMSE_h.idxmin(), "sigma"])
print("ERRO EM FUNÇÃO DE SIGMA (todas as estações)\n" + curva.iloc[::2].to_string(index=False))

# ── Validação cruzada: sigma escolhido em 4/5 das estações, testado no 1/5 restante
rng = np.random.default_rng(0); dobra = rng.permutation(len(est)) % 5; teste = np.zeros(len(est)); escolhidos = []
for k in range(5):
    tr, te = dobra != k, dobra == k
    s_k = min(SIGMAS, key=lambda s: ((PRED[s][tr] - obs[tr]) ** 2).mean()); escolhidos.append(float(s_k))
    teste[te] = PRED[s_k][te]
frios = obs >= 50
resumo = pd.DataFrame([
    {"metodo": "senoidal sem variabilidade (sigma = 0)", "grupo": "todas", **metricas(PRED[0.0], obs)},
    {"metodo": f"senoidal com sigma = {melhor}", "grupo": "todas", **metricas(PRED[melhor], obs)},
    {"metodo": "validação cruzada (5 partes)", "grupo": "todas", **metricas(teste, obs)},
    {"metodo": "senoidal sem variabilidade (sigma = 0)", "grupo": "50 h ou mais", **metricas(PRED[0.0][frios], obs[frios])},
    {"metodo": f"senoidal com sigma = {melhor}", "grupo": "50 h ou mais", **metricas(PRED[melhor][frios], obs[frios])},
    {"metodo": "validação cruzada (5 partes)", "grupo": "50 h ou mais", **metricas(teste[frios], obs[frios])}]).round(2)
resumo.to_csv(f"{V4}/frio_sigma_resumo_v4.csv", index=False)
json.dump({"sigma": melhor, "sigmas_validacao_cruzada": escolhidos, "n_estacoes": int(len(est))}, open(f"{V4}/frio_sigma_v4.json", "w"), indent=1)
print(f"\nSigma com menor erro: {melhor} °C | escolhidos na validação cruzada: {escolhidos}\n\n" + resumo.to_string(index=False))

est["frio_sigma"] = PRED[melhor].round(1); est.to_csv(f"{V4}/validacao_frio_inmet_v4.csv", index=False)
print("\nAs 15 estações mais frias:\n" + est.sort_values("frio_obs", ascending=False).head(15)[["nome", "uf", "alt", "frio_obs", "frio_senoide", "frio_sigma"]].to_string(index=False))

fig, ax = plt.subplots(figsize=(7.5, 7.5)); lim = max(obs.max(), PRED[melhor].max()) * 1.08
ax.plot([0, lim], [0, lim], color="#52514e", linewidth=1, zorder=1)
ax.scatter(obs, PRED[0.0], s=42, color="#eb6834", edgecolor="white", linewidth=0.8, label="Sinusoidal cycle only", zorder=2)
ax.scatter(obs, PRED[melhor], s=42, color="#2a78d6", edgecolor="white", linewidth=0.8, label=f"With day-to-day variability (sigma = {melhor} °C)", zorder=3)
ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal"); ax.grid(color="#f0efec", linewidth=1); ax.set_axisbelow(True)
ax.set_xlabel("Observed chilling hours, INMET stations (h/yr)", fontsize=13); ax.set_ylabel("Estimated chilling hours (h/yr)", fontsize=13)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=12, loc="upper left")
fig.savefig(f"{V4}/figuras/fig6_validacao_frio.png", dpi=300, bbox_inches="tight", facecolor="white"); plt.show()
print("\nCélula 8 concluída.")
