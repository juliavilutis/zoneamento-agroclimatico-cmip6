# ═══ CÉLULA 7: validar as horas de frio com as estações automáticas do INMET ═══
import os, io, zipfile, requests, numpy as np, pandas as pd, xarray as xr
import matplotlib.pyplot as plt
from scipy import stats as sp

DIR = "/content/drive/MyDrive/era5_brasil"
V4  = f"{DIR}/v4"
TMP = "/content/tmp_inmet"; os.makedirs(TMP, exist_ok=True); os.makedirs(f"{V4}/figuras", exist_ok=True)
ANOS = range(2008, 2021)
UFS  = {"RS", "SC", "PR", "SP", "MG", "RJ", "ES", "MS"}
MIN_VALIDO, MIN_ANOS = 0.90, 5                   # fração mínima de horas válidas no ano; mínimo de anos por estação

# ── 1. Horas de frio observadas (0 < T <= 7,2 °C), por estação e ano ────────────
def ler_estacao(bruto):
    texto = bruto.decode("latin-1"); linhas = texto.split("\n", 15)
    ini = next(i for i, l in enumerate(linhas[:15]) if "BULBO SECO" in l.upper())    # linha com os nomes das colunas
    meta = {l.split(";")[0].upper(): l.split(";")[1].strip() for l in linhas[:ini] if ";" in l}
    achar = lambda p: next((v for k, v in meta.items() if p in k), "")
    num = lambda v: float(v.replace(",", ".")) if v not in ("", "F") else np.nan
    cab = linhas[ini].split(";"); col = next(i for i, c in enumerate(cab) if "BULBO SECO" in c.upper())
    t = pd.read_csv(io.StringIO(texto), sep=";", skiprows=ini, usecols=[col], decimal=",", dtype=str).iloc[:, 0]
    t = pd.to_numeric(t.str.replace(",", "."), errors="coerce"); t[t < -90] = np.nan
    return {"codigo": achar("CODIGO") or achar("CÓDIGO"), "nome": achar("ESTAC") or achar("ESTAÇ"), "uf": achar("UF"),
            "lat": num(achar("LATITUDE")), "lon": num(achar("LONGITUDE")), "alt": num(achar("ALTITUDE")),
            "valido": t.notna().mean(), "frio": float(((t > 0) & (t <= 7.2)).sum())}

registros, mostrou = [], False
for ano in ANOS:
    z = f"{TMP}/{ano}.zip"
    if not os.path.exists(z):
        with requests.get(f"https://portal.inmet.gov.br/uploads/dadoshistoricos/{ano}.zip", stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(z, "wb") as f:
                for parte in r.iter_content(1 << 22): f.write(parte)
    n = 0
    with zipfile.ZipFile(z) as arq:
        for nome in arq.namelist():
            base = os.path.basename(nome).upper().split("_")
            if not nome.upper().endswith(".CSV") or len(base) < 4 or base[2] not in UFS: continue
            bruto = arq.read(nome)
            if not mostrou: print("Exemplo de cabeçalho:\n" + bruto.decode("latin-1")[:420] + "\n"); mostrou = True
            try: registros.append({"ano": ano, **ler_estacao(bruto)}); n += 1
            except Exception as ex: print("  não lido:", nome, "|", ex)
    print(f"{ano}: {n} estações lidas")

reg = pd.DataFrame(registros); reg.to_csv(f"{V4}/frio_inmet_anual_v4.csv", index=False)
ok = reg[reg.valido >= MIN_VALIDO].copy(); ok["frio"] = ok.frio / ok.valido      # corrige pelas horas faltantes
est = ok.groupby("codigo").agg(nome=("nome", "last"), uf=("uf", "last"), lat=("lat", "last"), lon=("lon", "last"), alt=("alt", "last"),
                               anos=("ano", "count"), frio_obs=("frio", "mean")).reset_index()
est = est[(est.anos >= MIN_ANOS) & est.lat.notna()]

# ── 2. Horas de frio estimadas pela climatologia, no pixel de cada estação ──────
ds = xr.open_dataset(f"{DIR}/climatologia_completa_brasil.nc", engine="netcdf4", decode_timedelta=False)
lat, lon = ds.lat.values, ds.lon.values
tm, tn, tx = ds.tmean.values, ds.tmin.values, ds.tmax.values
DIAS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; H = np.cos(2 * np.pi * (np.arange(24) + 0.5) / 24)
senoide = np.zeros(tn.shape[1:]); normal = np.zeros(tn.shape[1:])
for m in range(12):
    for h in H:
        temp = (tn[m] + tx[m]) / 2 + (tx[m] - tn[m]) / 2 * h
        senoide += DIAS[m] * ((temp <= 7.2) & (temp > 0))
    sig = (tm[m] - tn[m]) / 2; sig = np.where(sig <= 0, 1.0, sig)                 # método original do artigo
    normal += DIAS[m] * 24 * (sp.norm.cdf(7.2, tm[m], sig) - sp.norm.cdf(0.0, tm[m], sig))
valido = np.where(~np.isnan(tn[0]) & ~np.isnan(tm[0]))
def no_pixel(grade, la, lo):
    q = ((lat[valido[0]] - la) ** 2 + (lon[valido[1]] - lo) ** 2).argmin(); return float(grade[valido[0][q], valido[1][q]])
est["frio_senoide"] = [no_pixel(senoide, a, b) for a, b in zip(est.lat, est.lon)]
est["frio_normal"]  = [no_pixel(normal,  a, b) for a, b in zip(est.lat, est.lon)]
est = est.round(1).sort_values("frio_obs", ascending=False); est.to_csv(f"{V4}/validacao_frio_inmet_v4.csv", index=False)

# ── 3. Estatísticas e figura ────────────────────────────────────────────────────
def metricas(sub, col):
    e = sub[col] - sub.frio_obs
    return {"vies_h": e.mean(), "erro_abs_medio_h": e.abs().mean(), "RMSE_h": np.sqrt((e ** 2).mean()), "r": np.corrcoef(sub[col], sub.frio_obs)[0, 1]}
linhas = []
for rotulo, sub in (("todas as estações", est), ("estações com 50 h ou mais", est[est.frio_obs >= 50])):
    for nome, col in (("senoidal (novo)", "frio_senoide"), ("distribuição normal (original)", "frio_normal")):
        linhas.append({"grupo": rotulo, "n": len(sub), "metodo": nome, **metricas(sub, col)})
tab = pd.DataFrame(linhas).round(2); tab.to_csv(f"{V4}/validacao_frio_resumo_v4.csv", index=False)
print(f"\n{len(est)} estações com pelo menos {MIN_ANOS} anos válidos\n\n" + tab.to_string(index=False))
print("\nAs 15 estações mais frias:\n" + est.head(15)[["nome", "uf", "alt", "anos", "frio_obs", "frio_senoide", "frio_normal"]].to_string(index=False))

fig, ax = plt.subplots(figsize=(7.5, 7.5)); lim = max(est.frio_obs.max(), est.frio_senoide.max()) * 1.08
ax.plot([0, lim], [0, lim], color="#52514e", linewidth=1, zorder=1)
ax.scatter(est.frio_obs, est.frio_normal, s=42, color="#eb6834", edgecolor="white", linewidth=0.8, label="Normal distribution (original)", zorder=2)
ax.scatter(est.frio_obs, est.frio_senoide, s=42, color="#2a78d6", edgecolor="white", linewidth=0.8, label="Sinusoidal diurnal cycle (new)", zorder=3)
ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal"); ax.grid(color="#f0efec", linewidth=1); ax.set_axisbelow(True)
ax.set_xlabel("Observed chilling hours, INMET stations (h/yr)", fontsize=13); ax.set_ylabel("Estimated chilling hours (h/yr)", fontsize=13)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=12, loc="upper left")
fig.savefig(f"{V4}/figuras/fig6_validacao_frio.png", dpi=300, bbox_inches="tight", facecolor="white"); plt.show()
print("\nCélula 7 concluída. Arquivos em", V4)
