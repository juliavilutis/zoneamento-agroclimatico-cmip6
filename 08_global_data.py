# ═══ CÉLULA 1 de 3: baixar dados globais (WorldClim v2.1 e CROPGRIDS v1.08) ═══
import os, zipfile, requests, numpy as np, rasterio

DIR = "/content/drive/MyDrive/era5_brasil"
V4  = f"{DIR}/v4"
TMP = "/content/tmp_v4"
os.makedirs(f"{V4}/cropgrids", exist_ok=True); os.makedirs(TMP, exist_ok=True)

def baixar(url, destino):
    if os.path.exists(destino): return
    print("Baixando", os.path.basename(destino), "...")
    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(destino, "wb") as f:
            for parte in r.iter_content(1 << 22): f.write(parte)

# ── WorldClim v2.1, 10 minutos: média, mínima, máxima e chuva mensais ───────────
def ler_worldclim(var):
    z = f"{TMP}/wc_{var}.zip"
    baixar(f"https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_{var}.zip", z)
    with zipfile.ZipFile(z) as arq: arq.extractall(f"{TMP}/wc")
    meses = []
    for m in range(1, 13):
        with rasterio.open(f"{TMP}/wc/wc2.1_10m_{var}_{m:02d}.tif") as r:
            a = r.read(1).astype("float32"); a[a < -1e30] = np.nan; meses.append(a)
    return np.stack(meses)                      # (12, 1080, 2160), latitude decrescente

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

if not os.path.exists(f"{V4}/worldclim_vars.npz"):
    tavg, tmin, tmax, prec = (ler_worldclim(v) for v in ("tavg", "tmin", "tmax", "prec"))
    np.savez_compressed(f"{V4}/worldclim_vars.npz",
                        tmean=tavg.mean(0), tmin=tmin.min(0), tmax=tmax.max(0),
                        precip=prec.sum(0), frio=horas_frio(tmin, tmax))
print("OK:", f"{V4}/worldclim_vars.npz")

# ── CROPGRIDS v1.08 (Tang et al., 2024): área colhida global por cultura ────────
CULT_CG = ["hop", "pistachio", "hazelnut", "almond"]
if not all(os.path.exists(f"{V4}/cropgrids/CROPGRIDSv1.08_{c}.nc") for c in CULT_CG):
    z = f"{TMP}/cropgrids.zip"
    baixar("https://ndownloader.figshare.com/files/44950942", z)      # 807 MB
    with zipfile.ZipFile(z) as arq:
        for c in CULT_CG:
            nome = f"CROPGRIDSv1.08_NC_maps/CROPGRIDSv1.08_{c}.nc"
            with arq.open(nome) as origem, open(f"{V4}/cropgrids/CROPGRIDSv1.08_{c}.nc", "wb") as dest:
                dest.write(origem.read())
    os.remove(z)
print("OK:", sorted(os.listdir(f"{V4}/cropgrids")))

# ── Coordenadas dos municípios brasileiros (para a validação do lúpulo) ─────────
for arq in ("municipios.csv", "estados.csv"):
    baixar(f"https://raw.githubusercontent.com/kelvins/municipios-brasileiros/main/csv/{arq}", f"{V4}/{arq}")
print("Célula 1 concluída.")
