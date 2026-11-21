# ═══ PASSO 5: MIROC6 (tas e pr) para a análise de um modelo no passo 10 ═══
# Consolida as células de download do MIROC6 do notebook original (mesmos pedidos ao CDS).
# O ensemble de seis modelos, que inclui o MIROC6, é baixado no passo 11.
# Requer a licença do CMIP6 aceita no portal do CDS.
import os, time, zipfile, shutil, glob, requests

KEY  = "SEU-TOKEN-AQUI"                          # token pessoal do CDS (não publique o seu)
BASE = "https://cds.climate.copernicus.eu/api"
HEAD = {"PRIVATE-TOKEN": KEY, "Content-Type": "application/json"}
DIR  = "/content/drive/MyDrive/era5_brasil/cmip6"
TMP  = "/content/tmp_miroc6"
os.makedirs(DIR, exist_ok=True); os.makedirs(TMP, exist_ok=True)

VARIAVEIS = {"tas": "near_surface_air_temperature", "pr": "precipitation"}
PEDIDOS = {"historical_1991_2014": ("historical", range(1991, 2015))}
for ssp, cod in (("ssp245", "ssp2_4_5"), ("ssp585", "ssp5_8_5")):
    PEDIDOS[f"{ssp}_2041_2060"] = (cod, range(2041, 2061))
    PEDIDOS[f"{ssp}_2071_2090"] = (cod, range(2071, 2091))

def enviar(var, exp, anos):
    corpo = {"inputs": {"model": "miroc6", "experiment": exp, "variable": VARIAVEIS[var],
                        "temporal_resolution": "monthly", "year": [str(a) for a in anos],
                        "month": [f"{m:02d}" for m in range(1, 13)],
                        "area": [6, -74, -34, -28], "data_format": "netcdf"}}
    r = requests.post(f"{BASE}/retrieve/v1/processes/projections-cmip6/execution", json=corpo, headers=HEAD, timeout=60)
    r.raise_for_status()
    return r.json()["jobID"]

def esperar_e_baixar(job, destino, limite_min=60):
    for _ in range(limite_min * 2):
        status = requests.get(f"{BASE}/retrieve/v1/jobs/{job}", headers=HEAD, timeout=60).json().get("status")
        if status == "successful":
            res = requests.get(f"{BASE}/retrieve/v1/jobs/{job}/results", headers=HEAD, timeout=60).json()
            pasta = f"{TMP}/{os.path.basename(destino)}"; shutil.rmtree(pasta, ignore_errors=True); os.makedirs(pasta)
            with requests.get(res["asset"]["value"]["href"], stream=True, timeout=600) as r:
                r.raise_for_status()
                with open(f"{pasta}/dados.zip", "wb") as f:
                    for parte in r.iter_content(1 << 20): f.write(parte)
            with zipfile.ZipFile(f"{pasta}/dados.zip") as z: z.extractall(pasta)
            shutil.copy2(sorted(glob.glob(f"{pasta}/**/*.nc", recursive=True))[0], destino)
            shutil.rmtree(pasta, ignore_errors=True)
            print("  salvo:", destino); return True
        if status in ("failed", "rejected", "dismissed"):
            print("  pedido", status); return False
        time.sleep(30)
    print("  tempo esgotado; rode o script de novo"); return False

for var in VARIAVEIS:
    for nome, (exp, anos) in PEDIDOS.items():
        destino = f"{DIR}/{var}_{nome}.nc"
        if os.path.exists(destino): print("Já existe:", destino); continue
        print("Pedindo", var, nome, "...")
        esperar_e_baixar(enviar(var, exp, anos), destino)
print("Passo 5 concluído.")
