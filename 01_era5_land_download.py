# ═══ PASSO 1: baixar o ERA5-Land mensal (1991-2020) para o retângulo do Brasil ═══
# Consolida as células de download do notebook original: mesmos pedidos ao CDS,
# com envio e espera no mesmo script. Requer conta no Copernicus CDS e a licença
# do ERA5-Land aceita no portal. Rode no Colab com o Google Drive montado:
#   from google.colab import drive; drive.mount("/content/drive")
import os, time, zipfile, shutil, glob, requests

KEY  = "SEU-TOKEN-AQUI"                          # token pessoal do CDS (não publique o seu)
BASE = "https://cds.climate.copernicus.eu/api"
HEAD = {"PRIVATE-TOKEN": KEY, "Content-Type": "application/json"}
DIR  = "/content/drive/MyDrive/era5_brasil"
TMP  = "/content/tmp_era5"
os.makedirs(DIR, exist_ok=True); os.makedirs(TMP, exist_ok=True)

BASE_PAYLOAD = {
    "product_type": "monthly_averaged_reanalysis",
    "year": [str(a) for a in range(1991, 2021)],
    "month": [f"{m:02d}" for m in range(1, 13)],
    "time": "00:00",
    "area": [6, -74, -34, -28],                  # norte, oeste, sul, leste
    "data_format": "netcdf",
}
VARIAVEIS = {                                    # nome do arquivo: variável do ERA5-Land
    "t2m_mean": "2m_temperature",
    "precip":   "total_precipitation",
    "et":       "potential_evaporation",
}

def enviar(variavel):
    r = requests.post(f"{BASE}/retrieve/v1/processes/reanalysis-era5-land-monthly-means/execution",
                      json={"inputs": {**BASE_PAYLOAD, "variable": variavel}}, headers=HEAD, timeout=60)
    r.raise_for_status()
    return r.json()["jobID"]

def esperar_e_baixar(job, destino, limite_min=120):
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

for nome, variavel in VARIAVEIS.items():
    destino = f"{DIR}/brasil_{nome}_1991_2020.nc"
    if os.path.exists(destino): print("Já existe:", destino); continue
    print("Pedindo", variavel, "...")
    esperar_e_baixar(enviar(variavel), destino)
print("Passo 1 concluído.")
