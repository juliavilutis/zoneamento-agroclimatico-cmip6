# ═══ CÉLULA 4: baixar o ensemble CMIP6 (6 modelos, histórico + 3 cenários, 4 variáveis) ═══
import os, json, time, zipfile, shutil, glob, requests, xarray as xr

KEY  = "SEU-TOKEN-AQUI"
BASE = "https://cds.climate.copernicus.eu/api"
HEAD = {"PRIVATE-TOKEN": KEY, "Content-Type": "application/json"}
ENS  = "/content/drive/MyDrive/era5_brasil/cmip6/ensemble"
TMP  = "/content/tmp_ens"
MAX_ATIVOS = 6                                   # pedidos simultâneos na fila do CDS

MODELOS = ["miroc6", "gfdl_esm4", "ipsl_cm6a_lr", "mpi_esm1_2_lr", "mri_esm2_0", "ukesm1_0_ll"]
VARIAVEIS = {"tas": "near_surface_air_temperature", "tasmin": "daily_minimum_near_surface_air_temperature",
             "tasmax": "daily_maximum_near_surface_air_temperature", "pr": "precipitation"}
FUTURO = [str(a) for a in list(range(2041, 2061)) + list(range(2071, 2091))]
EXPERIMENTOS = {"historical": [str(a) for a in range(1991, 2015)],
                "ssp2_4_5": FUTURO, "ssp3_7_0": FUTURO, "ssp5_8_5": FUTURO}

os.makedirs(ENS, exist_ok=True); os.makedirs(TMP, exist_ok=True)
ESTADO = f"{ENS}/_pedidos.json"                  # guarda os pedidos em andamento, para retomar se o Colab cair
pedidos = json.load(open(ESTADO)) if os.path.exists(ESTADO) else {}
salvar_estado = lambda: json.dump(pedidos, open(ESTADO, "w"), indent=1)

tarefas = [(m, e, v) for m in MODELOS for e in EXPERIMENTOS for v in VARIAVEIS]
destino = lambda m, e, v: f"{ENS}/{m}/{v}_{e}.nc"
chave   = lambda m, e, v: f"{m}|{e}|{v}"
falhas  = {}

def enviar(m, e, v):
    corpo = {"inputs": {"model": m, "experiment": e, "variable": VARIAVEIS[v], "temporal_resolution": "monthly",
                        "year": EXPERIMENTOS[e], "month": [f"{i:02d}" for i in range(1, 13)],
                        "area": [6, -74, -34, -28], "data_format": "netcdf"}}
    r = requests.post(f"{BASE}/retrieve/v1/processes/projections-cmip6/execution", json=corpo, headers=HEAD, timeout=60)
    if r.status_code == 201: return r.json()["jobID"], None
    return None, f"{r.status_code}: {r.text[:200]}"

def baixar(job, m, e, v):
    res = requests.get(f"{BASE}/retrieve/v1/jobs/{job}/results", headers=HEAD, timeout=60).json()
    pasta = f"{TMP}/{m}_{e}_{v}"; shutil.rmtree(pasta, ignore_errors=True); os.makedirs(pasta)
    with requests.get(res["asset"]["value"]["href"], stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(f"{pasta}/dados.zip", "wb") as f:
            for parte in r.iter_content(1 << 20): f.write(parte)
    with zipfile.ZipFile(f"{pasta}/dados.zip") as z: z.extractall(pasta)
    arquivos = sorted(glob.glob(f"{pasta}/**/*.nc", recursive=True))
    dados = xr.concat([xr.open_dataset(a)[v] for a in arquivos], "time").sortby("time")   # junta se vier em mais de um arquivo
    os.makedirs(os.path.dirname(destino(m, e, v)), exist_ok=True)
    dados.to_dataset(name=v).to_netcdf(destino(m, e, v))
    shutil.rmtree(pasta, ignore_errors=True)

inicio = time.time()
while True:
    faltam = [t for t in tarefas if not os.path.exists(destino(*t)) and chave(*t) not in falhas]
    if not faltam: break
    ativos = [t for t in faltam if chave(*t) in pedidos]
    for t in faltam:                              # preenche as vagas livres da fila
        if len(ativos) >= MAX_ATIVOS: break
        if chave(*t) in pedidos: continue
        job, erro = enviar(*t)
        if job: pedidos[chave(*t)] = job; ativos.append(t); salvar_estado()
        else:
            print("Não enviado:", chave(*t), "|", erro)
            if erro.startswith(("401", "403")): raise SystemExit("Token recusado ou licença do CMIP6 não aceita no portal do CDS.")
            if erro.startswith("429"): break      # fila cheia: espera a próxima rodada
            falhas[chave(*t)] = erro
    for t in ativos:                              # confere o andamento
        job = pedidos[chave(*t)]
        try: status = requests.get(f"{BASE}/retrieve/v1/jobs/{job}", headers=HEAD, timeout=60).json().get("status")
        except Exception as ex: status = f"erro de conexão ({ex})"
        if status == "successful":
            try: baixar(job, *t); print(f"[{(time.time() - inicio) / 60:5.1f} min] salvo: {t[0]}/{t[2]}_{t[1]}.nc")
            except Exception as ex: falhas[chave(*t)] = f"download: {ex}"; print("Falha no download:", chave(*t), ex)
            pedidos.pop(chave(*t), None); salvar_estado()
        elif status in ("failed", "rejected", "dismissed"):
            falhas[chave(*t)] = status; pedidos.pop(chave(*t), None); salvar_estado(); print("Pedido", status, ":", chave(*t))
    time.sleep(30)

prontos = sum(os.path.exists(destino(*t)) for t in tarefas)
print(f"\nConcluído: {prontos} de {len(tarefas)} arquivos em {ENS}")
if falhas:
    print("Falharam (rode a célula de novo para tentar outra vez):")
    for k, v in falhas.items(): print("  ", k, "|", v)
