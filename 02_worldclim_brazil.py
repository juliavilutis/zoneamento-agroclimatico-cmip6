# ═══ PASSO 2: WorldClim v2.1 (10 min) - Tmin e Tmax mensais recortadas para o Brasil ═══
# Consolida as células de download e recorte do notebook original.
import os, zipfile, requests, numpy as np, xarray as xr, rasterio
from rasterio.windows import from_bounds

DIR  = "/content/drive/MyDrive/era5_brasil"
BBOX = {"left": -74, "bottom": -34, "right": -28, "top": 6}
os.makedirs(DIR, exist_ok=True)

for var in ["tmin", "tmax"]:
    pasta = f"{DIR}/worldclim_{var}"
    if not os.path.exists(f"{pasta}/wc2.1_10m_{var}_12.tif"):
        print(f"Baixando WorldClim {var} ...")
        z = f"/tmp/worldclim_{var}.zip"
        with requests.get(f"https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_{var}.zip",
                          stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(z, "wb") as f:
                for parte in r.iter_content(1 << 20): f.write(parte)
        os.makedirs(pasta, exist_ok=True)
        with zipfile.ZipFile(z) as arq: arq.extractall(pasta)
        os.remove(z)

    meses = []
    for mes in range(1, 13):
        with rasterio.open(f"{pasta}/wc2.1_10m_{var}_{mes:02d}.tif") as src:
            janela = from_bounds(**BBOX, transform=src.transform)
            dados = src.read(1, window=janela).astype(np.float32)
            dados[dados <= -9999] = np.nan                         # nodata do WorldClim
            dados[dados < -200] = np.nan
            meses.append(dados)
            if mes == 1:
                t = src.window_transform(janela); alt, larg = dados.shape
                lons = np.array([t.c + t.a * (i + 0.5) for i in range(larg)])
                lats = np.array([t.f + t.e * (j + 0.5) for j in range(alt)])

    da = xr.DataArray(np.stack(meses), dims=["month", "lat", "lon"],
                      coords={"month": list(range(1, 13)), "lat": lats, "lon": lons},
                      name=var, attrs={"units": "degrees_C", "source": "WorldClim v2.1 10min"})
    da.to_netcdf(f"{DIR}/worldclim_{var}_brasil.nc")
    print("Salvo:", f"{DIR}/worldclim_{var}_brasil.nc")
print("Passo 2 concluído.")
