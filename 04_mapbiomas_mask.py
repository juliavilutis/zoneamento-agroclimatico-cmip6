# ═══ PASSO 4: máscara de uso agrícola (MapBiomas Coleção 9, 2023) no grid de 0,1° ═══
# Célula do notebook original. O mapa de cobertura de 2023 (30 m, GeoTIFF) precisa ser
# baixado manualmente na plataforma do MapBiomas (https://brasil.mapbiomas.org) e salvo em DIR.
# Cada célula de 0,1° recebe a classe majoritária dos pixels de 30 m que contém.
# Demora de 15 a 30 minutos.
import numpy as np, xarray as xr, rasterio
from rasterio.windows import from_bounds
from scipy import stats as sp_stats

DIR = "/content/drive/MyDrive/era5_brasil"
MAPBIOMAS_TIF = f"{DIR}/2023_coverage_lclu_1-1-1_2242342e-4c61-44e4-83ac-2a2e7e23b148.tif"   # ajuste ao nome do seu arquivo

# Pastagem, lavouras temporárias e perenes, mosaico de usos, silvicultura e outras áreas não vegetadas
CLASSES_AGRICOLAS = {9, 15, 20, 21, 25, 35, 39, 40, 41, 46, 47, 48, 62, 75}

clim = xr.open_dataset(f"{DIR}/climatologia_completa_brasil.nc", engine="netcdf4")
lat, lon = clim["lat"].values, clim["lon"].values
res = abs(lat[1] - lat[0])

classes = np.zeros((len(lat), len(lon)), dtype=np.int16)
mascara = np.zeros((len(lat), len(lon)), dtype=np.uint8)
with rasterio.open(MAPBIOMAS_TIF) as src:
    nodata = src.nodata or 0
    for i, la in enumerate(lat):
        if i % 20 == 0: print(f"  {i}/{len(lat)} linhas", end="\r")
        for j, lo in enumerate(lon):
            try:
                d = src.read(1, window=from_bounds(lo - res / 2, la - res / 2, lo + res / 2, la + res / 2, src.transform))
                d = d[(d != nodata) & (d > 0)]
                if len(d) == 0: continue
                moda = int(sp_stats.mode(d, keepdims=True).mode[0])
                classes[i, j] = moda
                mascara[i, j] = 1 if moda in CLASSES_AGRICOLAS else 0
            except Exception:
                continue

print(f"\nCélulas agrícolas: {int(mascara.sum()):,}")
coords = {"lat": lat, "lon": lon}
xr.Dataset({"mascara": xr.DataArray(mascara, dims=["lat", "lon"], coords=coords),
            "classe":  xr.DataArray(classes, dims=["lat", "lon"], coords=coords)}
           ).to_netcdf(f"{DIR}/mapbiomas_mask_era5.nc")
print("Salvo:", f"{DIR}/mapbiomas_mask_era5.nc")
