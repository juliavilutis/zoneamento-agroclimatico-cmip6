# ═══ PASSO 3: climatologia de referência no grid do ERA5-Land (0,1°) ═══
# Consolida as duas células do notebook original que geraram climatologia_completa_brasil.nc.
# Variáveis usadas nos passos seguintes: tmean e precip (ERA5-Land, 1991-2020) e
# tmin e tmax (WorldClim v2.1, 1970-2000, interpoladas bilinearmente para o grid do ERA5-Land).
# As variáveis et, deficit e chilling_hours são mantidas por fidelidade ao arquivo original,
# mas não são usadas: as horas de frio são recalculadas nos passos 7 a 12.
import os, numpy as np, xarray as xr
from scipy import stats as sp_stats

DIR = "/content/drive/MyDrive/era5_brasil"

ds_tmean  = xr.open_dataset(f"{DIR}/brasil_t2m_mean_1991_2020.nc", engine="netcdf4")
ds_precip = xr.open_dataset(f"{DIR}/brasil_precip_1991_2020.nc",   engine="netcdf4")
ds_et     = xr.open_dataset(f"{DIR}/brasil_et_1991_2020.nc",       engine="netcdf4")

tmean  = ds_tmean["t2m"] - 273.15                 # K para °C
precip = ds_precip["tp"] * 1000 * 30              # m/dia para mm/mês (30 dias por mês, como no processamento original)
et     = np.abs(ds_et["pev"]) * 1000 * 30         # pev é negativo no ERA5 (convenção de fluxo)

tmean_clim  = tmean.groupby("valid_time.month").mean("valid_time")
precip_clim = precip.groupby("valid_time.month").mean("valid_time")
et_clim     = et.groupby("valid_time.month").mean("valid_time")
renomear = lambda d: d.rename({"latitude": "lat", "longitude": "lon"})
tmean_clim, precip_clim, et_clim = renomear(tmean_clim), renomear(precip_clim), renomear(et_clim)

tmin = xr.open_dataset(f"{DIR}/worldclim_tmin_brasil.nc", engine="netcdf4")["tmin"]
tmax = xr.open_dataset(f"{DIR}/worldclim_tmax_brasil.nc", engine="netcdf4")["tmax"]
tmin_i = tmin.interp(lat=tmean_clim.lat, lon=tmean_clim.lon, method="linear")
tmax_i = tmax.interp(lat=tmean_clim.lat, lon=tmean_clim.lon, method="linear")

# Horas de frio da versão original (distribuição normal); não usadas nos resultados finais
dias = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
ch = []
for i, m in enumerate(range(1, 13)):
    tm, tn = tmean_clim.sel(month=m).values, tmin_i.sel(month=m).values
    s = np.where((tm - tn) / 2 <= 0, 1.0, (tm - tn) / 2)
    ch.append(dias[i] * 24 * (sp_stats.norm.cdf(7.2, tm, s) - sp_stats.norm.cdf(0.0, tm, s)))
ch = np.stack(ch)

coords = {"month": range(1, 13), "lat": tmean_clim.lat, "lon": tmean_clim.lon}
ds = xr.Dataset({
    "tmean": tmean_clim, "tmin": tmin_i, "tmax": tmax_i, "precip": precip_clim,
    "et": et_clim, "deficit": et_clim - precip_clim,
    "chilling_hours": xr.DataArray(ch, dims=["month", "lat", "lon"], coords=coords),
    "chilling_hours_anual": xr.DataArray(ch.sum(0), dims=["lat", "lon"],
                                         coords={"lat": tmean_clim.lat, "lon": tmean_clim.lon}),
})
out = f"{DIR}/climatologia_completa_brasil.nc"
ds.to_netcdf(out)
print(f"Salvo: {out} ({os.path.getsize(out) / 1e6:.1f} MB)")
