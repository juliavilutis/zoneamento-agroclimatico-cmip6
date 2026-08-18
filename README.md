# Zoneamento Agroclimático de Culturas de Alto Valor sob Cenários CMIP6

Pipeline computacional em Python para mapeamento de aptidão climática de culturas especiais de alto valor econômico (lúpulo, avelã, amêndoa, pistache, açafrão e trufas) no território brasileiro, sob clima presente (baseline 1991–2020) e projeções futuras CMIP6.

## Funcionalidades
- Download e processamento de dados **ERA5-Land** via API do Copernicus Climate Data Store (temperatura, precipitação, evapotranspiração potencial)
- Integração com **WorldClim v2.1** (temperaturas mínimas e máximas mensais) com reamostragem bilinear para grid comum de 0,1°
- Estimativa de **horas de frio** por modelo probabilístico (Richardson et al., 1974) via CDF normal
- Projeções **CMIP6 (MIROC6)** com correção de viés pelo método delta
- Funções de pertinência **fuzzy trapezoidais** para índice de aptidão multicritério
- Filtragem de uso e cobertura do solo com **MapBiomas Coleção 9**
- Geração de mapas e figuras para publicação (300 DPI)

## Stack
`Python` · `xarray` · `rasterio` · `geopandas` · `numpy` · `scipy` · `matplotlib` · Google Colab

## Autora
Julia Vilutis Bronzeri — Programa de Pós-Graduação em Agricultura e Ambiente, UFSCar Araras.
Pipeline vinculado a manuscrito científico e a resumo submetido ao AGU Annual Meeting 2026.

*Os dados de entrada são públicos (ERA5-Land, WorldClim, CMIP6, MapBiomas); nenhum dado proprietário é distribuído neste repositório.*
