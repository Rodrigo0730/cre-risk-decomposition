# Spatial risk decomposition of London commercial real estate

Empirical application of spatial independent component analysis (spFOBI and spJADE) to commercial real estate transactions in Greater London, 2023. This repository is part of the author's Master's thesis.

The methods are introduced in Morales-Martínez, R., Nordhausen, K. and Ruiz, A. M. (2026). Spatial independent component analysis for heteroskedastic random fields. *Statistics & Probability Letters*, 238, 110858. https://doi.org/10.1016/j.spl.2026.110858

## Data
HM Land Registry Price Paid Data 2023, restricted to Greater London transactions of type "Other" or category "B". Postcodes are geocoded with postcodes.io. The data is not included; to rebuild it, run from the repository root:
1. `src/download.R`: downloads and filters the Price Paid Data into `data/processed/london_cre_2023.csv`
2. `src/geodata.R`: geocodes the postcodes into `data/processed/london_cre_2023_spatial.gpkg`

## Analysis
- `notebooks/data_prep.Rmd`: builds the location-level variables (median price, liquidity, local volatility) and `data/processed/londonCRE.gpkg`
- `notebooks/eda.Rmd`: exploratory spatial data analysis (maps, Moran's I)
- `notebooks/analysis.Rmd`: spatial kernels, SBSS, spFOBI and spJADE; maps of the independent components, loadings heatmap, Moran's I and kurtosis tables
- `notebooks/test.Rmd`: additional checks
- `notebooks/sar_spjade.R`, `notebooks/sar_spjade_n1000.R`: SAR model combined with spJADE (5,000- and 1,000-location samples)
- `src/spICA.R`: spFOBI and spJADE implementation used by the notebooks

Figures: `IC*_spfobi.pdf`, `IC*_spjade.pdf`, `loadings_heatmap_comparison.pdf`, `notebooks/combined_esda_maps.pdf`.
