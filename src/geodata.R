library(tidyverse)
library(httr)
library(jsonlite)
library(sf)

cat("Loading processed CRE dataset...\n")
cre_data <- read_csv("data/processed/london_cre_2023.csv", show_col_types = FALSE)

unique_postcodes <- cre_data %>%
  filter(!is.na(postcode)) %>%
  pull(postcode) %>%
  unique()

cat(sprintf("Found %d unique postcodes. Beginning bulk geocoding...\n", length(unique_postcodes)))

fetch_coords <- function(postcode_batch) {
  res <- POST("https://api.postcodes.io/postcodes",
    body = list(postcodes = postcode_batch),
    encode = "json"
  )

  parsed <- fromJSON(content(res, "text", encoding = "UTF-8"), simplifyVector = FALSE)

  results <- parsed$result

  if (is.null(results)) {
    return(NULL)
  }

  tibble(
    postcode = map_chr(results, "query"),
    longitude = map_dbl(results, ~ if (is.null(.x$result) || is.null(.x$result$longitude)) NA_real_ else as.numeric(.x$result$longitude)),
    latitude = map_dbl(results, ~ if (is.null(.x$result) || is.null(.x$result$latitude)) NA_real_ else as.numeric(.x$result$latitude))
  )
}

batches <- split(unique_postcodes, ceiling(seq_along(unique_postcodes) / 100))

cat("Querying postcodes.io API\n")
coords_list <- lapply(batches, fetch_coords)
coords_df <- bind_rows(coords_list)

cre_spatial_df <- cre_data %>%
  left_join(coords_df, by = "postcode") %>%
  filter(!is.na(longitude) & !is.na(latitude))

cre_sf <- st_as_sf(cre_spatial_df, coords = c("longitude", "latitude"), crs = 4326)

spatial_dest <- "data/processed/london_cre_2023_spatial.gpkg"
st_write(cre_sf, spatial_dest, delete_dsn = TRUE, quiet = TRUE)

cat(sprintf("Dataset saved to %s.\n", spatial_dest))
