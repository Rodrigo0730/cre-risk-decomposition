library(tidyverse)

# direct URL for the 2023 Price Paid Data
url <- "http://prod.publicdata.landregistry.gov.uk.s3-website-eu-west-1.amazonaws.com/pp-2023.csv"
raw_dest <- "data/raw/pp-2023.csv"

cat("Downloading HM Land Registry Data...\n")
download.file(url, destfile = raw_dest, method = "auto")

# column headers
ppd_headers <- c(
  "transaction_id", "price", "date_of_transfer", "postcode",
  "property_type", "old_new", "duration", "paon", "saon",
  "street", "locality", "town_city", "district", "county",
  "ppd_category", "record_status"
)

cat("Loading data into memory...\n")
raw_data <- read_csv(raw_dest, col_names = ppd_headers, show_col_types = FALSE)

# "O" = Other (often commercial/land)
# "B" = Additional Price Paid entry (includes corporate/commercial transactions)
cat("Filtering for Greater London CRE...\n")
london_cre_data <- raw_data %>%
  filter(county == "GREATER LONDON") %>%
  filter(property_type == "O" | ppd_category == "B") %>%
  mutate(date_of_transfer = as.Date(date_of_transfer, format = "%Y-%m-%d"))

processed_dest <- "data/processed/london_cre_2023.csv"
write_csv(london_cre_data, processed_dest)

cat(sprintf(
  "Filtered dataset saved to %s. It contains %d transactions.\n",
  processed_dest, nrow(london_cre_data)
))
