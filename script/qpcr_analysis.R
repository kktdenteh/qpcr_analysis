# load libraries
library(tidyverse)
library(ggpubr)


# Import the qpcr data
# We use 'skip' to skip the first 27 lines as they are just metadata
qpcr_data <- read.csv(
  file = "data/qpcr_data.csv",
  skip = 27
)

# View the qpcr data
View(qpcr_data)
glimpse(qpcr_data)

# Tidy the data
# We select only the columns that we need for the analysis
qpcr_data %>%
  janitor::clean_names() %>%
  select(well, sample_name, ct)
# Remove rows where sample_name is missing and ct value is undetermined
