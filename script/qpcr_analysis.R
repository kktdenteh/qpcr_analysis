# load libraries
library(tidyverse)
library(ggpubr)


# Import the qpcr data
# We use 'skip' to skip the first 27 lines as they are just metadata
qpcr_data <- read.csv(file = "data/qpcr_data.csv",
                      skip = 27)

# View the qpcr data
View(qpcr_data)
glimpse(qpcr_data)

# Tidy the data
# We select only the columns that we need for the analysis
tidy_data <- qpcr_data %>% 
  janitor::clean_names() %>% 
  select(well, sample_name, ct) %>% 
# Remove rows where sample_name is missing or NTC and ct value is undetermined
  filter(sample_name != "" & sample_name != "NTC", ct != "Undetermined") %>% 
# Check data types and make ct a numeric variable not character
  mutate(ct = as.numeric(ct)) %>% 
# Separate the column 'well' into row and column number
  separate(well, into = c("row", "column"), sep = 1, convert = TRUE)

tidy_data %>% 
  glimpse()

# Now we need a primer key to join our data on.
primer_key <- data.frame(
  row = c("A", "B", "C", "D", "E", "F", "G", "H"),
  primer = c(rep("HK", 4), rep("test", 4))
)

# Join the tidy_data and primer key based on 'row'
full_data <- tidy_data %>% 
  left_join(primer_key, by = "row")

full_data %>% 
# We will simulate an image of our well plate
  ggplot(aes(x = factor(column), y = row, label = sample_name, fill = primer)) +
  geom_tile(color = "black") +
  geom_text() +
  theme_bw() +
  labs(title = "qPCR Plate",
       x = "Column",
       y = "Row")

  ggsave("images/qpcr_plate.png",
         width = 10,
         height = 7,
         units = "in",
         dpi = 300)

# Now we begin analysis
# We have 3 samples (two RNAi treatments and a control)
# (3 biological replicates and 3 technical replicates each)
# We need to calculate the mean ct value for technical replicates of each sample

full_data
summarised_data <- full_data %>% 
# Group by sample_name and primer used then calculate mean technical ct
  group_by(sample_name, primer) %>% 
  summarise(avg_ct = mean(ct))

summarised_data %>% 
  ggplot(aes(x = sample_name, y =avg_ct, color = primer)) +
  geom_point() + 
  labs(title = "Mean ct values across replicates",
       y = "Mean Ct value", x = NULL) +
  theme_pubr()

  ggsave(filename = "images/comparing_ct_values_across_bioreplicates.png",
         width = 10,
         height = 5,
         units = "in",
         dpi = 300)


  
# We need to calculate delta ct by subtracting each replicate ct from the rf ct
# We need to extract just samples with 'test' primer and just samples with 'HK'
# We need to get the sample name and replicate number as different columns
  
separated_data <- summarised_data %>% 
    separate(sample_name, into = c("sample", "replicate"), sep = "-")

separated_data

# we filter out the test primers and ref primers to rejoin them for delta ct
test_data <- separated_data %>% 
  filter(primer == "test")

test_data

ref_data <- separated_data %>% 
  filter(primer == "HK") %>% 
  rename(ref_ct = avg_ct)

ref_data
# we rejoin the ref and test data back but this time on sample name and replicate number
rejoined_data <- left_join(ref_data, test_data, by = c("sample", "replicate"))

rejoined_data

# And we then calculate delta ct (Target - housekeeping)
delta_ct <- rejoined_data %>% 
  mutate(delta_ct = avg_ct - ref_ct) %>% 
  select(sample, ref_ct, delta_ct)

delta_ct

delta_ct %>% 
  ggplot(aes(x = factor(sample), y = delta_ct, color = sample)) +
  geom_point() +
  labs(title = "Delta ct values among Samples",
       y = "Delta ct", 
       x = NULL) +
  theme(legend.position = "none") +
  theme_pubr()

ggsave(filename = "images/comparison_delta_ct.png",
       width = 10,
       height = 5,
       units = "in",
       dpi = 300)

# Now to calculate delta delta ct
# But we need to calculate mean delta ct first for control and RNAi samples
# So we group by sample name and summarise
delta_mean <- delta_ct %>% 
  group_by(sample) %>% 
  summarise(mean_delta = mean(delta_ct)) 

delta_mean

# we will pull the control ct value
control <- delta_mean %>% 
  filter(sample == "Control") %>% 
  pull(mean_delta)

control

# Now we subtract RNAi ct from control ct
delta_delta_ct  = delta_mean %>% 
  mutate(delta_delta = mean_delta - control) %>% 
  select(sample, delta_delta)

delta_delta_ct

# Now we make a scatterplot and compare
delta_delta_ct %>% 
  ggplot(aes(x = factor(sample), y = delta_delta)) +
  geom_point() +
  labs(title = "Delta Delta Ct amongst RNAI-1, RNAI-2, and the Control", 
       y = "Delta Delta Ct Value",
       x = NULL) + 
  theme_pubr()

ggsave(filename = "images/Delta_delta_ct_comparison.png",
       width = 10,
       height = 5,
       units = "in",
       dpi = 300)

# Relative Fold Change
#the amount of cDNA theoretically doubles every cycle.
rel_concentration <- delta_delta_ct %>% 
  mutate(rel_conc = 2 ^ - delta_delta)

# Primer efficiency 
rel_concentration %>% 
ggplot(aes(x = sample, y = rel_conc)) +
geom_point() +
theme_pubr()
