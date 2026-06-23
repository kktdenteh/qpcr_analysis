\# qPCR Delta-Delta Ct Analysis Pipeline



This repository contains an automated R pipeline for processing raw quantitative PCR (qPCR) data and calculating relative gene expression using the comparative $2^{-\\Delta\\Delta C\_t}$ method. 



As a final-year Biochemistry student learning data analytics, I built this workflow to bridge wet-lab laboratory concepts with programmatic data cleaning and analysis.



\## Workflow Overview



The script (`your\_script\_name.R`) executes a complete data pipeline using the `tidyverse` and `janitor` packages:



1\. \*\*Data Ingestion \& Cleaning:\*\* Skips raw machine metadata, standardizes column formats, filters out No-Template Controls (NTC) and "Undetermined" values, and splits well coordinates into rows and columns.

2\. \*\*Plate Mapping:\*\* Dynamically maps target and reference (`HK`) primers across the physical layout using a data-frame join.

3\. \*\*Biostatistical Calculation:\*\* 

&#x20;  \* Averages technical replicates for each sample and primer.

&#x20;  \* Calculates $\\Delta C\_t$ using standard convention: $C\_t(\\text{Target}) - C\_t(\\text{Housekeeping})$.

&#x20;  \* Normalizes against the control sample to find $\\Delta\\Delta C\_t$.

&#x20;  \* Computes final relative concentration ($2^{-\\Delta\\Delta C\_t}$).

4\. \*\*Data Visualization:\*\* Generates and saves clean `ggplot2` plots tracking individual steps of the analysis.



\## Generated Visualizations



The pipeline automatically processes data and saves the following assets to the `images/` directory:

\* `qpcr\_plate.png` – A simulated tile map representing the layout of the well plate.

\* `comparing\_ct\_values\_across\_bioreplicates.png` – Raw $C\_t$ value distributions.

\* `comparison\_delta\_ct.png` – Calculated $\\Delta C\_t$ values across treatments.

\* `Delta\_delta\_ct\_comparison.png` – Normalized $\\Delta\\Delta C\_t$ values against the control.



\## Dependencies



To run this pipeline, you will need the following libraries installed in your R environment:

```r

install.packages("tidyverse")

install.packages("ggpubr")

install.packages("janitor")

```



\## How to Run

1\. Clone this repository: `git clone https://github.com`

2\. Place your raw file in `data/qpcr\_data.csv`.

3\. Run the analysis script in RStudio.



