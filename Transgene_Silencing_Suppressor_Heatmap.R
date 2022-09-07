library(tidyverse)
library(zoo)
library(heatmaply)
library(plotly)
library(htmlwidgets)

raw_data <- read.csv("Transgene_Silencing_Suppressor_Chart.csv") %>%
  as_tibble()

raw_data

## create function to parse strings into usable values
revalue <- function(string) {
  ifelse(str_detect(string, "\\+"), 1,
          ifelse(str_detect(string, "-"), 0, NA))
}

## clean data
data <- raw_data %>%
  ## remove empty rows and populate Class column
  mutate(Class = ifelse(Class == "", NA, Class)) %>%
  filter(row_number() <= 191) %>%
  na.locf() %>%
  filter(Reference != "") %>%
  ## split columns with two values
  mutate("Exogenous Germline Transgene Silencing" = 
           ifelse(
             str_detect(Transgene.Germline.Desilencing, "Ex"),
             Transgene.Germline.Desilencing,
             NA),
         "Endogenous Germline Transgene Silencing" =
           ifelse(
             str_detect(Transgene.Germline.Desilencing, "Is"),
             Transgene.Germline.Desilencing,
             NA),
         "Enhanced Germline RNAi" =
           ifelse(
             str_detect(Enhanced.RNAi.phenotype..G...germline..S...somatic., "G"),
             Enhanced.RNAi.phenotype..G...germline..S...somatic.,
             NA),
         "Enhanced Somatic RNAi" =
           ifelse(
             str_detect(Enhanced.RNAi.phenotype..G...germline..S...somatic., "S"),
             Enhanced.RNAi.phenotype..G...germline..S...somatic.,
             NA),
         "Germline RNAi Defective" =
           ifelse(
             str_detect(RDE.phenotype..G...germline..S...somatic., "S"),
             RDE.phenotype..G...germline..S...somatic.,
             NA),
         "Somatic RNAi Defective" =
           ifelse(
             str_detect(RDE.phenotype..G...germline..S...somatic., "S"),
             RDE.phenotype..G...germline..S...somatic.,
             NA)) %>%
  ## remove unwanted rows and reorder chart
  select(Mutant..or.if.RNAi...predicted.gene.,
         Transgene.Somatic.Silencing,
         "Exogenous Germline Transgene Silencing",
         "Endogenous Germline Transgene Silencing",
         Suppressor.of.Transgene.Silencing.in.an.eri.1.Background,
         Suppressor.of.Transgene.Silencing.in.an.eri.6.7.Background,
         Suppressor.of.Transgene.Silencing.in.an.rrf.3.Background,
         Suppressor.of.Transgene.Silencing.in.a.tam.1..cc567..Background,
         "Enhanced Somatic RNAi",
         "Enhanced Germline RNAi",
         "Somatic RNAi Defective",
         "Germline RNAi Defective",
         Class,
         Description,
         Reference) %>%
  rename(Gene = Mutant..or.if.RNAi...predicted.gene.,
         "Somatic Transgene Silencing" = Transgene.Somatic.Silencing,
         "Suppressor of Transgene Silencing (eri-1)" = Suppressor.of.Transgene.Silencing.in.an.eri.1.Background,
         "Suppressor of Transgene Silencing (eri-6/7)" = Suppressor.of.Transgene.Silencing.in.an.eri.6.7.Background,
         "Suppressor of Transgene Silencing (rrf-3)" = Suppressor.of.Transgene.Silencing.in.an.rrf.3.Background,
         "Suppressor of Transgene Silencing (tam-1)" = Suppressor.of.Transgene.Silencing.in.a.tam.1..cc567..Background
         ) %>%
  mutate(across(
    "Somatic Transgene Silencing":"Germline RNAi Defective",
    revalue))
  
data

##check duplicates as we will be using the Gene column as a row name
sum(duplicated(data$Gene))
##extract duplicates

duplicated_data <- data %>%
  filter(duplicated(Gene)|duplicated(Gene, fromLast = TRUE))

duplicated_data

##we can remove the second rde but dpy-20 must be merged
##mutate replace had some issues on somatic transgene silencing column
##I used base r to circumvent these issues

fixed_data <- data[-c(23, 107),]
fixed_data[111, 2] = 1
fixed_data[111, 15] = "1, 3"

fixed_data

##plot data as a heatmap
heatmap_data <- fixed_data %>%
  select(-Class, -Description, -Reference) %>%
  column_to_rownames("Gene") %>%
  data.matrix()

heatmaply(heatmap_data, 
          Rowv = FALSE, 
          Colv = FALSE,
          showticklabels = c(TRUE, FALSE),
          ylab = "Gene",
          fontsize_col = 7)

##create custom hovertext dataframe keeping Class column to sort
label_data <- fixed_data %>%
  mutate(across(2:12,
                ~ paste("Results:", ifelse(.x == 1, "Positive",
                                                 ifelse(.x == 0, "Negative", 
                                                        "Experiment not done")),
                        "\nDescription:", .data$Description,
                        "\nReferences:", .data$Reference)))

label_data
##too big to make sense of and options in heatmaply are limited
##to better visualize subset data by Class (excluding uninteresting classes)

##heatmap of RNAi factors
rnaifactor_heatmap_data <- fixed_data %>%
  filter(Class == "RNAi Factors") %>%
  arrange(Gene) %>%
  select(-Class, -Description, -Reference) %>%
  column_to_rownames("Gene") %>%
  data.matrix()

rnaifactor_label_data <- label_data %>%
  filter(Class == "RNAi Factors") %>%
  arrange(Gene) %>%
  select(-Class, -Description, -Reference) %>%
  column_to_rownames("Gene")

rnaifactor_heatmap <- heatmaply(rnaifactor_heatmap_data,
          Rowv = FALSE,
          Colv = FALSE,
          hide_colorbar = TRUE,
          margins = c(100,50,NA,0),
          fontsize_row = 7,
          fontsize_col = 7,
          label_names = c("Gene", "Phenotype", "Coded Value"),
          custom_hovertext = rnaifactor_label_data,
          main = "RNAi Factors")

rnaifactor_heatmap

##heatmap of chromatin factors
chromatinfactor_heatmap_data <- fixed_data %>%
  filter(Class == "Chromatin Factors") %>%
  arrange(Gene) %>%
  select(-Class, -Description, -Reference) %>%
  column_to_rownames("Gene") %>%
  data.matrix()

chromatinfactor_label_data <- label_data %>%
  filter(Class == "Chromatin Factors") %>%
  arrange(Gene) %>%
  select(-Class, -Description, -Reference) %>%
  column_to_rownames("Gene")

chromatinfactor_heatmap <- heatmaply(chromatinfactor_heatmap_data,
          Rowv = FALSE,
          Colv = FALSE,
          hide_colorbar = TRUE,
          margins = c(100,50,NA,0),
          fontsize_row = 7,
          fontsize_col = 7,
          label_names = c("Gene", "Phenotype", "Coded Value"),
          custom_hovertext = chromatinfactor_label_data,
          main = "Chromatin Factors")

chromatinfactor_heatmap

##heatmap of RNA binding factors
rnabindingprocessing_heatmap_data <- fixed_data %>%
  filter(Class == "RNA Binding and Processing") %>%
  arrange(Gene) %>%
  select(-Class, -Description, -Reference) %>%
  column_to_rownames("Gene") %>%
  data.matrix()

rnabindingprocessing_label_data <- label_data %>%
  filter(Class == "RNA Binding and Processing") %>%
  arrange(Gene) %>%
  select(-Class, -Description, -Reference) %>%
  column_to_rownames("Gene")

rnabindingprocessing_heatmap <- heatmaply(rnabindingprocessing_heatmap_data,
          Rowv = FALSE,
          Colv = FALSE,
          hide_colorbar = TRUE,
          margins = c(100,50,NA,0),
          fontsize_row = 7,
          fontsize_col = 7,
          label_names = c("Gene", "Phenotype", "Coded Value"),
          custom_hovertext = rnabindingprocessing_label_data,
          main = "RNA Binding and Processing Factors")

rnabindingprocessing_heatmap

##bind relevent heatmaps together
plots <- subplot(rnaifactor_heatmap, chromatinfactor_heatmap, rnabindingprocessing_heatmap,
                 margin = 0.1,
                 widths = c(0.29,0.4,0.29)) %>%
  layout(title = NA,
         annotations = list(
           list(x = 0.10,  
                y = 1.0,  
                text = "RNAi Factors",  
                xref = "paper",  
                yref = "paper",   
                xanchor = "center",  
                yanchor = "bottom",  
                showarrow = FALSE),  
           list(x = 0.50,  
                y = 1,  
                text = "Chromatin Factors",  
                xref = "paper",  
                yref = "paper",  
                xanchor = "center",  
                yanchor = "bottom",  
                showarrow = FALSE),  
           list(x = 0.90,  
                y = 1,  
                text = "RNA Factors",  
                xref = "paper",  
                yref = "paper",  
                xanchor = "center",  
                yanchor = "bottom",  
                showarrow = FALSE)))

plots

saveWidget(plots,
           file = "Transgene_Suppressor_Chart_Widget.html")
