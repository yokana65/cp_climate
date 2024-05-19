---
title: "Compositional Data Analysis"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    number_sections: true
---

# Introduction

This is where you can write an introduction to your document.

# Data Analysis

## CHB 14

## KL 15

A first look of the data in KL_15:

```{r, eval=TRUE}
# we want to read them in R
data_KL15 <- read.table("data/Africa_NE_200/data/data_KL15_XRF_tuned.txt", header = TRUE, sep = "\t")

# summarise all attributes of our data
summary(data_KL15)

```

This is clearly a compositional data set, as the sum of all elements is constant. 

# Convert Matlab files


## Read the data

```{r, eval=TRUE}



```