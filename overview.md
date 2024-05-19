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

A first look of the data in KL_15:

```{r, eval=TRUE}
# we want to read them in R
data_KL15 <- read.table("data/Africa_NE_200/data/data_KL15_XRF_tuned.txt", header = TRUE, sep = "\t")

# summarise all attributes of our data
summary(data_KL15)

```