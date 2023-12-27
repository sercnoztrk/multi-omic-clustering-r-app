# Multi-Omic Clustering of TCGA-LUAD Patients

This app utilizes multi omic integration techniques to analyze
- Copy Number Variation (CNV)
- RNA Sequence
- Protein Sequence
- Clinicals

datasets from [TCGA](https://www.cancer.gov/ccg/research/genome-sequencing/tcga), after querying and preprocessing those in required to analyze.

## Querying TCGA Data

Preview of queried datasets

- CNV

<img src="/screenshots/query/Screenshot from 2020-12-05 18-40-29.png" width="49%">  <img src="/screenshots/query/Screenshot from 2020-12-05 18-41-14.png" width="49%">

- RNA

<img src="/screenshots/query/Screenshot from 2020-12-05 18-42-17.png" width="49%">  <img src="/screenshots/query/Screenshot from 2020-12-05 18-42-37.png" width="49%">

- Protein

<img src="/screenshots/query/Screenshot from 2020-12-05 18-43-31.png" width="49%">  <img src="/screenshots/query/Screenshot from 2020-12-05 18-43-57.png" width="49%">

- Clinicals

<img src="/screenshots/query/Screenshot from 2020-12-05 18-45-28.png" width="49%">  <img src="/screenshots/query/Screenshot from 2020-12-05 18-45-56.png" width="49%">

## Clustering

Integrated clustering by individual technique

- Bayesian Consensus Clustering

<img src="/screenshots/clustering/Screenshot from 2020-12-05 18-27-37.png" width="49%">  <img src="/screenshots/clustering/Screenshot from 2020-12-05 18-29-39.png" width="49%">

<img src="/screenshots/clustering/Screenshot from 2020-12-05 18-29-59.png" width="49%">  <img src="/screenshots/clustering/Screenshot from 2020-12-05 18-30-15.png" width="49%">

- iClusterPlus

<img src="/screenshots/clustering/Screenshot from 2020-12-07 16-02-23.png" width="49%">  <img src="/screenshots/clustering/Screenshot from 2020-12-07 16-18-05.png" width="49%">

<img src="/screenshots/clustering/Screenshot from 2020-12-07 16-18-21.png" width="49%">  <img src="/screenshots/clustering/Screenshot from 2020-12-07 16-18-35.png" width="49%">

- LRAcluster

<img src="/screenshots/clustering/Screenshot from 2020-12-07 17-00-10.png" width="49%">  <img src="/screenshots/clustering/Screenshot from 2020-12-07 17-18-45.png" width="49%">

<img src="/screenshots/clustering/Screenshot from 2020-12-07 17-19-03.png" width="49%">  <img src="/screenshots/clustering/Screenshot from 2020-12-07 17-19-22.png" width="49%">

- PINSPlus

<img src="/screenshots/clustering/Screenshot from 2020-12-07 21-34-23.png" width="49%">  <img src="/screenshots/clustering/Screenshot from 2020-12-07 21-35-12.png" width="49%">

<img src="/screenshots/clustering/Screenshot from 2020-12-07 21-35-38.png" width="49%">  <img src="/screenshots/clustering/Screenshot from 2020-12-07 21-35-53.png" width="49%">

- SNF

<img src="/screenshots/clustering/Screenshot from 2020-12-08 16-34-43.png" width="49%">  <img src="/screenshots/clustering/Screenshot from 2020-12-08 16-35-09.png" width="49%">

<img src="/screenshots/clustering/Screenshot from 2020-12-08 16-35-25.png" width="49%">  <img src="/screenshots/clustering/Screenshot from 2020-12-08 16-35-44.png" width="49%">

## Analysis

After clustering patients by integrating CNV, RNA, Protein datasets, every cluster of patients could be analyzed either with:

- Clinical Statistics

<img src="/screenshots/analysis/Clinical Stats.png" width="100%">

- Enrichment Analysis

<img src="/screenshots/analysis/Enrichment Analysis.png" width="100%">

- Gene Expression Heatmap

<img src="/screenshots/analysis/Gene Expression Heatmaps.png" width="100%">

- MAF Summary

<img src="/screenshots/analysis/Maf Summary.png" width="100%">

- Oncoplot

<img src="/screenshots/analysis/Oncoplot.png" width="100%">

- Permutation Test

<img src="/screenshots/analysis/Permutation Test.png" width="100%">

- Recurrent CNV

<img src="/screenshots/analysis/Recurrent CNV.png" width="100%">

- Survival Analysis

<img src="/screenshots/analysis/Survival Analysis.png" width="100%">


