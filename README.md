# Fryer et al., 2021


## Overview

These are the scripts used for preparing plate layouts, creating some figures
and analysing certain data. 
They're provided along with the raw data so that others can follow our 
workflow in the interests of transparency. These scripts were originally
built and used within a personal file system, so some parts have been edited for
anonymisation and clarification purposes, and we believe they should now run 
as intended (and functionally as we used them ourselves) from any file location.


## Benzylguanine quantification

The scripts in this folder were used to convert raw HPLC data into a
chromatogram figure ('HPLC-chromatogram-figure.R') and for quantitative 
analysis using peak integrations output from Agilent ChemStation
('Quantifying BG formation efficiency.R'). We provide the data used in
the script itself and in a CSV file, formatted for clarity. 


## Sample randomiser

This script produces randomised 96-well plate layouts, and was used during 
plate setup for the apoptosis work in our paper. It actually takes user 
input at runtime, but we've refactored these parts of the code so it should 
run as when we used it, but without requiring additional user input.


## Apoptosis FC analysis

The script here is responsible for handling statistics which were output from 
the Attune NxT (gates were set manually using control samples), deconvoluting
them using the randomised plate layouts generated earlier and then plotting. 