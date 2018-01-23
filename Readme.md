# SNP_Density_windows
Git repository created 10-28-16
Bradford Condon, University of Kentucky

## Overview
Project currently under construction.

The goal of this project is to generate sliding window reports for SNPs based on clade, and propose hypotheses regarding SNP patterns.



# Perl

## CallWindowClades.pm

This is the main script for generating sliding windows.


Required parameters:

* `-d|directory`  
* cladelist
	- a file that matches each strain to a list.
	- use `example_input_files/masterCladeList.txt`
* `e|lengthFile`
	- list of scaffold lengths.
	- see `example_input_files/B71_SZL_lengthlist_named.txt`.

## File listing:

### R

* helperFunctions.R
	* trimFileNamesFromSNPTables
	* readSNPReports
	* slidingWindowForSNP
	* slidingWindowForSNPoutputDF
		* debugging
	* analyzeSlidingWindows
		* still under construction

* LoadData.R
	* Script reads in SNP reports.  Designed to work locally


### ShinyDeploy

This directory contains everything needed to open a Docker container running a Shiny webapp to display plots.  To operate, replace sample data with output of these scripts use the Docker-compose tool

```
Docker-compose up
```