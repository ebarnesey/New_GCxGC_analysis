---
title: "Investigating Sources of Organic Aerosol in the Amazon"
author: "Emily Barnes"
output: html_document
---
# Final Project Proposal

##*Project Desctiption*: 
During the 'Wet Season' of early 2014 and the 'Dry Season' of mid 2015, a team of researchers from the Goldstein Lab traveled to a semi-remote site in the center of the Amazon Rainforest to participate in the GoAmazon Field Campaign.  Throughout both of these visits, researchers collected fine particulate in the air on quartz fiber filters, which were subsequently frozen and transported back to UC Berkeley.  In late 2018, I processed the majority of these filters on a *TD-GCxGC-EI-HRToF-MS*.  This instrument thermally desorbs organic components of the trapped aerosols (ranging in mass from C12-C36 n-alkane equivalents), separates them in two dimensions by both volatility and polarity, and generates a mass spectra for each separated species using an EI-HR-ToF mass spec. While a few of the separated compounds may be difinitively identified by matching their mass spectra to those of known standards processed on similar instruments, the vast majority have no authentic standards.  My challenge is to glean as much information as possible about the compounds that cannot be definitively identified, including VOC source group, chemical functional group, and sensitivity to anthropogenic influence. My ultimate goal is to better understand the differences between wet season and dry season amazon, and to contribute to answering the question of why dry season aerosol loadings are so much higher that wet season loadings.

## *Data Description*:
The output of the *TD-GCxGC-EI-HRToF-MS* is an .h5 file which contains the location (based on retention times in the two columns), response (known as 'blob volume'), and mass spectrum of each compound within the instrument's detection limitations.  Using a software called 'GC Image', I search all of these spectra against three libraries of mass spectra that I have generated from previous images.  These three libraries are 1) a library of the 27 internal standard compounds that I injected on each filter, 2) a library of all of the contamination comoounds identified by running field blanks using the full preparation protocol, and 3) a library of 1,000 + amazonian aerosol compounds identified in key periods of differing climatic conditions throughout the campaign.  My ultimate goal is to generate robust time series of each of these compounds which will allow me to identify its sources and sensitivity to differing oxidative conditions.  The output of the GC-Image analysis is a CSV file containing the following key pieces of information: 

* Compound Name (as assigned by library match)
* Compound Library (name of library from which identity was determined)
* Compound Linear Retention Index
    + position of the compound in the first dimension (volatility) relative to an alkane standard series
* Linear Retention Index of Library Entry
    + position of the **library** compound that was identified as a match to the current compound.  This is critical, as compounds are often misidentified based on solely mass spectrum, and ensuring that Library LRI and Compound LRI are close vastly improves the odds of making a correct match
* Library Match Factor
    + A measure of how closely the mass spectrum of the image compound matches that of the compound saved in the library
* Volume
    + A measure of response
  
I ran approximately 300 filters and will ultimately generate 300 CSV's, each containing ~1000-2500 compounds.  My data source will be a selected subset of these CSV's, along with a CSV tying each data file to the time the filter was collected and how many segments of filter were run.    

## *Questions/Analysis Tasks*
1. Determine which internal standard species are most reliable and best suited for normalization of analyte response

* The internal standard is a mixture of compounds spanning a representative spectrum of volatilities and polarities added to each filter such that the same amount of each compound is present in each instrument run.  The purpose of the internal standard is to correct the  percieved blob volume for something closer to mass of a particular analyte run, as the instrument is more sensitive to some areas of the 2D GC space than others, instrument condition changes over the course of an intensive run campaign, and high filter loading leads to better recovery of all species (a phenomenon known as *matrix effects*).  However, some of the internal standard species are frequently misidentified, frequently coellute with analyte blobs, or are for other reasons poor options for normalization.  I plan to assess which internal standards are appropriate choices for use in correcting analyte blob volumes by automating a process for importing files, correctly identifying internal standards, visualizing timelines of those internal standards, and assessing the the degree to which they correlate with other internal standards throughout the campaign.  

2. Automate normalization of analyte response by internal standard response 
* Once I have selected an appropriate subset of internal standards, I need to build a framework to normalize analyte blob volumes by one or more of the nearest internal standards.  This will be critical to achieving my goal of creating a time series for each of my 1000+ amazonian analytes of interest.  

3. Determine how to normalize response for different amounts of filter material analyzed
    
* In order to keep the majority of my analytes above detection limits for clean days and below an overload threshold for polluted days, I was forced to run different amounts of filter for different days.  Ideally, amount of compound run (eg. number of filter punches) and response should be linearly correlated- however, that does not currently appear to be the case.  To enable me to compare images generated by different numbers of filter segments, it is critical that I determine how to properly normalize response (indicated by 'blob volume') by number of punches for different compounds.  

## *Skills*:

* 5+ dplyr verbs
* Writing and using custom R functions
* Use of regular expressions



