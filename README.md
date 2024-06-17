# Python-MetaMaps Identity Plotter
This is a script that generates a PDF full of plots from output data produced by MetaMaps.

## Description
A multipage PDF will be generated using the MetaMaps output files ending with '.EM.contigCoverage' and '.EM.lengthAndIdentitiesPerMappingUnit'.
The first page will end in '.identitiesAndCoverage.pdf', and will be written to the same directory that was specified for the previous two files.
A second pdf page (containing organisms flagged as outliers) ending in '.outliers.identitiesAndCoverage.pdf' will be written there as well.

* Overview of plots: https://docs.google.com/document/d/1IQLHs_48OcQaLuEPkGXnO1k3fYRw_zdLG7iZNPRMCno/edit?usp=sharing
```
usage: plot_identities_em.py [-h] [-f MIN_FREQUENCY] [-t MIN_TRIM_MEAN] [-p TRIM_PROPORTION] [-s] [-v] classification_file_prefix

Plot MetaMaps Identity Results

positional arguments:
  classification_file_prefix
                        Prefix of the classification file

options:
  -h, --help            show this help message and exit
  -f MIN_FREQUENCY, --min-frequency MIN_FREQUENCY
                        Minimum frequency of taxon label to plot (default: 0.0)
  -t MIN_TRIM_MEAN, --min-trim-mean MIN_TRIM_MEAN
                        Minimum coverage trim mean value to consider a taxon ID as an outlier. Outliers will be written to a separate PDF where trim_mean <= min_trim_mean (default: 0.0)
  -p TRIM_PROPORTION, --trim-proportion TRIM_PROPORTION
                        Proportion of sorted coverage data to trim from both ends for outlier detection (default: 0.03)
  -s, --sort            An additional sorted PDF containing both outlier and non-outlier data will be written. Takes longer to generate. (default: False)
  -v, --verbose         Verbose mode (default: False)
```

Data filtering is done in two steps.
* Frequency Filtering (Using mapping unit occurrence frequency)
  *  First, any mapping unit or contig that does not occur frequently enough will not be plotted (except for on the genome-wide coverage over all contigs plot)
* Trim Mean Filtering (Using coverages)
  * Next, outlier organism detections will be found by sorting the data, cutting out a percentage/proportion as specified, and then flagging all organisms whose coverage trim-mean is less than or equal to the specified value

## Dependencies

* Python 3
* Matplotlib
* Pandas
* Scipy
* Numpy


## Authors

Kamran Haq

## Acknowledgments

This script is based on plotIdentities_EM.R by the MetaMaps authors at https://github.com/DiltheyLab/MetaMaps
