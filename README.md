# CRARs-2024-Metagenomics
This is a collection of scripts for generating reports on metagenomics data produced by the analysis of the MetaMaps and MTSV mapping algorithms.

**Prior to running the diversity scripts alpha_diversity_report.py and beta_diversity_report.py, you need to generate a pickle file**, see pickle_filtered_reads.py

## pickle_filtered_reads.py
```
usage: pickle_filtered_reads.py [-h] [-mtsv MTSV_FILE] [-mtsvl MTSV_LOOKUP_FILE] [-meta META_MAPS_FILE] [-metaref META_MAPS_REFERENCE_FILE] [-mtsvref MTSV_REFERENCE_FILE] [-s SEED]
                                [-f MIN_FREQUENCY] [-a ALIAS] [-C] [-B] [-r RARE] [-S]

Takes either MTSV or MetaMaps reads and loads into a pickle file for analysis.

options:
  -h, --help            show this help message and exit
  -mtsv MTSV_FILE, --mtsv-file MTSV_FILE
                        MTSV File (default: None)
  -mtsvl MTSV_LOOKUP_FILE, --mtsv-lookup-file MTSV_LOOKUP_FILE
                        MTSV Lookup File (default: None)
  -meta META_MAPS_FILE, --meta-maps-file META_MAPS_FILE
                        MetaMaps File (default: None)
  -metaref META_MAPS_REFERENCE_FILE, --meta-maps-reference-file META_MAPS_REFERENCE_FILE
                        MetaMaps file to filter MTSV with (default: None)
  -mtsvref MTSV_REFERENCE_FILE, --mtsv-reference-file MTSV_REFERENCE_FILE
                        MTSV file to filter MetaMaps reads with (default: None)
  -s SEED, --seed SEED  Seed for random number generator (default: 0)
  -f MIN_FREQUENCY, --min-frequency MIN_FREQUENCY
                        Minimum frequency of taxon label to plot (default: 0.0)
  -a ALIAS, --alias ALIAS
                        Alias for file (default: None)
  -C, --clear           DANGER! Clears all files in reads directory (default: False)
  -B, --sig-bin         Enable sig-bin filtering (METAMAPS ONLY) (default: False)
  -r RARE, --rare RARE  Rareify to given read count (default: None)
  -S, --SRS             Enable SRS. Requires R to be set up in the current system's environment (default: False)
```

## plot_identities_em.py
* Overview of plots: https://docs.google.com/document/d/1IQLHs_48OcQaLuEPkGXnO1k3fYRw_zdLG7iZNPRMCno/edit?usp=sharing
```
usage: plot_identities_em.py [-h] [-f MIN_FREQUENCY] [-t MIN_TRIM_MEAN] [-p TRIM_PROPORTION] [-o OUTPUT] [-I] [-S] [-v] classification_file_prefix

Plot MetaMaps Identity Results

positional arguments:
  classification_file_prefix
                        Prefix of the classification file

options:
  -h, --help            show this help message and exit
  -f MIN_FREQUENCY, --min-frequency MIN_FREQUENCY
                        Minimum frequency of taxon label to plot (default: 0.0)
  -o OUTPUT, --output OUTPUT
                        Output file name (default: )
  -I, --ignore-ids      Ignore ids in the file .ignoreids (default: False)
  -S, --skip-coverage-filter
                        Skip the coverage filtering step. Saves time if you already have a .ignoreids file. Will not generate an outlier pdf. (default: False)
  -v, --verbose         Verbose mode (default: False)
```

Data filtering is done in two steps.
* Frequency Filtering (Using mapping unit occurrence frequency)
  *  First, any mapping unit or contig that does not occur frequently enough will not be plotted (except for on the genome-wide coverage over all contigs plot)
* SigBin Filtering (Using coverages)
  * Next, outlier organism detections will by identifying organisms with few concentrated histogram bins

## Dependencies

* Python 3
* Matplotlib
* Pandas
* Scipy
* Numpy


## Authors

Kamran Haq

## Acknowledgments

The script plot_identities_em.py is based on plotIdentities_EM.R by the MetaMaps authors at https://github.com/DiltheyLab/MetaMaps
