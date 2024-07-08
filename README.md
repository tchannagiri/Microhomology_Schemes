# Microhomology Schemes

## Overview

The purpose of these scripts are to generate schemes for the (1) human plasmid and (2) yeast constructs that are used in the double-strand break repair experiments by the Storici Lab at Georgia Tech ([https://storicilab.gatech.edu/](https://storicilab.gatech.edu/)). This script is developed by the Math-Bio group at the University of South Florida ([https://knot.math.usf.edu/](https://knot.math.usf.edu/)). The publication this respository was made for is at [https://www.biorxiv.org/content/10.1101/2022.11.01.514688v1](https://www.biorxiv.org/content/10.1101/2022.11.01.514688v1).

## Installation/Usage

* Installation: Download this repository. No further installation is needed.

* Usage: Open the file `html/index.html` in a web browser to allow navigating to the different schemes.  Use the dropdowns at the top to change the type of schemes displayed and its visual aesthetics. Use the download button to download the schemes as SVG files.

## Python Scripts

* `make_plasmid_data.py`: Generates JavaScript data files in `html/plasmid` using the CSV files in `csv/plasmid`.
* `make_plasmid_2_data.py`: Generates JavaScript data files in `html/plasmid_2` using the CSV files in `csv/plasmid_2`.
* `make_yeast_data.py`: Generates JavaScript data files in `html/yeast` using the CSV files in `csv/yeast`.

## HTML/Javascript

* `html/FileSaver.js-XXX`: Library used for allowing SVG downloads in the browser ([https://github.com/eligrey/FileSaver.js](https://github.com/eligrey/FileSaver.js)).
* `html/JSZip`: Library used for zipping multiple SVG files for download ([https://stuk.github.io/jszip/](https://stuk.github.io/jszip/)).
* `html/index.html`: Master index that allows navigating to the different schemes.
* `html/schemes.js`: JavaScript code for plotting the schemes based on the JavsScript data files.
* `html/plasmid`, `html/plasmid_2`, `html/yeast`: Files for displaying the specific schemes and for housing the specific data for each construct.

## SVG/PDF

* `saved_svgs`: The SVGs of each construct saved for convenience. The same data that can be obtained with the download buttons.
* `final_pdfs`: The hand-edited SVGs and PDFs used for the supplementary figures of the final publication. These were edited/created using Inkscape ([https://github.com/inkscape/inkscape](https://github.com/inkscape/inkscape)).

## Maintainer

This respository is maintained by Tejasvi Channagiri ([https://github.com/tchannagiri](https://github.com/tchannagiri)).
