.. cfDNAPro documentation master file, created by
   sphinx-quickstart on Wed Aug  5 17:28:19 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

--------------
About cfDNAPro
--------------


Liquid Biopsy Cell-free DNA Fragment Size Profiler

Currently cfDNAPro is only compatible with insert sizes metrics files exported by picard tools:
`CollectInsertSizeMetrics <http://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics>`__. 
More features will be added in coming versions.



Introduction
------------

cfDNA fragment size metrics are important features for utilizing liquid
biopsy in tumor early detection, diagnosis, therapy personlization and
monitoring. Analyzing and visualizing insert size metrics could be time
intensive.

This package intends to simplify this exploration process, and it offers
two sets of functions for data characterization and data visualization.

Installation
------------

.. code:: R

    if (!require(devtools)) install.packages("devtools")
    library(devtools)
    devtools::install_github("hw538/cfDNAPro", build_vignettes = TRUE)

A Quick Example
---------------

If your insert metrics size txt files from multiple cohorts were stored
in their own folders under *path/to/main/folder*, and you would like to
visualize the mode size of each sample, the simpliest way is:

.. code:: R

    library(cfDNAPro)
    path <- "path/to/main/folder"
    myplot <- callMode(path = path) %>% plotMode()

Vignettes
---------

To see the vignettes in Rstudio (you have to indicate
``build_vignettes = TRUE`` during aforementioned installation step), use
the command:

.. code:: R

    browseVignettes("cfDNAPro")

Citation
--------


Please cite package ‘cfDNAPro’ in publications:

Haichao Wang and Christopher G. Smith (2019). cfDNAPro: Liquid Biopsy
Cell-free DNA Fragment Size Profiler. R package version 0.0.4.
https://github.com/hw538/cfDNAPro


.. toctree::
   :caption: Getting started
   :name: getting_started
   :hidden:
   :maxdepth: 1

   tutorial/getting_started


.. toctree::
   :caption: Reference 
   :name: ref
   :hidden:
   :maxdepth: 1