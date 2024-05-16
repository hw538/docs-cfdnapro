*****************************************
Extracting and Visualising cfDNA Features
*****************************************

.. contents:: Table of Contents

Importing BAM Files for Comprehensive Fragmentomic Feature Extraction and Visualization
=======================================================================================

``cfDNAPro`` can also directly process BAM files,
specifically those generated using Illumina sequencing
technology.
The tool extracts paired-end sequencing information applying the following criteria:
(a) the reads are flagged as'proper pair',
(b) the mapping quality score is at least 30,
and (c) the CIGAR string excludes any 'I' (insertions) or 'D' (deletions).
