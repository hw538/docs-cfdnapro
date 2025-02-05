.. _cfdnapro_functions:

cfDNAPro Functions and Parameters
=================================

readBam
========

The ``readBam()`` function reads a BAM file and returns a curated ``GRanges`` object. It processes the BAM file according to various user-defined parameters, such as genome label, strand mode, and sequence names to keep.

### Parameters

- **genome_label** (str): Specifies the genome used in the alignment. Accepted values are `"hg19"`, `"hg38"`, or `"hg38-NCBI"`. Default is `"hg19"`. It loads corresponding genome packages for Homo sapiens sequences.
- **bamfile** (str): The path to the BAM file.
- **curate_start_and_end** (bool): If `TRUE`, the start and end coordinates of alignments are curated. Default is `TRUE`.
- **outdir** (str or NA): Path to save the RDS file. If `NA`, no file is saved.
- **strand_mode** (int): Defines strand mode; 1 means the strand of the pair is taken from the first alignment. Default is `1`.
- **chromosome_to_keep** (vector of str or bool): A character vector containing seqnames to retain in the ``GRanges`` object. Default is `paste0("chr", 1:22)`. If `FALSE`, no filtering is applied.
- **use_names** (bool): Whether to assign read names to the ``GRanges`` object. Default is `TRUE`.
- **galp_flag** (Rsamtools ``ScanBamFlag``): Specifies the flags for scanning the BAM file.
- **galp_what** (vector of str): The fields to return from the BAM file, such as "cigar", "mapq", "isize", "seq", and "qual".
- **galp_tag** (vector of str): Specifies optional fields (tags) to retrieve from the BAM file.
- **galp_mapqFilter** (int): The minimum mapping quality to include a read. Default is `40`.
- **galp_bqFilter** (int): The minimum base quality at the mutation locus. Default is `20`.
- **mutation_file** (str or NULL): An optional file containing mutation loci for mutational annotation.
- **mut_fragments_only** (bool): If `TRUE`, only retrieves alignments overlapping mutation loci. Default is `FALSE`.
- **...**: Additional arguments passed to or from other methods.

### Returns

A curated ``GRanges`` object containing the genomic alignments.

---
