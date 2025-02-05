.. _cfdnapro_functions:

cfDNAPro Functions and Parameters
=================================

readBam
--------

The ``readBam()`` function reads a BAM file and returns a curated ``GRanges`` object. It processes the BAM file according to various user-defined parameters, such as genome label, strand mode, and sequence names to keep.

**Parameters**:

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

**Returns**:

A curated ``GRanges`` object containing the genomic alignments.

callTrinucleotide
-----------------

The ``callTrinucleotide()`` function processes a GRanges object, summarizing cfDNA fragment information for each target mutation locus. It annotates each mutation locus with the number and type of supporting fragments according to their read-pair overlap status. The median fragment length is also annotated for each read-pair overlap type. Additionally, the function calculates the locus-based consensus mutation by selecting the most frequent mismatch type, with priority given to concordant read-pair mutations (CO_MUT), followed by single-read mutations (SO_MUT). The consensus mutations are used to derive the trinucleotide substitution types (SBS96).

**Parameters**:

- **frag_obj_mut** (GRanges): A ``GRanges`` object containing fragment and mutation data.

**Returns**:

A dataframe with summarized mutational and trinucleotide data.

plotTrinucleotide
-----------------

The ``plotTrinucleotide()`` function processes and plots trinucleotide data. It first applies specified filters and transformations to the data, then generates a visual representation of the results. The function handles data normalization, exclusion, and retention based on provided column names, and it creates detailed plots with options for customization of plot aesthetics.

**Parameters**:

- **trinuc_df** (DataFrame): The dataframe containing trinucleotide data.
- **exclude_if_type_present** (vector): Mutation locus read-pair overlap types (e.g., CO_MUT, SO_MUT, CO_REF, SO_REF, DO, SO_OTHER, CO_OTHER) whose non-zero presence triggers exclusion of loci. For example, `c("DO")` will exclude any loci that contain even a single discordant read-pair overlap.
- **retain_if_type_present** (vector): Mutation locus read-pair overlap types that must be present to retain those loci. For example, `c("CO")` will retain loci that contain even a single concordant read-pair overlap.
- **remove_type** (vector): Mutation locus read-pair overlap types (e.g., SO_MUT, CO_REF, DO) to set to 0 across all loci in the dataframe.
- **normalize_counts** (bool): If `TRUE`, normalizes SBS counts so they sum to 1. Default is `TRUE`.
- **show_overlap_type** (bool): If `TRUE`, displays read-pair overlap types. Default is `TRUE`.
- **ylim** (numeric): Limits for the y-axis in the plot. Default is `c(0, 0.5)`.
- **plot_title** (str): The title for the plot. Default is `"Trinucleotide Profile"`.
- **y_axis_title** (str): The title for the y-axis. Default is `"Percentage of Single Base Substitutions"`.
- **draw_x_axis_labels** (bool): Whether to draw x-axis labels. Default is `TRUE`.
- **draw_y_axis_labels** (bool): Whether to draw y-axis labels. Default is `TRUE`.
- **draw_y_axis_title** (bool): Whether to display a title for the y-axis. Default is `TRUE`.
- **output_file** (str): The path to the output PDF file. Default is `"./trinucleotide_profile.pdf"`.
- **ggsave_params** (list): A list of parameters to be passed to ``ggplot2::ggsave()``. This list can include any arguments accepted by ``ggsave()``. Example: `list(width = 17, height = 6, units = "cm", device = "pdf")`.

**Returns**:

A trinucleotide SBS plot object and an optional PDF file.
