***********************************************
Multi-Feature Integration and Visualisation 
***********************************************

.. contents:: Table of Contents

------------------------------------------

Each type of cfDNA fragment plot (length, motif, and CNV)
can also incorporate mutational data by setting ``integrate_mut = TRUE``.
This integration is feasible if the GRanges object was created through
the execution of ``readBam(mutation_file = /path/to/file.tsv)`` or
``readBam(call_mutations = TRUE)``.

The default behavior of the ``readBam()`` function is to process all DNA fragments
within the BAM file. However, by setting ``mut_fragments_only = TRUE``,
the function will only analyse fragments that overlap with the specified mutation loci.
This option reduces the computational load and may be adequate for users who do not
require information from non-overlapping fragments for motif and length plots.

.. code:: R

  # Process all fragments present within the BAM file without mutational annotation
  readBam(bamfile = "path/to/bamfile.bam")

  # Process all fragments present within the BAM file with additional mutation annotation
  readBam(bamfile = "path/to/bamfile.bam",
          mutation_file = "/path/to/mutation_file.tsv")
  
  # Process fragments that overlap loci indicated in the mutation file
  readBam(bamfile = "path/to/bamfile.bam",
          mutation_file = "/path/to/mutation_file.tsv",
          mut_fragments_only = TRUE)
  
  # Process fragments that overlap loci generated during the pileup
  readBam(bamfile = "path/to/bamfile.bam",
          call_mutations = TRUE,
          mut_fragments_only = TRUE)
  

Fragment Length and Fragment Mutations
================================================

.. code:: R

  # First, call the cfDNA fragment length object
  length_obj <- callLength(gr_ob, integrate_mut = TRUE, mut_plot_type = "normalised")

  # Plot normalised MUT and REF fragment lengths
  plotLength(length_obj,
             output_file = "/path/to/plot.png",
             ggsave_params = list())

The ``mut_plot_type`` parameter allows the user to display the normalised counts of fragments as fractions when using ``mut_plot_type = normalised``.
The default and the normalised mutation fragment plots use the overlapping fragments supporting the references base as REF.
If ``mut_plot_type = ref_outer`` or ``mut_plot_type = ref_outer_normalised`` then REF fragements that do not overlap the mutation locus will be used in plotting.
However, this requires that the user generates the initial GRanges object with mut_only = FALSE so that outer fragment are included. 


Plots will be automatically saved if you specify a path and
filename using the ``output_file`` argument.
Additionally, users can customize the plot file using the ``ggsave_params`` parameter.


.. image:: static/cfDNA_plasma_length_mut.png
  :width: 500
  :height: 440
  :align: center
  :alt: mut_length

|

Fragment Motif and Fragment Mutations
================================================


The motifs can be plotter similarly.

.. code:: R

  # First, call the cfDNA fragment length object
  motif_obj <- callLength(gr_obj,
                          integrate_mut = TRUE
                          mut_plot_type = "normalised")

  # Plot normalised MUT and REF fragment lengths
  plotMotif(motif_obj,
            output_file = "/path/to/plot.png",
            ggsave_params = list())


.. image:: static/cfDNA_plasma_motif_mut.png
  :width: 800
  :height: 200
  :align: center
  :alt: mut_motif

|

Fragment Copy-Number and Fragment Mutations
================================================

You can also plot CNV with integrated mutational information for each SNV within genes of interest.
This requires that the gene of interest includes SNVs listed in the mutation file or those processed during the pileup in the ``readBam()`` function.
The plot will then display total counts of all SNVs within that gene, including both MUT and REF fragments, as an additional annotation for the specified genes.

.. code:: R

  # First, call the cfDNA fragment length object
  legth_obj <- callLength(gr_obj, genome_label = "hg38")

  # Process fragments that overlap loci indicated in the mutation file
  frag_obj <- readBam(bamfile = "path/to/bamfile.bam",
                      mutation_file = "/path/to/mutation_file.tsv",
                      mut_fragments_only = TRUE) 

  # Plot normalised MUT and REF fragment lengths
  plotCNV(motif_obj,
          frag_obj_mut =  frag_obj
          output_file = "/path/to/plot.png",
          gene_to_highlight = list("ENTREZID" = NULL,
                                                "ENSEMBL" = NULL,
                                                "SYMBOL" = c("BRAF", "PIK3CA")),
          ggsave_params = list())


.. image:: static/cfDNA_plasma_cnv_mut.png
  :width: 800
  :height: 180
  :align: center
  :alt: mut_cnv

|
