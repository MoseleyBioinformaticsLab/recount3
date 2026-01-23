Tutorial
========

This tutorial walks through common use cases for the recount3 package.

Downloading Annotation Files
----------------------------

.. code-block:: python

   from recount3 import R3Annotations, R3Resource

   # Human gene annotations
   annotations = R3Annotations(
       organism="human",
       genomic_unit="gene", 
       annotation_file_extension="G026"
   )
   resource = R3Resource(annotations)
   resource.download(cache_mode="enable")
   df = resource.load()

   print(f"Loaded {len(df)} annotation records")

Working with Count Data
-----------------------

.. code-block:: python

   from recount3 import search_count_files_gene_or_exon, R3ResourceBundle

   # Find and download count files for a project
   resources = search_count_files_gene_or_exon(
       organism="human",
       data_source="sra",
       genomic_unit="gene",
       project="SRP107565",
       annotation_file_extension="G026"
   )

   bundle = R3ResourceBundle(resources)
   bundle.load()
   
   # Stack count matrices
   combined_df = bundle.stack_count_matrices()

Using BigWig Files
------------------

.. code-block:: python

   from recount3 import search_bigwig_files

   # Get BigWig coverage files
   bw_resources = search_bigwig_files(
       organism="mouse",
       data_source="sra", 
       project="DRP001299",
       sample="DRR014697"
   )

   for resource in bw_resources:
       bw_file = resource.load()
       # Access coverage data
       coverage = bw_file.stats("chr1", 1000, 2000)