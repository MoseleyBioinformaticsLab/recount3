==========================
Test data: recount3 mirror
==========================

A local mirror of a small, real subset of the recount3 raw-files tree.
Tests set ``RECOUNT3_URL`` to the ``file://`` URL for the mirror root at
``tests/data/recount3_mirror/recount3/``.

What is mirrored
================

The mirror includes three human SRA projects. The gene and exon counts use
the ``G026`` annotation, and the corresponding annotation files are included.
The mirror also contains project discovery files for human and mouse SRA.

============= ============================================================
Project       Contents
============= ============================================================
``DRP004088`` ``exon_sums``, ``gene_sums``, junctions (ID, MM, RR), metadata
``SRP009615`` BigWig for sample ``SRR387777`` only
``SRP014565`` ``exon_sums``, ``gene_sums``, junctions (ID, MM, RR), metadata
============= ============================================================

All files were trimmed from real upstream recount3 data.
