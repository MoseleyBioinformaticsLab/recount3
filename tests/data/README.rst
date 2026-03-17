===========================
Test data: recount3 mirror
===========================

A local mirror of a small, real subset of the recount3 raw-files tree.
Tests point ``RECOUNT3_URL``.

Mirror root: ``tests/data/recount3_mirror/recount3/``

What is mirrored
================

Three SRA projects (human, ``G026`` annotation):

========== =======================================================
Project    Contents
========== =======================================================
SRP014565  gene_sums, exon_sums, junctions (MM+ID+RR), metadata
DRP004088  gene_sums, exon_sums, junctions (MM+ID+RR), metadata
SRP009615  BigWig for sample SRR387777 only
========== =======================================================

All files were trimmed from real upstream recount3 data.