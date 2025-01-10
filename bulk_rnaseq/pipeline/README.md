**01. download_genome.pbs**                    : download genome annotation (GRCh38)

**02. quality_check.pbs**                      : quality check of sequencing reads using FastQC

**03. trimming.pbs**                           : trimming of sequencing reads using BBDuk suite

**04. quality_check_post_trimming.pbs**        : quality check after the trimming step

**05. adapter_removal.pbs**                    : removal of adapter sequences

**06. quality_check_post_adapter.pbs**         : quality check after the adapter removal step 

**07. generate_genome_indexes.pbs**            : building a STAR index file/create a genome index


**infer_library_strandness.pbs**: verify the type of sequencing library considerinf infer_experiment.py (avialable at https://rseqc.sourceforge.net/#infer-experiment-py)
