RNA-sequencing of SUDHL-5 cell line treated with  three RiBi inhibitors (Doxorubucin, ActD and CX5461) in the absence or presence of BCL-2.
Final libraries, prepared in triplicate, were sequenced on an Illumina Hi-Seq 2000with paired-end 2 x 51-bp read lengths.
In total, the RNA-seq cohort inlcudes 24 samples.
The below bioinformatic sequencing pipepline has been applied for each sample.

**01. download_genome.pbs**                    : download genome annotation (GRCh38)

**02. quality_check.pbs**                      : quality check of sequencing reads using FastQC

**03. trimming.pbs**                           : trimming of sequencing reads using BBDuk suite

**04. quality_check_post_trimming.pbs**        : quality check after the trimming step

**05. adapter_removal.pbs**                    : removal of adapter sequences

**06. quality_check_post_adapter.pbs**         : quality check after the adapter removal step 

**07. generate_genome_indexes.pbs**            : building a STAR index file/create a genome index


**infer_library_strandness.pbs**: verify the type of sequencing library considerinf infer_experiment.py (avialable at https://rseqc.sourceforge.net/#infer-experiment-py)
