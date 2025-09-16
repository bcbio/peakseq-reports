# info params


coldata_fn='/path/to/nf-core/samplesheet.csv'
# This folder is in the nf-core output directory inside multiqc folder
multiqc_data_dir='/path/to/nf-core/output/multiqc/narrow_peak/multiqc_data/'
# This folder is in the nf-core output directory, maybe is broadPeak instead of narrowPeak
peaks_dir = '/path/to/nf-core/output/bowtie2/merged_library/macs3/narrow_peak/'
# This folder is in the nf-core output directory, maybe is broadPeak instead of narrowPeak, also includes antibody name
counts_fn = '/path/to/nf-core/output/bowtie2/merged_library/macs3/narrow_peak/consensus/deseq2/consensus_peaks.mLb.clN.rds'
# This folder is in the nf-core output directory,
inserts_dir <- '/path/to/nf-core/output/bowtie2/merged_library/picard_metrics/'
