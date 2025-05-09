qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2 --input-path sample_list --output-path seqs/seqs.qza
qiime demux summarize --i-data seqs/seqs.qza --o-visualization seqs/seqs.qzv

qiime dada2 denoise-paired --i-demultiplexed-seqs seqs/seqs.qza --p-trim-left-f 26 --p-trim-left-r 26 --p-trunc-len-f 296 --p-trunc-len-r 295 \
    --o-representative-sequences dada2/rep_seqs.qza --o-table dada2/table.qza --o-denoising-stats dada2/stats.qza --p-n-threads 100
qiime tools export --input-path dada2/table.qza --output-path dada2 && biom convert -i dada2/feature-table.biom -o dada2/table.tsv --to-tsv

qiime feature-table filter-features --i-table dada2/table.qza --p-min-frequency 1 --p-min-samples 2 --o-filtered-table dada2/table_f.qza
qiime tools export --input-path dada2/table_f.qza --output-path dada2/ && biom convert -i dada2/feature-table.biom -o dada2/table_f.tsv --to-tsv
qiime feature-table summarize --i-table dada2/table_f.qza --o-visualization dada2/table_f.qzv --m-sample-metadata-file sample_group
qiime feature-table filter-seqs --i-data dada2/rep_seqs.qza --i-table dada2/table_f.qza --o-filtered-data dada2/rep_seqs_f.qza
qiime tools export --input-path dada2/rep_seqs_f.qza --output-path dada2/

qiime phylogeny align-to-tree-mafft-fasttree --i-sequences dada2/rep_seqs_f.qza --o-alignment tree/aligned.qza \
    --o-masked-alignment tree/masked-aligned.qza --o-tree tree/unrooted-tree.qza --o-rooted-tree tree/rooted-tree.qza --p-n-threads 56
qiime tools export --input-path tree/rooted-tree.qza  --output-path tree/ && mv tree/tree.nwk tree/rooted_tree.nwk

qiime diversity core-metrics --i-table dada2/table_f.qza --p-sampling-depth 27279 --m-metadata-file sample_group --p-n-jobs 56 \
    --o-rarefied-table div/table_rarefied.qza --output-dir div/core_div
qiime tools export --input-path div/table_rarefied.qza --output-path div/ && biom convert -i div/feature-table.biom -o div/table_rarefied.tsv --to-tsv

qiime feature-classifier classify-sklearn --i-classifier /data/database/silva/338F_806R/silva-138-99-338F_806R_classifier.qza --i-reads dada2/rep_seqs_f.qza \
    --o-classification tax/taxonomy.qza


