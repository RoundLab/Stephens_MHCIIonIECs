!/bin/bash

# Requires:
# - Singularity 3.1.1
# - Filled in paths, environmental variables and parameters below.

# Set QIIME2 container image location to use with Singularity
ContainerImage=docker://qiime2/core:2019.4

#### Directory paths #####
# Will be made if not already existing. All 3 must be listed.
SCRATCH=
# Main results directory for core files (table, tree, taxonomy) and analyses. ResultDir is 16S/results/18684R_IlealSepCo_Homeo/ for this run.
ResultDir=
# Working directory for all artifact files including intermediate input and trimmed sequence artifact files.
WRKDIR=
##########################
#### User Input Files ####
# Use full paths
# Manifest format (for this QIIME2 version) should be a CSV with headers: sample-id,absolute-filepath,direction
MANIFEST=
# A metadata map file. Used the file 16S/metadata/map_18684R_IlealSepCo_Homeo.txt for this run.
MAP=
# A classifier artifact. Used Greengenes 13_8 trained on v3-v4 region targeted by listed primers.
# Note that this seq run uses different primers (that target the same regions) than the other 3 16S seq runs in this project.
CLASSIFIER=
##########################
#### QIIME 2 Param #######
Nproc=16
JoinedTrimLength=392
FPrimerSeqToTrim=TGCCTACGGGNBGCASCAG
RPrimerSeqToTrim=GCGACTACNVGGGTATCTAATCC
MinMergeLength=189
##########################
#### Setup ###############
mkdir -p $SCRATCH
mkdir -p ${SCRATCH}/tmp_XDG
mkdir -p ${WRKDIR}
mkdir -p ${ResultDir}/q2_viz
##########################
#### SINGULARITYENV ###########
XDG_RUNTIME_DIR=${SCRATCH}/tmp_XDG
export SINGULARITYENV_UHOME=${HOME}
export SINGULARITYENV_MANIFEST=${MANIFEST}
export SINGULARITYENV_SCRATCH=${SCRATCH}
export SINGULARITYENV_WRKDIR=${WRKDIR}
export SINGULARITYENV_XDG_RUNTIME_DIR=${XDG_RUNTIME_DIR}
export SINGULARITYENV_Nproc=${Nproc}
export SINGULARITYENV_FPrimerSeqToTrim=${FPrimerSeqToTrim}
export SINGULARITYENV_RPrimerSeqToTrim=${RPrimerSeqToTrim}
export SINGULARITYENV_MinMergeLength=${MinMergeLength}
export SINGULARITYENV_MAP=${MAP}
export SINGULARITYENV_CLASSIFIER=${CLASSIFIER}
###############################

## Import and process sequences through table creation, taxonomy calls and phylogeny creation.

cd ${SCRATCH}

# Part 1: Import

echo "TIME: START import = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ${MANIFEST} \
--output-path ${SCRATCH}/inseqs-demux.qza \
--input-format PairedEndFastqManifestPhred33

singularity exec ${ContainerImage} qiime demux summarize \
  --i-data ${SCRATCH}/inseqs-demux.qza \
  --o-visualization ${SCRATCH}/inseq-demux-summary.qzv  
echo "TIME: END import = `date +"%Y-%m-%d %T"`"

# Copy visualization results to ResultDir and artifacts to working dir
cp ${SCRATCH}/*.qzv ${ResultDir}/q2_viz
cp ${SCRATCH}/inseqs-demux.qza ${WRKDIR}/ # Generally, it is not useful and too big to copy full input sequence artifacts

# Part 2: Clean, denoise, filter chimeras, create table and phylogeny, call taxonomy

echo "TIME: START trim, merge = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime cutadapt trim-paired \
--i-demultiplexed-sequences ${SCRATCH}/inseqs-demux.qza \
--o-trimmed-sequences PE-demux_trim.qza \
--p-front-f ${FPrimerSeqToTrim} \
--p-front-r ${RPrimerSeqToTrim} \
--p-cores $Nproc

singularity exec ${ContainerImage} qiime vsearch join-pairs \
--i-demultiplexed-seqs PE-demux_trim.qza \
--o-joined-sequences PE-demux_trim_join.qza \
--p-minmergelen ${MinMergeLength} \
--verbose \

singularity exec ${ContainerImage} qiime demux summarize \
  --i-data PE-demux_trim_join.qza \
  --o-visualization PE-demux_trim_join.qzv
  
#At this stage plots should be inspected to infer the JoinedTrimLength for deblur. 

singularity exec ${ContainerImage} qiime quality-filter q-score-joined \
--i-demux PE-demux_trim_join.qza \
--o-filtered-sequences PE-demux_trim_join_filt.qza \
--o-filter-stats PE-demux_trim_join_filt_stats.qza \
--p-min-quality 10 
echo "TIME: END trim, merge = `date +"%Y-%m-%d %T"`"

echo "TIME: START denoise = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime deblur denoise-16S \
--i-demultiplexed-seqs PE-demux_trim_join_filt.qza \
--p-trim-length $JoinedTrimLength \
--p-jobs-to-start $Nproc \
--o-table table.qza \
--o-representative-sequences repseq.qza \
--o-stats table_stats.qza

singularity exec ${ContainerImage} qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv
echo "TIME: END denoise = `date +"%Y-%m-%d %T"`"

# Chimera filtering using "include borderline chimeras parameters"

echo "TIME: START chimera filter = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime vsearch uchime-denovo \
--i-table table.qza \
--i-sequences repseq.qza \
--output-dir uchime_out
  
singularity exec ${ContainerImage} qiime feature-table filter-features \
--i-table table.qza \
--m-metadata-file uchime_out/chimeras.qza \
--p-exclude-ids \
--o-filtered-table table_nochim.qza
  
singularity exec ${ContainerImage} qiime feature-table filter-seqs \
--i-data repseq.qza \
--m-metadata-file uchime_out/chimeras.qza \
--p-exclude-ids \
--o-filtered-data repseq_nochim.qza
  
singularity exec ${ContainerImage} qiime feature-table summarize \
--i-table table_nochim.qza \
--o-visualization table_nochim.qzv
--m-sample-metadata-file ${MAP}

singularity exec ${ContainerImage} qiime feature-table tabulate-seqs \
--i-data repseq_nochim.qza \
--o-visualization repseq_nochim.qzv
echo "TIME: END chimera filter = `date +"%Y-%m-%d %T"`"

echo "TIME: START phylogeny = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences repseq_nochim.qza \
--o-alignment aligned_repseq_nochim.qza \
--o-masked-alignment masked_aligned_repseq_nochim.qza \
--o-tree tree_unroot.qza \
--o-rooted-tree tree_root.qza
echo "TIME: END phylogeny = `date +"%Y-%m-%d %T"`"

echo "TIME: START taxonomy = `date +"%Y-%m-%d %T"`"
singularity exec ${ContainerImage} qiime feature-classifier classify-sklearn \
--i-classifier ${CLASSIFIER} \
--i-reads repseq_nochim.qza \
--o-classification taxonomy.qza \
--p-n-jobs ${Nproc}

singularity exec ${ContainerImage} qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

singularity exec ${ContainerImage} qiime taxa barplot \
--i-table table_nochim.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file ${MAP} \
--o-visualization table_taxbarplots.qzv
echo "TIME: END taxonomy = `date +"%Y-%m-%d %T"`"

cp *.qzv ${ResultDir}/q2_viz/; cp repseq_nochim.qza ${ResultDir}/; cp table_nochim.qza ${ResultDir}/; cp taxonomy.qza ${ResultDir}/; cp tree_root.qza ${ResultDir}/; 
cp *.qza ${WRKDIR}/
rm -R ${SCRATCH}

## Analysis Section
ContainerImage=docker://qiime2/core:2020.6
cd $ResultDir
mkdir filt_tables
mkdir ancom

# Remove a few ASVs that were not assigned to even a phylum

singularity exec ${ContainerImage} qiime taxa filter-table --i-table table_nochim.qza --i-taxonomy taxonomy.qza --o-filtered-table filt_tables/Ileal.qza --p-include p__

# Break up by housing and assess difference in separately housed and cohoused groups:
singularity exec ${ContainerImage} qiime feature-table filter-samples --i-table filt_tables/Ileal.qza --o-filtered-table filt_tables/Ileal_Sep.qza --m-metadata-file $MAP --p-where "[HousingByGenotype]='Sep'"
singularity exec ${ContainerImage} qiime feature-table filter-samples --i-table filt_tables/Ileal.qza --o-filtered-table filt_tables/Ileal_Co.qza --m-metadata-file $MAP --p-where "[HousingByGenotype]='Co'"
singularity exec ${ContainerImage} qiime composition add-pseudocount --i-table  filt_tables/Ileal_Co.qza --o-composition-table ancom/Ileal_Co_pc.qza
singularity exec ${ContainerImage} qiime composition add-pseudocount --i-table  filt_tables/Ileal_Sep.qza --o-composition-table ancom/Ileal_Sep_pc.qza

singularity exec ${ContainerImage} qiime composition ancom --i-table ancom/Ileal_Co_pc.qza --o-visualization ancom/Ileal_Co_pc.qzv --m-metadata-file $MAP --m-metadata-column Genotype
singularity exec ${ContainerImage} qiime composition ancom --i-table ancom/Ileal_Sep_pc.qza --o-visualization ancom/Ileal_Sep_pc.qzv --m-metadata-file $MAP --m-metadata-column Genotype

# Assess alpha diversity differences
singularity exec ${ContainerImage} qiime diversity core-metrics-phylogenetic --i-table filt_tables/Ileal_Sep.qza --output-dir corediv_3k_Sep --p-sampling-depth 3000 --p-n-jobs-or-threads 10 --i-phylogeny tree_root.qza --m-metadata-file $MAP

for adiv in corediv_3k_Sep/shannon_vector.qza corediv_3k_Sep/observed_features_vector.qza corediv_3k_Sep/evenness_vector.qza; do singularity exec ${ContainerImage} qiime diversity alpha-group-significance --i-alpha-diversity $adiv --o-visualization ${adiv%_vector.qza}_GroupSig.qzv --m-metadata-file $MAP; done

singularity exec ${ContainerImage} qiime diversity core-metrics-phylogenetic --i-table filt_tables/Ileal_Co.qza --output-dir corediv_3k_Co --p-sampling-depth 3000 --p-n-jobs-or-threads 10 --i-phylogeny tree_root.qza --m-metadata-file $MAP

for adiv in corediv_3k_Co/shannon_vector.qza corediv_3k_Co/observed_features_vector.qza corediv_3k_Co/evenness_vector.qza; do singularity exec ${ContainerImage} qiime diversity alpha-group-significance --i-alpha-diversity $adiv --o-visualization ${adiv%_vector.qza}_GroupSig.qzv --m-metadata-file $MAP; done


