#! bin/bash
# ccQTL (cell-type composition QTL) run: params.mode = 'ccqtl' is set in ccqtl_healthy.config


# Remove old logs but not the most previous run
# rm -r *html.*;
# rm .nextflow.log.*;
# rm flowchart.png.*;
# rm trace.txt.*;

# TODO: User edit this section to fit workstation
# export NF_PATH="/software/hgi/installs/anaconda3/envs/nextflow/bin/"
export EQTL_REPO="/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/Dart_BTNLs/QTLight"
export CONFIG="${EQTL_REPO}/conf/extra_confs/sanger-cd_healthy.config"
export INPUT="${EQTL_REPO}/sample_input/input_cd_healthy.nf"
export CCQTL_CONFIG="${EQTL_REPO}/conf/extra_confs/ccqtl_healthy.config"

# Add Singularity to path
# PATH=$PATH:/software/singularity/v3.10.0/bin
module load cellgen/singularity
module load nextflow-23.10.0

# Nextflow settings
export NXF_OPTS="-Xms25G -Xmx25G"
# Uncomment this if get strange bus errors
# export NXF_OPTS="${NXF_OPTS} -Dleveldb.mmap=false" # No resume functionality
export NXF_HOME=${EQTL_REPO}
export NXF_WORK="${NXF_HOME}/work"
export NXF_TEMP="${NXF_HOME}/nextflow_temp"
export NXF_CONDA_CACHEDIR="${NXF_HOME}/nextflow_conda"
export NXF_SINGULARITY_CACHEDIR="${NXF_HOME}/cache_singularity"

# Farm specific settings
#export QT_QPA_PLATFORM='offscreen'
export LSB_DEFAULT_JOBGROUP="/${USER}/nf"
# export LSB_DEFAULTGROUP="sc-eqtl-ibd"
# export LSB_DEFAULTGROUP="cnv_15x"
export LSB_DEFAULTGROUP="team152"
export PYTHONHASHSEED=0

mkdir -p $NXF_TEMP
# mkdir results
# cp ${INPUT} results/
# cp ${EQTL_REPO}/conf/extra_confs/* results/

# To stop TclError
#export QT_QPA_PLATFORM='offscreen'

# Run nextflow
nextflow run ${EQTL_REPO} \
-profile sanger \
-resume \
-c ${CONFIG} \
-c ${INPUT} \
-c ${CCQTL_CONFIG}

# Run in a week queue
# tmux new -s eqtl
# bsub -Is -R "select[mem=30000] rusage[mem=30000]" -M 30000 -q oversubscribed bash
# bash run_ccQTL.sh

