#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
PROGNAME=$(basename $SRCNAME)
VERSION="1.4.0"
NUM_POS_ARGS="3"

source "$(dirname ${SRCNAME})/snpnet_misc.sh"

############################################################
# functions
############################################################

show_default_helper () {
    cat ${SRCNAME} | grep -n Default | tail -n+3 | awk -v FS=':' '{print $1}' | tr "\n" "\t"
}

show_default () {
    cat ${SRCNAME} \
        | tail -n+$(show_default_helper | awk -v FS='\t' '{print $1+1}') \
        | head  -n$(show_default_helper | awk -v FS='\t' '{print $2-$1-1}')
}

usage () {
cat <<- EOF
	$PROGNAME (version $VERSION)
	Run snpnet and compute PRS for all individuals in the plink2 pgen file.
	(we will add evaluation, plot, etc. in the next update of the pipeline)

	Usage: $PROGNAME [options] genotype_pfile phenotype_name results_dir
	  genotype_pfile  The plink2 pgen/pvar.zst/psam file.
	  phenotype_name  The name of the phenotype. We assume the phenotype is stored with the same column name
	  results_dir     The results directory. The script will write the following files:
	                   - snpnet.tsv         The BETAs for genotypes
	                   - snpnet.covars.tsv  The BETAs for covariates (when specified)

	Options:
	  --snpnet_dir       Specify the directory of the snpnet package
	  --nCores     (-t)  Number of CPU cores
	  --memory     (-m)  The memory amount

	Default configurations for snpnet (please use the options above to modify them):
	  snpnet_dir=${snpnet_dir}
EOF
    show_default | awk -v spacer="  " '{print spacer $0}'
}

############################################################
# tmp dir
############################################################
tmp_dir_root="$LOCAL_SCRATCH"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

############################################################
# parser start
############################################################
## == Default parameters (start) == ##
nCores=4
mem=30000
## == Default parameters (end) == ##

declare -a params=()
for OPT in "$@" ; do
    case "$OPT" in
        '-h' | '--help' )
            usage >&2 ; exit 0 ;
            ;;
        '-v' | '--version' )
            echo $VERSION ; exit 0 ;
            ;;
        '--snpnet_dir' )
            snpnet_dir=$2 ; shift 2 ;
            ;;
        '-t' | '--cpus' | '--nCores' )
            nCores=$2 ; shift 2 ;
            ;;
        '-m' | '--mem' | '--memory' )
            mem=$2 ; shift 2 ;
            ;;
        '--'|'-' )
            shift 1 ; params+=( "$@" ) ; break
            ;;
        -*)
            echo "$PROGNAME: illegal option -- '$(echo $1 | sed 's/^-*//')'" 1>&2 ; exit 1
            ;;
        *)
            if [[ $# -gt 0 ]] && [[ ! "$1" =~ ^-+ ]]; then
                params+=( "$1" ) ; shift 1
            fi
            ;;
    esac
done

if [ ${#params[@]} -lt ${NUM_POS_ARGS} ]; then
    echo "${PROGNAME}: ${NUM_POS_ARGS} positional arguments are required" >&2
    usage >&2 ; exit 1 ;
fi

genotype_pfile="${params[0]}"
phenotype_name="${params[1]}"
results_dir=$(readlink -f ${params[2]})

############################################################
if [ ! -d ${results_dir} ] ; then
    echo "The specified results dir (${results_dir}) does not exist!" >&2
    exit 1
fi

prevIter="$(find_prevIter ${results_dir})"
rdata="${results_dir}/results/output_iter_${prevIter}.RData"

Rscript "$(dirname ${SRCNAME})/export_betas.R" ${rdata}

ln -sf ${rdata%.RData}.tsv        ${results_dir}/snpnet.iter_${prevIter}.tsv
ln -sf ${rdata%.RData}.covars.tsv ${results_dir}/snpnet.iter_${prevIter}.covars.tsv

ln -sf ${results_dir}/snpnet.iter_${prevIter}.tsv        ${results_dir}/snpnet.tsv
ln -sf ${results_dir}/snpnet.iter_${prevIter}.covars.tsv ${results_dir}/snpnet.covars.tsv

plink_score ${results_dir} ${phenotype_name} ${genotype_pfile} ${nCores} ${mem}

rm ${results_dir}/snpnet.tsv
rm ${results_dir}/snpnet.covars.tsv
