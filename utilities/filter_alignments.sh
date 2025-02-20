# load necessary modules
module load hpc-env/13.1
module load parallel/20230822-GCCcore-13.1.0
module load Python/3.11.3-GCCcore-13.1.0


threshold=0.2

aligroove_files=$(find ~/master_input/AliGROOVE_Matrices/*matrix.txt -type f -printf "%f\n" | grep -v -x -f ~/master_input/forbidden_loci.txt)
aligroove_count=$(echo $aligroove_files | wc -w)


filter_function(){
    aligroove_matrix=${1}
    # get ID for current locus
    locus_id=$(echo "${aligroove_matrix##*/}" | cut -d'-' -f1)

    echo Locus ID: $locus_id
    alignment_file=~/master_input/all_hits_aligned_renamed/$locus_id-renamed.fasta

    # filters alignments for current locus for all exons with a score under the given threshold
    python3 ~/genotree/2-2_filter_alignments.py -d ~/master_input/AliGROOVE_Matrices/$aligroove_matrix -a $alignment_file -l $locus_id -t $threshold
}

export -f filter_function

for aligroove_matrix in $aligroove_files; do
    filter_function $aligroove_matrix
done