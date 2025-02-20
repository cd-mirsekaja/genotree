
threshold=0.35

alignment_files=$(find ~/master_input/all_hits_aligned_renamed/*-renamed.fasta -type f -printf "%f\n" | grep -v -x -f ~/master_input/forbidden_loci.txt | head)

for file in alignment_files; do
    # get ID for current locus
    locus_id=$(echo "${file##*/}" | cut -d'-' -f1)
    aligroove_matrix=~/master_input/AliGROOVE_Matrices/AliGROOVE_seqsim_matrix_$locus_id-renamed.txt
    echo Locus ID: $locus_id
    # filters alignments for current locus for all exons with a score under the given threshold
    python3 ~/genotree/2-2_filter_alignments.py -d $aligroove_matrix -a $file -l $locus_id -t $threshold
done