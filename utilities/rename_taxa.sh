for file in ~/master_input/all_hits_aligned_2/*-aligned.fasta; do
    echo editing $file
    python3 rename_script.py -f $file
done