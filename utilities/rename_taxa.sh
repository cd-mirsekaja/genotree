for file in ~/master_input/all_hits_aligned_original/*-aligned.fasta; do
    echo editing $file
    python3 rename_script.py -f $file -d ~/master_input/all_hits_aligned_renamed
done

