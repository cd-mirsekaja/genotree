for file in ~/master_input/all_hits_aligned/*-aligned.fasta; do
    echo editing $file
    #sed -i .edit -e 's/|/-/g' $file
    python3 rename_script.py -f $file
done