for file in ./human_subsample_7/*
do
  python assignment3.py -i "$file" -o protein_sequence_output.csv -w 9 
done
