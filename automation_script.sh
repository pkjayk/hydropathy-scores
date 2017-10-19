# loop thru each file in directory
for file in ./human_subsample_7/*
do
  # run command for file
  python assignment3.py -i "$file" -o protein_sequence_output.csv -w 9 
done
