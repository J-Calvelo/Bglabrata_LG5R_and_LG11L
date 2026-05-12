#!/usr/bin/env bash

# This script takes the coordinates for the final work regions and extracts the sequences from the region with the intent to conduct dotplot analysis.

region_input_file=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Selected_Regions_input/Selected_Regions_Additional_peaks_from_iM_Zhang_2024.in
storage_sequences=Region_Full_Sequences

if [ -d $storage_sequences ]
then
  rm -rf $storage_sequences
fi
mkdir $storage_sequences

count=2
walker=$(grep -c . $region_input_file)

while [ $count -le $walker ]
do
  regionID=$(sed -n $count"p" $region_input_file | awk -F "\t" '{print $1}')
  AssemblyID=$(sed -n $count"p" $region_input_file  | awk -F "\t" '{print $2}')
  Sequence_File=$(sed -n $count"p" $region_input_file  | awk -F "\t" '{print $3}')
  raw_location=$(sed -n $count"p" $region_input_file | awk -F "\t" '{print $4}')

  echo "Working on: "$regionID": "$AssemblyID

  chr=$(echo $raw_location | awk -F ":" '{print $1}' )
  coords=$(echo $raw_location | tr "-" ":" | tr "+" ":" | awk -F ":" '{print $2":"$3}' )
  strand=$(echo $raw_location | grep -o . | tail -n 1)

  rename_seq=$(echo ">"$regionID"_"$AssemblyID)
  echo $regionID" "$AssemblyID" "$chr" "$coords" "$Sequence_File | tr " " "\t" >> $storage_sequences"/Work_resgistry.fasta"

  if [ $strand == "+" ]
  then
    seqkit grep -p $chr $Sequence_File |  seqkit seq -u -w 60 | seqkit subseq -r $coords | sed "s/>.*/$rename_seq/" >>  $storage_sequences"/"$AssemblyID"_"$regionID".fasta"
  else
    seqkit grep -p $chr $Sequence_File | seqkit seq -u -r -p -t DNA -v -w 60 | seqkit subseq -r $coords | sed "s/>.*/$rename_seq/" >> $storage_sequences"/"$AssemblyID"_"$regionID".fasta"
  fi
  count=$(($count + 1))
done
