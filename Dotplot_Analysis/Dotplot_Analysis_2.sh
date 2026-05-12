#!/usr/bin/env bash

# This script checks fasta files in a folder matching one of the regions
# and runs the steps necesary for blastn2dotplots.

# You need to specify:
# 1) Folder with the sequences to analyze
# 2) Ids of the regions
# 3) Min alignment lenght and Identity
# 4) File with regions to highlight. They need  to include: 1) Assembly, 2) Region ID, 3) Start, 4) End, 5) Color. With no header
# The script will parse it when adecaute as the blastn2dotplots format.
# 5) Folder where to store the inputs for submission to a cluster. The script will try to delete this folder before making a new one. So please please please take it into account.

work_folder=/nfs5/IB/Blouin_Lab/users/javierc/scripts/Make_Tables/Blast2Dotplot_runs
work_regions=$(echo "LG10_Left LG5_Right" )
aln_len=1000
min_ident=90
highlights_file=Highlights.in
cluster_inputs_dir=/nfs5/IB/Blouin_Lab/users/javierc/scripts/Make_Tables/Dotplot_Cluster

check_high () {
  check_add_highlights=$(awk -F "\t" -v seq=$current_seq -v region=$region '{if ($1==seq && $2==region) print}' $work_folder"/"$highlights_file  | grep -c .)
  if [ $check_add_highlights -gt 0 ]
  then
    run_highlithgs=TRUE
    seq_name=$(seqkit seq -n $work_folder"/"$current_seq"_"$region".fasta")
    awk -F "\t" -v seq=$current_seq -v region=$region -v name=$seq_name '{if ($1==seq && $2==region) print name" "$3" "$4" "$5}' $work_folder"/"$highlights_file | tr " " "\t" >> $storage"/Higlight_info.in"
  fi
}

remove_old_highlights () {
  if [ ! -d $storage ]
  then
    mkdir $storage
  fi

  if [ -f $storage"/Higlight_info.in" ]
  then
    rm -f $storage"/Higlight_info.in"
  fi
}

###############################################################################

if [ -d $cluster_inputs_dir ]
then
  rm -rf $cluster_inputs_dir
fi
mkdir $cluster_inputs_dir

for region in $work_regions
do
  echo $region
  if [ ! -d $work_folder"/"$region ]
  then
    mkdir $work_folder"/"$region
  fi


  echo "Already Done:" > $work_folder"/"$region"/Already_done.temp"

  remove=$(echo "_"$region".fasta")
  assemblies=$(ls $work_folder"/"*"_"$region".fasta" | sed "s/.*\///" | sed "s/$remove//" )

  for seq1 in $assemblies
  do
    for seq2 in $assemblies
    do
      run_highlithgs=FALSE
      check=$(grep -w -F -c $seq1"_x_"$seq2 $work_folder"/"$region"/Already_done.temp")
      if [ $check -eq 0 ]
      then
        echo "Working on: "$seq1" vs "$seq2
        if [ $seq1 == $seq2 ]
        then
          storage=$work_folder"/"$region"/"$seq1"_Self_check"
          remove_old_highlights

          current_seq=$seq1
          check_high

          seqkit seq -n $work_folder"/"$seq1"_"$region".fasta" > $storage"/Input_i1.in"
          echo "makeblastdb -in "$work_folder"/"$seq1"_"$region".fasta -dbtype nucl -parse_seqids -input_type fasta -out "$storage"/Ref_"$seq1 >> "Dotplot_Cluster/1_Blast_DB.sh"
          echo "blastn -query "$work_folder"/"$seq1"_"$region".fasta -db "$storage"/Ref_"$seq1" -outfmt '6 std qlen slen' -out "$storage"/Blast_Results.in" >> "Dotplot_Cluster/2_Blastn.sh"

          if [ $run_highlithgs == TRUE ]
          then
            echo "blastn2dotplots -i1 "$storage"/Input_i1.in --blastn "$storage"/Blast_Results.in --min_alignlen "$aln_len" --min_identity "$min_ident" --line_width 2 --out "$storage"/Dotplot_"$seq1"_vs_self_Aln-"$aln_len"_Min_Ident-"$min_ident" --highlight "$storage"/Higlight_info.in" >> "Dotplot_Cluster/3_Dotplot.sh"
          else
            echo "blastn2dotplots -i1 "$storage"/Input_i1.in --blastn "$storage"/Blast_Results.in --min_alignlen "$aln_len" --min_identity "$min_ident" --line_width 2 --out "$storage"/Dotplot_"$seq1"_vs_self_Aln-"$aln_len"_Min_Ident-"$min_ident >> "Dotplot_Cluster/3_Dotplot.sh"
          fi

          echo $seq1"_x_"$seq2 > $work_folder"/"$region"/Already_done.temp"
        else
          storage=$work_folder"/"$region"/"$seq1"_vs_"$seq2
          remove_old_highlights

          current_seq=$seq1
          check_high
          current_seq=$seq2
          check_high

          seqkit seq -n $work_folder"/"$seq1"_"$region".fasta" > $storage"/Input_i1.in"
          seqkit seq -n $work_folder"/"$seq2"_"$region".fasta" > $storage"/Input_i2.in"

          echo "makeblastdb -in "$work_folder"/"$seq1"_"$region".fasta -dbtype nucl -parse_seqids -input_type fasta -out "$storage"/Ref_"$seq1 >> $cluster_inputs_dir"/1_Blast_DB.sh"
          echo "blastn -query "$work_folder"/"$seq2"_"$region".fasta -db "$storage"/Ref_"$seq1" -outfmt '6 std qlen slen' -out "$storage"/Blast_Results.in" >> $cluster_inputs_dir"/2_Blastn.sh"

          if [ $run_highlithgs == TRUE ]
          then
            echo "blastn2dotplots -i1 "$storage"/Input_i1.in -i2 "$storage"/Input_i2.in --blastn "$storage"/Blast_Results.in --min_alignlen "$aln_len" --min_identity "$min_ident" --line_width 2 --out "$storage"/Dotplot_"$seq1"_vs_"$seq2"_Aln-"$aln_len"_Min_Ident-"$min_ident" --highlight "$storage"/Higlight_info.in" >> $cluster_inputs_dir"/3_Dotplot.sh"
          else
            echo "blastn2dotplots -i1 "$storage"/Input_i1.in -i2 "$storage"/Input_i2.in --blastn "$storage"/Blast_Results.in --min_alignlen "$aln_len" --min_identity "$min_ident" --line_width 2 --out "$storage"/Dotplot_"$seq1"_vs_"$seq2"_Aln-"$aln_len"_Min_Ident-"$min_ident >> $cluster_inputs_dir"/3_Dotplot.sh"
          fi
        fi
      fi
    done
  done
  echo "----------------------------------------------------------------"
done
