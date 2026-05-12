#!/usr/bin/env bash

# Input information
main_genome=iM_Zhang_2024
other_anotation_transfers=$(echo "BS90_bu_et_al_qtl Bg_2024 iM_bu_et_al_qtl")
annotation_transfers=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Annotation_Transfer_Liftoff/Individual_Transfers
genome_input_table=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Input_assemblies.in
targetp_ids_guide=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Ref_Genes_Annotation/TargetP
signap6_folder=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Ref_Genes_Annotation/SignalP6
tmhmm_data_folder=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Ref_Genes_Annotation/TMHMM

window_lenght=2000000
window_movement=500000
liftoff_coverage_cut_off=0.80

# Run_information
work_dir=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation
original_annotation_gff=$(awk -F "\t" -v genome_ID=$main_genome '{if ($1==genome_ID) print $4}' $genome_input_table)
original_genome_fasta=$(awk -F "\t" -v genome_ID=$main_genome '{if ($1==genome_ID) print $3}' $genome_input_table)
tmmm_information_file=$(echo $targetp_ids_guide"/"$main_genome"_tmhmm.txt")

# Modules
RUN_GET_DATA=FALSE
RUN_SLIDING_WINDOW=TRUE

# Functions
get_tmhmm () {
  count_gene=1
  walker=$(grep -c . $current_gff)

  echo "GeneID ProtID Prot_num Chr Start End Strand Liftoff_Coverage N_TMHMMs SP_Warning N_Work_TMHMMs SignalP6_Hit CS_SP TargetP_Hit CS_mTP" | tr " " "\t" >> $work_dir"/TMHMM_Density_Ref_"$main_genome"/Genes_TMHMM_from_"$current_anot".txt"
  while [ $count_gene -le $walker ]
  do
    gene_id=$(sed -n $count_gene"p" $current_gff | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | sed "s/ID=gene-//")
    awk -F "\t" -v gen=$gene_id '{if ($1==gen) print $2}' $targetp_ids_guide"/"$current_anot"_TargetP.txt" > $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/temp_prot_ids.tmp"
    check_prot_id=$(grep -c . $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/temp_prot_ids.tmp")
    if [ $check_prot_id -gt 0 ]
    then
      location=$(sed -n $count_gene"p" $current_gff | awk -F "\t" '{print $1" "$4" "$5" "$7}')
      if [ $original == "TRUE" ]
      then
        coverage=1
      else
        coverage=$(sed -n $count_gene"p" $current_gff | awk -F "\t" '{print $9}' | tr ";" "\n" | grep "coverage=" | awk -F "=" '{print $2}')
      fi

      protein_ids=$(cat $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/temp_prot_ids.tmp")
      num_prot=1
      for prot in $protein_ids
      do
        echo $current_anot": "$prot
        prot_num_print=$(echo $num_prot"_of_"$check_prot_id)
        grep $prot $tmhmm_data_folder"/"$current_anot"_tmhmm.txt" > $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/tmhmm_data.tmp"

        targetP_hit=$(awk -F "\t" -v protID=$prot '{if ($2==protID) print $3}' $targetp_ids_guide"/"$current_anot"_TargetP.txt" | awk '{if ($1=="mTP") print "X"; else print "."}')
        if [ $targetP_hit == "X" ]
        then
          targetP_pos=$(awk -F "\t" -v protID=$prot '{if ($2==protID) print $7}' $targetp_ids_guide"/"$current_anot"_TargetP.txt"| awk -F "." '{print $1}' | sed "s/CS pos: //")
        else
          targetP_pos=NA
        fi

        signalp_hit=$(awk -F "\t" -v protID=$prot '{if ($1==protID) print $2}' $signap6_folder"/"$current_anot"/prediction_results.txt" | awk '{if ($1=="SP") print "X"; else print "."}')
        if [ $signalp_hit == "X" ]
        then
          signalp_pos=$(awk -F "\t" -v protID=$prot '{if ($1==protID) print $5}' $signap6_folder"/"$current_anot"/prediction_results.txt" | awk -F "." '{print $1}' | sed "s/CS pos: //")
        else
          signalp_pos=NA
        fi

        raw_tmms=$(grep "Number of predicted TMHs:" $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/tmhmm_data.tmp" | awk -F ":" '{print $2}' | tr -d " ")
        check_tmm_too_close=FALSE
        if [ $raw_tmms -eq 0 ]
        then
          work_tmms=0
        else
          if [ $targetP_hit == "X" ]
          then
            boundary=$(echo $targetP_pos | awk -F "-" '{print $2}')
            check_tmm_too_close=$(grep -v "#" $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/tmhmm_data.tmp" | grep -w TMhelix | head -n 1 | awk -v boundary=$boundary '{if ($4<=boundary) print "TRUE"; else print "FALSE"}')
          elif [ $signalp_hit == "X" ]
          then
            boundary=$(echo $signalp_pos | awk -F "-" '{print $2}')
            check_tmm_too_close=$(grep -v "#" $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/tmhmm_data.tmp" | grep -w TMhelix | head -n 1 | awk -v boundary=$boundary '{if ($4<=boundary) print "TRUE"; else print "FALSE"}')
          fi

          if [ $check_tmm_too_close == TRUE ]
          then
            work_tmms=$(($raw_tmms - 1))
          else
            work_tmms=$raw_tmms
          fi
        fi

        num_prot=$(($num_prot + 1))
        echo $gene_id" "$prot" "$prot_num_print" "$location" "$coverage" "$raw_tmms" "$check_tmm_too_close" "$work_tmms" "$signalp_hit" "$signalp_pos" "$targetP_hit" "$targetP_pos | tr " " "\t" >> $work_dir"/TMHMM_Density_Ref_"$main_genome"/Genes_TMHMM_from_"$current_anot".txt"
      done
    fi
    count_gene=$(($count_gene + 1))
  done
}

if [ $RUN_GET_DATA == TRUE ]
then
  # Original Annotation information
  if [ -d $work_dir"/TMHMM_Density_Ref_"$main_genome ]
  then
    rm -rf $work_dir"/TMHMM_Density_Ref_"$main_genome
  fi

  mkdir $work_dir"/TMHMM_Density_Ref_"$main_genome
  mkdir $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp"
  seqkit fx2tab -n -i -l $original_genome_fasta > $work_dir"/TMHMM_Density_Ref_"$main_genome"/genome.len"
  awk -F "\t" '{if ($3=="gene") print}' $original_annotation_gff > $work_dir"/TMHMM_Density_Ref_"$main_genome"/original_genes.gff"

  original=TRUE
  current_anot=$main_genome
  awk -F "\t" '{if ($3=="gene") print}' $original_annotation_gff > $work_dir"/TMHMM_Density_Ref_"$main_genome"/original_genes.gff"
  current_gff=$work_dir"/TMHMM_Density_Ref_"$main_genome"/original_genes.gff"
  get_tmhmm

  for alt in $other_anotation_transfers
  do
    original=FALSE
    current_anot=$alt
    awk -F "\t" '{if ($3=="gene") print}' $annotation_transfers"/"$alt"_transfered_to_"$main_genome"_annotation.gff_polished" > $work_dir"/TMHMM_Density_Ref_"$main_genome"/Transfer_genes_from_"$alt"_genes.gff"
    current_gff=$work_dir"/TMHMM_Density_Ref_"$main_genome"/Transfer_genes_from_"$alt"_genes.gff"
    get_tmhmm
  done
fi

################################################################################
################################################################################

if [ $RUN_SLIDING_WINDOW == "TRUE" ]
then
  if [ ! -d $work_dir"/TMHMM_Density_Ref_"$main_genome ]
  then
    echo "ABORT!! No folder: "$work_dir"/TMHMM_Density_Ref_"$main_genome
    exit
  fi

  if [ -d $work_dir"/TMHMM_Density_Ref_"$main_genome"/Window_Measurements_"$window_lenght"_jump_"$window_movement ]
  then
    rm -rf $work_dir"/TMHMM_Density_Ref_"$main_genome"/Window_Measurements_"$window_lenght"_jump_"$window_movement
  fi

  mkdir $work_dir"/TMHMM_Density_Ref_"$main_genome"/Window_Measurements_"$window_lenght"_jump_"$window_movement

  work_anots=$(echo $main_genome" "$other_anotation_transfers)
  chr_ids=$(awk -F "\t" '{print $1}' $work_dir"/TMHMM_Density_Ref_"$main_genome"/genome.len" | sort -u)

  echo "WindowID Chr Start End First_Gene Last_Gene N_Genes N_genes_TMHMM N_Genes_Single_TMHMM Per_genes_TMHMM Per_Genes_Single_TMHMM Warning_Truncated_Window Annot_source" | tr " " "\t" >> $work_dir"/TMHMM_Density_Ref_"$main_genome"/Window_Measurements_"$window_lenght"_jump_"$window_movement"/TMHMM_on_window_anot.txt"
  echo "GeneID Chr Start End N_TMHMMs WindowID Annot_source" | tr " " "\t" > $work_dir"/TMHMM_Density_Ref_"$main_genome"/Window_Measurements_"$window_lenght"_jump_"$window_movement"/Gene_Registry.txt"

  N_current_window=1
  for current_anot in $work_anots
  do
    for chrID in $chr_ids
    do
      chr_len=$(awk -F "\t" -v chr=$chrID '{if ($1==chr) print $2}' $work_dir"/TMHMM_Density_Ref_"$main_genome"/genome.len")
      start_coord=1
      end_coord=$window_lenght

      awk -F "\t" -v chr=$chrID -v lifftof_cov=$liftoff_coverage_cut_off '{if ($4==chr && $8 >= lifftof_cov) print}' $work_dir"/TMHMM_Density_Ref_"$main_genome"/Genes_TMHMM_from_"$current_anot".txt" > $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_TMHMM_raw_data.txt"

      if [ -f  $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_TMHMM_filter_data.txt" ]
      then
        rm -f $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_TMHMM_filter_data.txt"
      fi

      traverse_genes=$(awk -F "\t" '{print $1}'  $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_TMHMM_raw_data.txt" | sort -u)
      for gene in $traverse_genes
      do
        echo "Checking gene multi-proteins: "$gene
        awk -F "\t" -v gen=$gene '{if ($1==gen) print }' $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_TMHMM_raw_data.txt"  | sort -n -k 11 | tail -n 1 >> $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_TMHMM_filter_data.txt"
      done

      if [ -f $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_TMHMM_filter_data.txt" ]
      then
        while [ $start_coord -le $chr_len ]
        do
          warning=$(echo $end_coord" "$chr_len | awk -F " " '{if ($2<$1) print "X"; else print "."}')

          echo "Working on window: "$N_current_window"("$start_coord" to "$end_coord")"
          awk -F "\t" -v start=$start_coord -v end=$end_coord '{if ($6>=start && $5<=end) print}' $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_TMHMM_filter_data.txt" > $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_window_data.txt"

          N_genes=$(awk -F "\t" '{print $1}' $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_window_data.txt" | grep -c .)
          if [ $N_genes -gt 0 ]
          then
            First_gene=$(head -n 1 $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_window_data.txt" | awk -F "\t" '{print $1}')
            Last_gene=$(tail -n 1 $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_window_data.txt" | awk -F "\t" '{print $1}')
            N_genes_TMHMM=$(awk -F "\t" '{if ($11 > 0) print $1}' $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_window_data.txt" | grep -c .)
            N_genes_TMHMM_1=$(awk -F "\t" '{if ($11 == 1) print $1}' $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_window_data.txt" | grep -c .)
            Per_genes_TMHMM=$(echo $N_genes" "$N_genes_TMHMM | awk -F " " '{printf "%.3f\n", $2/$1}' | awk '{print $1*100"%"}')
            Per_genes_TMHMM_1=$(echo $N_genes" "$N_genes_TMHMM_1 | awk -F " " '{printf "%.3f\n", $2/$1}' | awk '{print $1*100"%"}')
          else
            First_gene=NA
            Last_gene=NA
            N_genes_TMHMM=0
            N_genes_TMHMM_1=0
            Per_genes_TMHMM=0%
            Per_genes_TMHMM_1=0%
          fi

          awk -F "\t" -v windowID=$N_current_window -v anot=$current_anot '{print $1"\t"$4"\t"$5"\t"$6"\t"$11"\t"windowID"\t"anot}' $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_window_data.txt"  >> $work_dir"/TMHMM_Density_Ref_"$main_genome"/Window_Measurements_"$window_lenght"_jump_"$window_movement"/Gene_Registry.txt"
          echo $N_current_window" "$chrID" "$start_coord" "$end_coord" "$First_gene" "$Last_gene" "$N_genes" "$N_genes_TMHMM" "$N_genes_TMHMM_1" "$Per_genes_TMHMM" "$Per_genes_TMHMM_1" "$warning" "$current_anot | tr " " "\t" >> $work_dir"/TMHMM_Density_Ref_"$main_genome"/Window_Measurements_"$window_lenght"_jump_"$window_movement"/TMHMM_on_window_anot.txt"

          start_coord=$(($start_coord + $window_movement))
          end_coord=$(($start_coord + $window_lenght - 1))

          if [ $warning == "X" ]
          then
            start_coord=$(($chr_len + $chr_len))
          fi

          N_current_window=$(($N_current_window + 1))
        done
        rm -f $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_window_data.txt"
        rm -f $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_TMHMM_filter_data.txt"
      fi
      rm -f $work_dir"/TMHMM_Density_Ref_"$main_genome"/Temp/Current_TMHMM_raw_data.txt"
    done
  done
fi
