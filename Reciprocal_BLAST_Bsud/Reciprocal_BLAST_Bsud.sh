#!/usr/bin/env bash

work_directory=/nfs5/IB/Blouin_Lab/users/javierc/Bsudanica/Reciprocal_Blast
Bglabrata_prot=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Backpup_Annotation_Runs/Bg_2024_GCF_947242115.1_protein.faa
Bglabrata_prot_ID_guide=Bg_Prot_gruide.txt
Bglabrata_gff=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Annotation_Transfer_Liftoff/Individual_Transfers/Bg_2024_transfered_to_iM_Zhang_2024_annotation.gff_polished

Bsudanica_prot=Bsudanica_prot.faa
Bsudanica_gff=Bsudanica_polyA_annotation.gff
Bsudanica_prot_ID_guide=Bsud_Prot_gruide.txt

blast_min_identity=40
blast_min_coverage=70

run_reciprocal_check () {
  awk -F "\t" -v ident=$blast_min_identity -v qcov=$blast_min_coverage '{if ($3>=ident && $13>=qcov) print}' $current_blast > $work_directory"/Best_reciprocal_hits/"$out_reciprocal"_filtered.blastp"

  awk  -F "\t" '{print $1"\t"$2"\t"$12*100}' $work_directory"/Best_reciprocal_hits/"$out_reciprocal"_filtered.blastp" >>  $work_directory"/Best_reciprocal_hits/"$out_reciprocal"_filtered.temp"

  check_genes=$(awk -F "\t" '{print $1}' $work_directory"/Best_reciprocal_hits/"$out_reciprocal"_filtered.temp" | sort -u)

  for gen in $check_genes
  do
    echo "Current: "$out_reciprocal" || "$gen
    best_bitscore=$(awk -F "\t" -v gen=$gen '{if ($1==gen) print $3} '$work_directory"/Best_reciprocal_hits/"$out_reciprocal"_filtered.temp" | sort -n | tail -n 1)
    best_hits=$(awk -F "\t" -v gen=$gen -v bitscore=$best_bitscore '{if ($1==gen && $3==bitscore) print $2}' $work_directory"/Best_reciprocal_hits/"$out_reciprocal"_filtered.temp")

    echo $gen","$best_hits | tr " " ";" | tr "," "\t" >>  $work_directory"/Best_reciprocal_hits/"$out_reciprocal"_hit_to_hit.txt"

    for hit in $best_hits
    do
      echo "XX"$gen"XX XX"$hit"XX" | tr " " "\n" | sort | tr "\n" " " | sed "s/ $/\n/" | sed "s/ /_xxx_/" >> $work_directory"/Best_reciprocal_hits/"$out_reciprocal"_hit_to_hit.work"
    done
  done

  rm -f $work_directory"/Best_reciprocal_hits/"$out_reciprocal"_filtered.temp"
}

RUN_BLASTS=FALSE
RUN_CROSS_BLAST=FALSE
RUN_FINAL_Registry=FALSE
RUN_INTERESTING_LINES=TRUE

if [ $RUN_BLASTS == "TRUE" ]
then
  if [ -d $work_directory"/Ref_Blasts" ]
  then
    rm -rf $work_directory"/Ref_Blasts"
    rm -rf $work_directory"/Blast_Results"
  fi

  if [ -f Run_BLast_cluster.sh ]
  then
    rm -f Run_BLast_cluster.sh
  fi

  mkdir $work_directory"/Ref_Blasts"
  mkdir $work_directory"/Blast_Results"

  makeblastdb -in $Bglabrata_prot -dbtype prot -parse_seqids -input_type fasta -out $work_directory"/Ref_Blasts/Bglabrata_ref"
  makeblastdb -in $Bsudanica_prot -dbtype prot -parse_seqids -input_type fasta -out $work_directory"/Ref_Blasts/Bsudanica_ref"

  echo "blastp -query "$Bglabrata_prot" -db "$work_directory"/Ref_Blasts/Bsudanica_ref -outfmt \"6 std qcovs\" -out "$work_directory"/Blast_Results/Bglabrata_vs_Bsudanica.blastp -num_threads "'$NPROCS' >> Run_BLast_cluster.sh
  echo "blastp -query "$Bsudanica_prot" -db "$work_directory"/Ref_Blasts/Bglabrata_ref -outfmt \"6 std qcovs\" -out "$work_directory"/Blast_Results/Bsudanica_vs_Bglabrata.blastp -num_threads "'$NPROCS' >> Run_BLast_cluster.sh
fi

if [ $RUN_CROSS_BLAST == TRUE ]
then
  if [ -d $work_directory"/Best_reciprocal_hits" ]
  then
    rm -rf $work_directory"/Best_reciprocal_hits"
  fi
  mkdir $work_directory"/Best_reciprocal_hits"

  current_blast=$work_directory"/Blast_Results/Bglabrata_vs_Bsudanica.blastp"
  out_reciprocal=Bglabrata_vs_Bsudanica
  run_reciprocal_check

  current_blast=$work_directory"/Blast_Results/Bsudanica_vs_Bglabrata.blastp"
  out_reciprocal=Bsudanica_vs_Bglabrata
  run_reciprocal_check

  cat $work_directory"/Best_reciprocal_hits/"*"_hit_to_hit.work" | sort | uniq -c | awk '{if ($1>1) print $2}' >> $work_directory"/Best_reciprocal_hits/Reciprocal_1to1.txt"

fi


if [ $RUN_FINAL_Registry == TRUE ]
then
  out_reciprocal=Bglabrata_vs_Bsudanica
  guide_1=$Bglabrata_prot_ID_guide
  guide_2=$Bsudanica_prot_ID_guide

  if [ -f $work_directory"/Best_reciprocal_hits/Temp_gen_reciprocal.txt" ]
  then
    rm -rf $work_directory"/Best_reciprocal_hits/Temp_gen_reciprocal.txt"
  fi

  prots_on_bglabrata=$(awk -F "\t" '{print $1}' $work_directory"/Best_reciprocal_hits/"$out_reciprocal"_hit_to_hit.txt")
  for prots_1_id in $prots_on_bglabrata
  do
    prot_pairs=$(grep -F "XX"$prots_1_id"XX" $work_directory"/Best_reciprocal_hits/Reciprocal_1to1.txt")
    check_prot=$(echo $prot_pairs | grep -c .)

    if [ $check_prot -gt 0 ]
    then
      for pair in $prot_pairs
      do
        remove=$(echo "XX"$prots_1_id"XX")
        prots_2_id=$(echo $pair | sed "s/$remove//" | sed "s/_xxx_//" | sed "s/^XX//" | sed "s/XX$//")
        echo $pair
        echo $prots_1_id" --- "$prots_2_id

        gen1=$(awk -F "\t" -v prot=$prots_1_id '{if ($2==prot) print $1}' $guide_1 )
        gen2=$(awk -F "\t" -v prot=$prots_2_id '{if ($2==prot) print $1}' $guide_2)
        echo $gen1" "$gen2 | tr " " "\t" >> $work_directory"/Best_reciprocal_hits/Temp_gen_reciprocal.txt"
      done
    fi
  done
  sort -u  $work_directory"/Best_reciprocal_hits/Temp_gen_reciprocal.txt" >  $work_directory"/Best_reciprocal_hits/Gene_to_Gene_IDs.txt"

  awk -F "\t" '{print $1}'  $work_directory"/Best_reciprocal_hits/Gene_to_Gene_IDs.txt" | sort | uniq -u > $work_directory"/Best_reciprocal_hits/unique1.tmp"
  awk -F "\t" '{print $2}'  $work_directory"/Best_reciprocal_hits/Gene_to_Gene_IDs.txt" | sort | uniq -u > $work_directory"/Best_reciprocal_hits/unique2.tmp"

  grep -w -F -f  $work_directory"/Best_reciprocal_hits/unique1.tmp" $work_directory"/Best_reciprocal_hits/Gene_to_Gene_IDs.txt" | grep -w -F -f $work_directory"/Best_reciprocal_hits/unique2.tmp" > $work_directory"/Best_reciprocal_hits/Reciprocal_Genes.txt"

  rm -fr $work_directory"/Best_reciprocal_hits/unique1.tmp"
  rm -fr $work_directory"/Best_reciprocal_hits/unique2.tmp"
  rm -rf $work_directory"/Best_reciprocal_hits/Temp_gen_reciprocal.txt"
fi

################################################################################

# The next bit is hardcoded because I'm angy
# that's right, angy

if [ $RUN_INTERESTING_LINES == "TRUE" ]
then
  Genes_LG5R=Gene_Model_LG5R.txt
  # Genes_LG11L=Gene_Model_LG11L.txt

  genes_chr5_bg=Gene_List_Chr5_Bg.in
  genes_bsud=Gene_List_Bsud.in

  genes_bg_chr5=$(awk -F "\t" '{print $1}' $genes_chr5_bg)
  echo "Bglabrata_GenID Gene_Model Start End Strand Bsucanica_GenID Linkage_Group" | tr " " "\t" > Bglabrata_to_Sudanica.out

  for gen in $genes_bg_chr5
  do
    echo $gen
    # LG5 Gene Gene_Model_LG5R
    coordinates=$(grep -w -F $gen $genes_chr5_bg | awk -F "\t" '{print $2" "$3" "$4}')

    check_gene_model=$(grep -w -c $gen $Genes_LG5R)
    if [ $check_gene_model -gt 0 ]
    then
      gene_Model=$(grep -w $gen $Genes_LG5R | awk -F "\t" '{print $2}')
    else
      gene_Model=NA
    fi

    check_reciprocal_sudanica=$(grep -w -F -c $gen $work_directory"/Best_reciprocal_hits/Reciprocal_Genes.txt")
    if [ $check_reciprocal_sudanica -gt 0 ]
    then
      sudanica_gen=$(grep -w -F $gen $work_directory"/Best_reciprocal_hits/Reciprocal_Genes.txt" | awk -F "\t" '{print $2}')
      linkage_group=$(awk -F "\t" -v gen=$sudanica_gen '{if ($1==gen) print $3}' $genes_bsud | tr "\n" ";" | sed "s/;$//")

      check_link=$(echo $linkage_group |  grep -c . )
      if [ $check_link -eq 0 ]
      then
        linkage_group=$(echo ".")
      fi
    else
      sudanica_gen=NA
      linkage_group=NA
    fi

    echo $gen" "$gene_Model" "$coordinates" "$sudanica_gen" "$linkage_group | tr " " "\t" >> Bglabrata_to_Sudanica.out
  done
fi
