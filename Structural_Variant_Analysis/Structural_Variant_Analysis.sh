#!/usr/bin/env bash

# This script takes the original unscaffolded genomes and aligns first Queries vs References, an then References against each other. Then it
# conducts a basic variant calling with paftools.js. Then it filters the variants found within a specific region, predicts their physiological
# impact given a particular annotation, and then counting by impact level, Gene Model and biological group.

# The quality of this approach degrades quickly with the fragmentation of the reference genome.
# If BS90_bu_et_al is to be used I recomend using one of the scaffolded version relative to the IM line.

# Before running set the following variables:

# 1) Work directory where to store all the information and the cluster submission folder to SLURM
work_directory=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/SNP_Search_on_Assemblies_Raw_Genomes
cluster_submission_folder=$work_directory"/Cluster_Submissions"

# 2) Genome_tables
# query_assemblies=$work_directory"/Query_assemblies_selection.txt"
query_assemblies=$work_directory"/Query_assemblies.txt"
reference_assemblies=$work_directory"/Reference_assemblies_selected.txt"

# 3) Regions of interest
# region_file=$work_directory"/Regions_Interest.txt"
region_file=$work_directory"/Regions_Interest_Extra_peaks.txt"

# 4) Groups of samples to calculate similarities in a group aware manner. For now it is intended to work with the BS90 and 1316 assemblies individually
bio_groups_file=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/SNP_Search_on_Assemblies_Raw_Genomes/Biologically_relevant_groups.txt

# 5) Annotation Refferences to consider
gene_annotation_sources=$(echo "BS90_bu_et_al_qtl iM_bu_et_al_qtl iM_Zhang_2024 Bg_2024")
total_annotations_considered=$(echo $gene_annotation_sources | tr " " "\n" | grep -c .)

# 6) snpEff configuration file
snpEff_config=/nfs5/IB/Blouin_Lab/users/javierc/data/snpEff_special_folder/custom_snpEff.txt

# 7) Annotation GFFs
self_gff_folder=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Annotation_Transfer_Liftoff/Temp_Files/Annotations
transfer_gff_folder=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Annotation_Transfer_Liftoff/Individual_Transfers

# 8) Gene model file
gene_model_file=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Interest_gene_reports/Gene_Model_Input.in

# 9) Summary Count notes
manual_notes_tables=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/SNP_Search_on_Assemblies_Raw_Genomes/Variation_Density_Mk3/Variant_evaluation_table__Manual

#################################################################################################
#################################################################################################

SLURM_MINIMAP2_SNP_RUN_Re=FALSE # maps the genomes assemblies and gets the basic VCF files
SLURM_ANNOT_VCFs=FALSE # Annotates variants and estimate potential impact
REFORMAT_OUTPUT=TRUE # Takes the annotated VCFs, reformats them and makes basic counts
SUMMARY_Variants=FALSE # Summarizes the number of variants per gene model

#################################################################################################
#################################################################################################

minimap_commands () {
  echo "minimap2 -c -x asm5 --cs -t "'$NPROCS'" "$ref_genome" "$query_genome" > "$work_directory"/Minimap_plus_Basecalling/Alignment/"$secondID"_vs_"$firstID"_aln.paf" >> $cluster_submission_folder"/1_"$submit_to_cluster"_minimap.txt"
  echo "sort -k6,6 -k8,8n "$work_directory"/Minimap_plus_Basecalling/Alignment/"$secondID"_vs_"$firstID"_aln.paf | paftools.js call -L 500 -l 100 -f "$ref_genome" -s "$secondID"_vs_"$firstID" - > "$work_directory"/Minimap_plus_Basecalling/VCF_Files/"$secondID"_vs_"$firstID".vcf" >> $cluster_submission_folder"/2_"$submit_to_cluster"_paftools.txt"
}

location_variant () {
  if [ $control_location_data == TRUE ]
  then
    registry_location_variant=$(echo REMOVE_ME)
    registry_gene_model=$(echo REMOVE_ME)

    for annot_source in $gene_annotation_sources
    do
      check_gene=$(awk -F "\t" -v chr=$region_chr -v coord=$coord '{if ($1==chr && $2<=coord && coord<=$3) print}' $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/Genes_from_"$annot_source"_reference_"$referenceID".txt" | grep -c .)
      if [ $check_gene -eq 0 ]
      then
        current_location=Intergenic
      elif [ $check_gene -gt 1 ]
      then
        current_location=Multiple_genes
      else
        gen_ID=$(awk -F "\t" -v chr=$region_chr -v coord=$coord '{if ($1==chr && $2<=coord && coord<=$3) print $4}' $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/Genes_from_"$annot_source"_reference_"$referenceID".txt" | tr ";" "\n" |  grep ^ID= | sed "s/ID=//")
        check_gene_model=$(awk -F "\t" -v gen=$gen_ID -v annot=$annot_source '{if ($4==gen && $5==annot) print $2}' $gene_model_file | grep -c .)

        if [ $check_gene_model -gt 0 ]
        then
          current_gene_model=$(awk -F "\t" -v gen=$gen_ID -v annot=$annot_source '{if ($4==gen && $5==annot) print $2}' $gene_model_file)
          current_gene_status=$(awk -F "\t" -v gen=$gen_ID -v annot=$annot_source '{if ($4==gen && $5==annot) print $3}' $gene_model_file)

          registry_gene_model=$(echo $registry_gene_model" "$current_gene_model)
          current_location=$(echo $gen_ID"("$current_gene_model"__"$current_gene_status")")
        else
          current_location=$(echo $gen_ID"(NA)")
        fi
      fi
      registry_location_variant=$(echo $registry_location_variant" "$current_location)
    done
    registry_location_variant=$(echo $registry_location_variant | sed "s/REMOVE_ME //")
    warning_multiple_gene_models=$(echo $registry_gene_model | tr " " "\n" | sort -u | grep -v -c REMOVE_ME | awk -F "\t" '{if ($1<2) print "FALSE"; else print "TRUE"}')
  fi
  control_location_data=FALSE
}

make_counts () {
  control_other_group=TRUE

  # Start The counts
  N_High=0
  N_Moderate=0
  N_Low=0
  N_Modifier=0

  count=1
  walker=$(grep -c . $work_directory"/Variation_Density_Mk3/Summary_Tables/Current_hits.tmp")

  while [ $count -le $walker ]
  do
    current_info=$(sed -n $count"p" $work_directory"/Variation_Density_Mk3/Summary_Tables/Current_hits.tmp" | tr "\t" ";")

    # Control for manual overrides
    check_issues=$(echo $current_info | awk -F ";" -v manual_override=$manual_override '{print $manual_override}' | grep -c .)
    check_issues_2=$(echo $current_info | awk -F ";" -v manual_override=$manual_override '{print $manual_override}')
    add_to_counts=TRUE

    if [ $check_issues -eq 0 ]
    then
      get_info_cols=$default_effect_cols
    elif [ $check_issues_2 == "NA" ] || [ $check_issues_2 != $model ]
    then
      get_info_cols=$(echo $current_info | tr ";" "\n" | grep -n $model"__" | awk -F ":" -v col_correction=$total_annotations_considered '{print $1+col_correction}')
    else
      get_info_cols=Skip
      add_to_counts=FALSE
      echo $model" "$current_info" Manual_Override" >> $work_directory"/Variation_Density_Mk3/Summary_Tables/Debug_"$region"_Registry.txt"
    fi

    if [ $add_to_counts == TRUE ]
    then
      for info in $get_info_cols
      do
        echo $current_info | awk -F ";" -v col_correction=$info '{print $col_correction}' >> $work_directory"/Variation_Density_Mk3/Summary_Tables/Current_counts.tmp"
      done

      check_low=$(grep -c LOW $work_directory"/Variation_Density_Mk3/Summary_Tables/Current_counts.tmp")
      check_modifier=$(grep -c MODIFIER $work_directory"/Variation_Density_Mk3/Summary_Tables/Current_counts.tmp")
      check_moderate=$(grep -c MODERATE $work_directory"/Variation_Density_Mk3/Summary_Tables/Current_counts.tmp")
      check_high=$(grep -c High $work_directory"/Variation_Density_Mk3/Summary_Tables/Current_counts.tmp")

      if [ $check_high -gt 0 ]
      then
        N_High=$(($N_High + 1))
        echo $model";"$current_info";High;"$count_identifier | tr ";" "\t" >> $work_directory"/Variation_Density_Mk3/Summary_Tables/Debug_"$region"_Registry.txt"
      elif [ $check_moderate -gt 0 ]
      then
        N_Moderate=$(($N_Moderate + 1))
        echo $model";"$current_info";Moderate;"$count_identifier | tr ";" "\t"  >> $work_directory"/Variation_Density_Mk3/Summary_Tables/Debug_"$region"_Registry.txt"
      elif [ $check_low -gt 0 ]
      then
        N_Low=$(($N_Low + 1))
        echo $model";"$current_info";Low;"$count_identifier   | tr ";" "\t" >> $work_directory"/Variation_Density_Mk3/Summary_Tables/Debug_"$region"_Registry.txt"
      elif [ $check_modifier -gt 0 ]
      then
        N_Modifier=$(($N_Modifier + 1))
        echo $model";"$current_info";Modifier;"$count_identifier | tr ";" "\t" >> $work_directory"/Variation_Density_Mk3/Summary_Tables/Debug_"$region"_Registry.txt"
      fi
      rm -f $work_directory"/Variation_Density_Mk3/Summary_Tables/Current_counts.tmp"
    fi

    count=$(($count + 1))
  done

  echo $count_identifier";"$model";"$N_High";"$N_Moderate";"$N_Low";"$N_Modifier
  echo $count_identifier" "$model" "$N_High" "$N_Moderate" "$N_Low" "$N_Modifier | tr " " "\t" >> $work_directory"/Variation_Density_Mk3/Summary_Tables/"$region"_summary.txt"
  rm -f $work_directory"/Variation_Density_Mk3/Summary_Tables/Current_hits.tmp"
}

run_snpeff_shenanigans () {
  for query in $current_assembly_seq
  do
    if [ ! $query == $referenceID ]
    then
      check_already_done=$(grep -w -F -c $regionID"_x_"$referenceID"_x_"$query $work_directory"/Variation_Density_Mk3/Already_Done.txt")
      if [ $check_already_done -eq 0 ]
      then
        echo "Run: "$regionID" --- "$query" bedtools intersect"
        bedtools intersect -header -a $work_directory"/Minimap_plus_Basecalling/VCF_Files/"$query"_vs_"$referenceID".vcf"  -b  $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/region.bed" >> $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/"$query"_vs_"$referenceID"_temp.vcf"
        echo "Run: "$regionID" --- "$query" bcftools filter"
        bcftools filter -o $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/"$query"_vs_"$referenceID".vcf" -e 'QUAL < 60' $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/"$query"_vs_"$referenceID"_temp.vcf"
        rm -f $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/"$query"_vs_"$referenceID"_temp.vcf"
        echo "Run: "$regionID" --- "$query" bgzip"
        bgzip $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/"$query"_vs_"$referenceID".vcf"
        echo "Run: "$regionID" --- "$query" bcftools index"
        bcftools index $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/"$query"_vs_"$referenceID".vcf.gz"

        for anot_source in $gene_annotation_sources
        do
          echo "snpEff -c "$snpEff_config" -v Bg_"$referenceID"_Anot_"$anot_source" "$work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/"$query"_vs_"$referenceID".vcf.gz > "$work_directory"/Variation_Density_Mk3/Annotation_VCF_Files/"$regionID"/"$query"_vs_"$referenceID"_anot_"$anot_source".vcf" >> $cluster_submission_folder"/"$submit_to_cluster"_snpEff.txt"
          echo "mv snpEff_summary.html "$work_directory"/Variation_Density_Mk3/Annotation_VCF_Files/"$regionID"/"$query"_vs_"$referenceID"_anot_"$anot_source"_snpEff_summary.html" >> $cluster_submission_folder"/"$submit_to_cluster"_snpEff.txt"
          echo "mv snpEff_genes.txt "$work_directory"/Variation_Density_Mk3/Annotation_VCF_Files/"$regionID"/"$query"_vs_"$referenceID"_anot_"$anot_source"_snpEff_genes.txt" >> $cluster_submission_folder"/"$submit_to_cluster"_snpEff.txt"
        done

        echo $regionID"_x_"$referenceID"_x_"$query >> $work_directory"/Variation_Density_Mk3/Already_Done.txt"

        chr=$(awk -F "\t" '{print $1}' $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/region.bed")
        start=$(awk -F "\t" '{print $2}' $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/region.bed")
        end=$(awk -F "\t" '{print $3}' $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/region.bed")

        awk -F "\t" -v chr=$chr -v start=$start -v end=$end '{if ($6==chr && $8<=end && start <= $9) print $6" "$8" "$9}' $work_directory"/Minimap_plus_Basecalling/Alignment/"$query"_vs_"$referenceID"_aln.paf" >> $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/Map_Coverage/"$query"_vs_"$referenceID"_mapped_sections.paf"
      fi
    fi
  done
}

################################################################################
################################################################################

if [ $SLURM_MINIMAP2_SNP_RUN_Re == "TRUE" ]
then
  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "SNP_CALL_"$current_date)

  if [ ! -d $work_directory"/Minimap_plus_Basecalling" ]
  then
    mkdir $work_directory"/Minimap_plus_Basecalling"
    mkdir $work_directory"/Minimap_plus_Basecalling/Alignment"
    mkdir $work_directory"/Minimap_plus_Basecalling/VCF_Files"

    echo "Already Processed:" >> $work_directory"/Minimap_plus_Basecalling/Already_Done.txt"
  fi

  reference_IDs=$(awk -F "\t" '{print $1}' $reference_assemblies | sed 1d)
  query_IDs=$(awk -F "\t" '{print $1}' $query_assemblies | sed 1d)

  for refID in $reference_IDs
  do
    firstID=$refID
    ref_genome=$(awk -F "\t" -v genomeID=$refID '{if ($1==genomeID) print $2}' $reference_assemblies)

    for queryID in $query_IDs
    do
      echo $queryID" vs "$refID
      check_done=$(grep -w -F -c $queryID"__"$refID $work_directory"/Minimap_plus_Basecalling/Already_Done.txt")
      if [ $check_done -eq 0 ]
      then
        query_genome=$(awk -F "\t" -v genomeID=$queryID '{if ($1==genomeID) print $2}' $query_assemblies)
        secondID=$queryID
        minimap_commands
        echo $queryID"__"$refID >> $work_directory"/Minimap_plus_Basecalling/Already_Done.txt"
      fi
    done
  done

  for Ref1 in $reference_IDs
  do
    for Ref2 in $reference_IDs
    do
      echo $Ref1" vs "$Ref2

      if [ ! $Ref1 == $Ref2 ]
      then
        check_done=$(grep -w -F -c $Ref1"__"$Ref2 $work_directory"/Minimap_plus_Basecalling/Already_Done.txt")
        if [ $check_done -eq 0 ]
        then
          firstID=$Ref1
          secondID=$Ref2

          ref_genome=$(awk -F "\t" -v genomeID=$Ref1 '{if ($1==genomeID) print $2}' $reference_assemblies)
          query_genome=$(awk -F "\t" -v genomeID=$Ref2 '{if ($1==genomeID) print $2}' $reference_assemblies)

          minimap_commands
          echo $Ref1"__"$Ref2>> $work_directory"/Minimap_plus_Basecalling/Already_Done.txt"
        fi
      fi
    done
  done
fi

################################################################################

if [ $SLURM_ANNOT_VCFs == "TRUE" ]
then
  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "VCF_Tools_"$current_date)

  region_IDs=$(sed 1d $region_file |awk -F "\t" '{print $2}' |  sort -u)

  if [ ! -d $work_directory"/Variation_Density_Mk3" ]
  then
    mkdir $work_directory"/Variation_Density_Mk3"
    mkdir $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files"
    mkdir $work_directory"/Variation_Density_Mk3/Annotation_VCF_Files"

    echo "Already Processed:" >> $work_directory"/Variation_Density_Mk3/Already_Done.txt"
  fi

  for regionID in $region_IDs
  do
    if [ ! -d $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID ]
    then
      mkdir $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID
      mkdir $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/Map_Coverage"
    fi

    if [ ! -d $work_directory"/Variation_Density_Mk3/Annotation_VCF_Files/"$regionID ]
    then
      mkdir $work_directory"/Variation_Density_Mk3/Annotation_VCF_Files/"$regionID
    fi

    awk -F "\t" -v region=$regionID '{if ($2==region) print $3" "$4" "$5}' $region_file | tr " " "\t" >> $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/region.bed"

    referenceID=$(awk -F "\t" -v region=$regionID '{if ($2==region) print $1}' $region_file)
    current_assembly_seq=$(awk -F "\t" '{print $1}' $query_assemblies | sed 1d)
    run_snpeff_shenanigans

    current_assembly_seq=$(awk -F "\t" '{print $1}' $reference_assemblies | sed 1d)
    run_snpeff_shenanigans
  done
fi

################################################################################

if [ $REFORMAT_OUTPUT == TRUE ]
then
  if [ ! -d $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw" ]
  then
    mkdir $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw"
  fi

  all_count_regionsIDs=$(sed 1d $bio_groups_file | awk -F "\t" '{print $5}')
  for count_regionID in $all_count_regionsIDs
  do
    echo "Working on... "$count_regionID
    regionID=$(awk -F "\t" -v count=$count_regionID '{if ($5==count) print $1}' $bio_groups_file)

    # DELETE ME!!!
    if [ -d $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID ]
    then
      rm -rf $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID
    fi

    if [ ! -d $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID ]
    then
      mkdir $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID
      mkdir $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files"

      N_Variant=1
      referenceID=$(awk -F "\t" -v region=$count_regionID '{if (region == $5) print $2}' $bio_groups_file)
      control_sus=$(awk -F "\t" -v region=$count_regionID '{if (region == $5) print $3}' $bio_groups_file | tr ";" " ")
      control_res=$(awk -F "\t" -v region=$count_regionID '{if (region == $5) print $4}' $bio_groups_file | tr ";" " ")

      region_chr=$(awk -F "\t" '{print $1}' $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/region.bed")
      region_start=$(awk -F "\t" '{print $2}' $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/region.bed")
      region_end=$(awk -F "\t" '{print $3}' $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/region.bed")

      total_sus=$(echo $control_sus | tr " " "\n" | grep -c .)
      total_res=$(echo $control_res | tr " " "\n" | grep -c .)

      all_assemblies=$(echo $control_sus" "$control_res | tr " " "\n" | sort)
      header_queries_variant=$(echo $all_assemblies | tr " " "\n" | awk '{print "Var_"$1}' |  tr "\n" " " | sed "s/ $/\n/")
      header_location_anot=$(echo $gene_annotation_sources | tr " " "\n" | awk '{print "Location_relative_to_"$1}' |  tr "\n" " " | sed "s/ $/\n/")
      header_effect_anot=$(echo $gene_annotation_sources | tr " " "\n" | awk '{print "Effect_on_"$1}' |  tr "\n" " " | sed "s/ $/\n/")

      echo "N_Variant Coord Ref Query "$header_queries_variant" "$header_location_anot" "$header_effect_anot" N_GrA_("$total_sus") N_GrB_("$total_res") Is_All_GrA Is_All_GrB All_samples Warn_Multiple_GeneModels Num_NA_GrA Num_NA_GrB Num_Ref_hits_GrA Num_Ref_hits_GrB" | tr " " "\t" >> $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Result_Table.txt"
      echo "Coord "$all_assemblies | tr " " "\t" >> $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Result_Variant_Heatmap.txt"

      echo "Preparing temporary files"
      # 0) Create smaller temp files
        # -> Variants in Query
      temp_anot_source=$(echo $gene_annotation_sources | awk -F " " '{print $1}')
      for query in $all_assemblies
      do
        echo "Variant Registry for: "$query
        grep -v "#" $work_directory"/Variation_Density_Mk3/Annotation_VCF_Files/"$regionID"/"$query"_vs_"$referenceID"_anot_"$temp_anot_source".vcf" | awk -F "\t" '{print $1"__"$2"__"$4"__"$5}' >> $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/"$query"_vs_"$referenceID"_variants.tmp"
      done

      cat $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/"*"_vs_"$referenceID"_variants.tmp" | sort -u | sed "s/__/\t/g" | sort -n -k2 | sed "s/\t/__/g" >> $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/All_variants.tmp"
      cat  $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/All_variants.tmp" | awk -F "__" '{print "__"$2"__"}' | uniq -d >> $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/Polymorphic_sites.tmp"
      all_individual_rariants=$(grep -v -f $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/Polymorphic_sites.tmp" $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/All_variants.tmp")

      rm -f $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/temp_all_variants.tmp"
      echo "Done Variant Registry"

      # Gene_locations
      echo "Preparing annotation..."
      for annot_source in $gene_annotation_sources
      do
        echo "Gene locations on... "$annot_source
        if [ $annot_source == $referenceID ]
        then
          gff_file=$self_gff_folder"/"$annot_source"_annotation.gff "
        else
          gff_file=$transfer_gff_folder"/"$annot_source"_transfered_to_"$referenceID"_annotation.gff_polished"
        fi
        awk -F "\t" -v chr=$region_chr -v start=$region_start -v end=$region_end '{if ($1==chr && $3=="gene" && $4<=end && start<=$5) print $1"\t"$4"\t"$5"\t"$9}' $gff_file >> $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/Genes_from_"$annot_source"_reference_"$referenceID".txt"
      done

      echo "Done with preparations. Now to do the work."
      for variant in $all_individual_rariants
      do
        echo "Working on... "$count_regionID" --- "$variant
        coord=$(echo $variant | awk -F "__" '{print $2}')
        ref_version=$(echo $variant | awk -F "__" '{print $3}')
        query_version=$(echo $variant | awk -F "__" '{print $4}')

        num_sus=0
        num_res=0

        NA_hits_sus=0
        NA_hits_res=0

        Ref_hits_sus=0
        Ref_hits_res=0

        poly_hits_sus=0
        poly_hits_res=0

        # 1) Check overlapping genes
        control_location_data=TRUE
        control_evaluation_data=TRUE
        registry_query_report=REMOVE_ME

        for query in $all_assemblies
        do
          #2) Check if variant present in the sample
          check_sample_group_sus=$(echo $control_sus | grep -w -F -c $query)
          check_sample_group_res=$(echo $control_res | grep -w -F -c $query)

          check_SNP_on_query=$(grep -w -F -c $variant $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/"$query"_vs_"$referenceID"_variants.tmp")
          # -> No: Check if its in alignment
              # -> Yes: Assume equal to reference
              # -> No: print NA
          if [ $check_SNP_on_query -eq 0 ]
          then
            check_SNP_covered=$(awk -F " " -v coord=$coord '{ if (coord >= $2 && coord <= $3) print}' $work_directory"/Variation_Density_Mk3/Filtered_VCF_Files/"$regionID"/Map_Coverage/"$query"_vs_"$referenceID"_mapped_sections.paf" | grep -c .)
            if [ $check_SNP_covered -eq 0 ]
            then
              NA_hits_sus=$(($NA_hits_sus + $check_sample_group_sus))
              NA_hits_res=$(($NA_hits_res + $check_sample_group_res))

              report_query=NA
            else
              report_query=$(echo ".")
              Ref_hits_sus=$(($Ref_hits_sus + $check_sample_group_sus))
              Ref_hits_res=$(($Ref_hits_res + $check_sample_group_res))
              location_variant
            fi
          else
            # 3) Count toward the right group
            report_query=X

            num_sus=$(($num_sus + $check_sample_group_sus))
            num_res=$(($num_res + $check_sample_group_res))

            location_variant

            # 4) Check for each annotation the effect on at least one gene
            if [ $control_evaluation_data == TRUE ]
            then
              registry_effect_per_anot=REMOVE_ME
              for annot_source in $gene_annotation_sources
              do
                awk -F "\t" -v chr=$region_chr -v coord=$coord -v ref_v=$ref_version -v query_v=$query_version '{if ($1==chr && $2==coord && $4==ref_v && $5==query_v) print $8}' $work_directory"/Variation_Density_Mk3/Annotation_VCF_Files/"$regionID"/"$query"_vs_"$referenceID"_anot_"$annot_source".vcf" | tr ";" "\n" | grep ^ANN= > $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/current_variant.tmp"

                check_warnings=$(grep -c WARNING_ $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/current_variant.tmp")
                check_low=$(grep -c \|LOW\| $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/current_variant.tmp")
                check_modifier=$(grep -c \|MODIFIER\| $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/current_variant.tmp")
                check_moderate=$(grep -c \|MODERATE\| $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/current_variant.tmp")
                check_high=$(grep -c \|HIGH\| $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Temp_Files/current_variant.tmp")

                if [ $check_warnings -gt 0 ]
                then
                  current_effect=Warn
                elif [ $check_high -gt 0 ]
                then
                  current_effect=High
                elif [ $check_moderate -gt 0 ]
                then
                  current_effect=MODERATE
                elif [ $check_low -gt 0 ]
                then
                  current_effect=LOW
                elif [ $check_modifier -gt 0 ]
                then
                  current_effect=MODIFIER
                fi
                registry_effect_per_anot=$(echo $registry_effect_per_anot" "$current_effect)
              done
              registry_effect_per_anot=$(echo $registry_effect_per_anot | sed "s/REMOVE_ME //")
              control_evaluation_data=FALSE
            fi
          fi
          registry_query_report=$(echo $registry_query_report" "$report_query)
        done
        registry_query_report=$(echo $registry_query_report | sed "s/REMOVE_ME //")


        heatmap_registry=$(echo $registry_query_report | sed "s/\./0/g" | sed "s/X/1/g")
        echo $coord" "$heatmap_registry | tr " " "\t" >> $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Result_Variant_Heatmap.txt"

        # check if it meets exclusivility groups
        All_sus=$(echo $num_sus" "$total_sus | awk -F " " '{if ($1==$2) print "TRUE"; else print "FALSE"}')
        All_res=$(echo $num_res" "$total_res | awk -F " " '{if ($1==$2) print "TRUE"; else print "FALSE"}')
        All_samples=$(echo $All_sus" "$All_res | awk -F " " '{if ($1==$2 && $1=="TRUE") print "TRUE"; else print "FALSE"}')

        echo $N_Variant" "$coord" "$ref_version" "$query_version" "$registry_query_report" "$registry_location_variant" "$registry_effect_per_anot" "$num_sus" "$num_res" "$All_sus" "$All_res" "$All_samples" "$warning_multiple_gene_models" "$NA_hits_sus" "$NA_hits_res" "$Ref_hits_sus" "$Ref_hits_res | tr " " "\t"  >> $work_directory"/Variation_Density_Mk3/Variant_evaluation_table__Raw/"$count_regionID"/Result_Table.txt"
        N_Variant=$(($N_Variant + 1))
      done
    fi
  done
fi

################################################################################

if [ $SUMMARY_Variants == TRUE ]
then
  if [ -d $work_directory"/Variation_Density_Mk3/Summary_Tables" ]
  then
    rm -rf $work_directory"/Variation_Density_Mk3/Summary_Tables"
  fi
  mkdir $work_directory"/Variation_Density_Mk3/Summary_Tables"

  summarize_regions=$(ls $manual_notes_tables"/"*"_manual_checks.txt" | sed "s/.*\///" | sed "s/_manual_checks.txt//")
  for region in $summarize_regions
  do
    column_All_Sus=$(head -n 1 $manual_notes_tables"/"$region"_manual_checks.txt" | grep -n "Is_All_GrA" | awk -F ":" '{print $1}')
    column_All_Res=$(head -n 1 $manual_notes_tables"/"$region"_manual_checks.txt" | grep -n "Is_All_GrB" | awk -F ":" '{print $1}')

    count_sus=$(head -n 1 $manual_notes_tables"/"$region"_manual_checks.txt" | grep -n "N_GrA_(" | awk -F ":" '{print $1}')
    count_res=$(head -n 1 $manual_notes_tables"/"$region"_manual_checks.txt" | grep -n "N_GrB_(" | awk -F ":" '{print $1}')

    NA_sus=$(head -n 1 $manual_notes_tables"/"$region"_manual_checks.txt" | grep -n "Num_NA_GrA" | awk -F ":" '{print $1}')
    NA_res=$(head -n 1 $manual_notes_tables"/"$region"_manual_checks.txt" | grep -n "Num_NA_GrB" | awk -F ":" '{print $1}')

    hit_ref_sus=$(head -n 1 $manual_notes_tables"/"$region"_manual_checks.txt" | grep -n "Num_Ref_hits_GrA" | awk -F ":" '{print $1}')
    hit_ref_res=$(head -n 1 $manual_notes_tables"/"$region"_manual_checks.txt" | grep -n "Num_Ref_hits_GrB" | awk -F ":" '{print $1}')
    manual_override=$(head -n 1 $manual_notes_tables"/"$region"_manual_checks.txt" | grep -n "Manual_Override" | awk -F ":" '{print $1}')

    echo "Working on: "$region
    echo "Count_Group Gene_Model N_High N_Moderate N_Low N_Modifier" | tr " " "\t" >> $work_directory"/Variation_Density_Mk3/Summary_Tables/"$region"_summary.txt"
    get_regionID=$(awk -F "\t" -v region=$region '{if ($2==region) print $6}' $region_file)
    gene_model_IDs=$(awk -F "\t" -v region=$get_regionID '{if ($1==region) print $2}' $gene_model_file | sort -u)
    default_effect_cols=$(head -n 1 $manual_notes_tables"/"$region"_manual_checks.txt" | tr "\t" "\n" | grep -n Effect_on_ | awk -F ":" '{print $1}')

    control_sus=$(awk -F "\t" -v region=$region '{if (region == $1) print $3}' $bio_groups_file | tr ";" " ")
    control_res=$(awk -F "\t" -v region=$region '{if (region == $1) print $4}' $bio_groups_file | tr ";" " ")
    Nsus=$(echo $control_sus | tr " " "\n" | grep -c .)
    Nres=$(echo $control_res | tr " " "\n" | grep -c .)

    treshold=$(awk -F "\t" -v region=$region '{if (region == $5) print $6}' $bio_groups_file)
    is_fraction=$(echo $treshold | awk '{if ($1>=1) print "FALSE"; else print "TRUE"}')

    for model in $gene_model_IDs
    do

      echo $model
      # count_identifier=Sus
      count_identifier=Pigmented
      if [ $is_fraction == "TRUE" ]
      then
        work_tresh=$( echo $Nsus" "$treshold | awk -F " " '{printf "%.0f\n", $1*$2}')
      else
        work_tresh=$treshold
      fi
      grep "("$model"__" $manual_notes_tables"/"$region"_manual_checks.txt" | awk -F "\t" -v treshold=$work_tresh -v check_col=$count_sus -v avoid_col=$count_res -v inner_hits_ref=$hit_ref_sus '{if ($check_col>=treshold && $avoid_col==0 && $inner_hits_ref==0) print}' |   awk -F "\t" -v other_na=$NA_res -v other_hit_ref=$hit_ref_res '{if ($other_na > 1 || $other_hit_ref >1) print}' > $work_directory"/Variation_Density_Mk3/Summary_Tables/Current_hits.tmp"
      make_counts

      # count_identifier=Res
      count_identifier=Albino
      if [ $is_fraction == "TRUE" ]
      then
        work_tresh=$( echo $Nres" "$treshold | awk -F " " '{printf "%.0f\n", $1*$2}')
      else
        work_tresh=$treshold
      fi

      grep "("$model"__" $manual_notes_tables"/"$region"_manual_checks.txt" | awk -F "\t" -v treshold=$work_tresh -v check_col=$count_res -v avoid_col=$count_sus -v inner_hits_ref=$hit_ref_res '{if ($check_col>=treshold && $avoid_col==0 && $inner_hits_ref==0) print}' | awk -F "\t" -v other_na=$NA_sus -v other_hit_ref=$hit_ref_sus '{if ($other_na > 1 || $other_hit_ref >1) print}' > $work_directory"/Variation_Density_Mk3/Summary_Tables/Current_hits.tmp"
      make_counts

      count_identifier=All_Samples
      grep "("$model"__" $manual_notes_tables"/"$region"_manual_checks.txt" | awk -F "\t" -v check_col=$column_All_Res -v check_col_2=$column_All_Sus '{if ($check_col=="TRUE" && $check_col_2=="TRUE") print}' > $work_directory"/Variation_Density_Mk3/Summary_Tables/Current_hits.tmp"
      make_counts
    done
  done
fi
