#!/usr/bin/env bash

# Access to the cluster is achieved with "hpcman", it works similar to samtools
# but with three levels. For submissions there is an alias "hqsub" to submit batch jobs.

# There are many ways to submit comands to the cluster but the easier way seems
# to be the use of the bash interpreter to execute a file, and then subbmit them by a pipe.

#  bash submit_jobs.sh | hqsub -q '*' -r - -p 8 -t array -

# -q: indicates the node to use
# -r: provided runmane. "-" makes an automatic name
# -p: number of requested
# -t: type of job. Batch (one line) or Array (multiple)

################################################################################

# This script takes assemblies classified as references and Blouin_work located
# on different folders and maps them to one another in order to find out interest
# regions on the references.
# In order to make use of the cluster infraestructure, this script does not runs
# most of the analysis. Instead it creates files for array jobs and then sorts
# out the output. When possible/advisable.

# To run it you need to edit the content of the following variables
# 1) Indicate where output files are to be located, batch files for the cluster.
# The folder must exist previosly

work_folder=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation
cluster_submission_folder=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Cluster_Submissions

# 2) Set path to a work table indicating names annotations for every genome file
# Because output names will be based on the fields provided it is important that:
#    a) No special characters. Just vanilla english letters and "_", "-" or "."
#    b) Cols tab separated
#    c) NO WHITE SPACES ANYWHERE!!! The script will let you know you mess up.
#
# The rable must have 4 columns:
#   - ID: An unique ID that identifies the assembly
#   - R/Q: Indication if it is a reference assembly (R) or one of the lab's assemblies (Q)
#   - Location: Full path to genome file (asumed to be compresed)
#   - Annotation: Full path to annotation files in GFF if available. NA otherwise
#
# Only assemblies on this table will be considered. The script will abort if
# one is missing.
# On its current implementation, the script assumes that the genomes are not compressed.
# The script also assumes good beheavior and that the table first line is a header.
# it does not uses that line, it can be anything and even blank. But the script
# ignores it at severa points

genome_table=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Input_assemblies.in

# 3) Table specifying the initial search for the region of interest.

# This must be a 7 column table separated with tabs:
# 1:Region - ID for the region of interest.
# 2:Assembly - ID for the assembly (must exist on Genome table)
# 3:Contig-Scaffold - ID for the coting-scaffold-chromosome where the region was identified
# 4:Anchor_point - Single base position that is confirmed to be close to the start of the region
# 5:Anchor_location - Expected location of the anchor point "Start", "Middle" and "End". The extension will be applied relative to this information
# 6:Extension - Aproximate extension of the region of interest
# 7:Margins - Desired error margins around the regions to consider
# 8:Primer_File - Full path to primer file with primers. Write "NA" if unavailable.

# The script also assumes good beheavior and that the table first line is a header.
# it does not uses that line, it can be anything and even blank. But the script
# ignores it at severa points

region_table=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Wanted_region_Final_LG5_And_10.in

# 3) Pick what block of code you are going to run right now.
# This script is modular, each one of the following variables makes a different
# set of commands. Everything needs to be run but not everything needs to be done
# right now. The separation allows for modifications of the script and partial
# re-runs. To control this set the relecant control variables, written in all
# caps, to "TRUE", anything else like "42", "bazinga" or even "true" will not
# count.
# Varuables that start with "SLURM" don't do actual analysis, they generate submision
# files for SLURM. They work a little differently than SLURM's documentation in
# the name of simplicity since they are not scripts, they are a text file with the
# commands to run.  To load them to the cluster use cat:

# cat [File to upload] | hpcman queue submit - -q '*' -t array -p [Nº Processors] --watch -r [Name for run] -f [Requested Memory]G

# Furthermore, given the high chance of future expansion of the study to more
# assemblies, many of these blocks can be run wither from scratch (any previous
# output is removed) or updated (only missing samples are processed). By default
# the script assumes you want to add on already existing runs. But you can change
# this by setting the variable CLEAN_BURN to "TRUE", this will remove all the
# saved for that block of code. Except the submission commands.

# Blocks of code with this feature end their name with _Re

# 4) Final region selection and processing.
# Once you have all the information requested and the data explored. The final step before annotating the genes themselves
# is to manually select and trim the regions that are going to be studied. This is going to gbe based mostly on the presence of gene markers
# in the target area. Corrected by the synteny blocks in cases where missasemblies are suspected to have happened

# The table provided here needs to have X fields

final_region_recovery=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Selected_Regions_input/Selected_Regions_Additional_peaks_Color_from_iM_Zhang_2024.in

# This is a tab delimited table with the following fields:
  # 1) RegionID: ID of the region
  # 2) Assembly_ID: Assembly ID
  # 3) Sequence_File: Full path for the genome file to use
  # 4) Subregions: Location of the refined region to analyze. It must follow the format: [SequenceID]":"[Start Coord]"-"[End Coord]":"[Strand symbol +/-]. Eg: iM_Zhang_2024_CM074214.1_RagTag:39536346-43254379:+
  # Multiple fragments can be set up toghether here by separating entries with ";"
  # 5) Is_Scaffolded: Boolean indication if the genome used is scaffolded by this pipeline or not.
  # 6) Scaffold_Reference: ID of the reference assembly used for scaffolding. Set to "NA" is none was used

# The file name itself needs to follow the format: [Something]"_from_"[Ref ID where the region was defined]".in". This is used to specify later if a gene signature was found in the original assembly or not

# 5) File with instructions on how to extract the loci ID from the protein sequence from a GFF. Necessary annoyance to account for different formats.
annotation_keystone=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/GFF_gene_ID_keystone.in
  # Assembly                Prefix   Subfix  Take
  # Bg_2024                 ID=cds-  ;       gene=
  # BS90_bu_et_al_qtl       ID=cds-  ;       locus_tag=
  # iM_bu_et_al_qtl         ID=cds-  ;       locus_tag=
  # iM_Zhang_2024           ID=cds-  ;       locus_tag=

# For this Bg_2024 the script will search for ID=cds-[PROT_ID]; and take the ID found on gene=

# 6) It was undersuded in our study, but if you have a combination combination of Interpro IDs that identified to be part of the same family
interpro_domain_combo=/nfs5/IB/Blouin_Lab/users/javierc/ANALYSIS/Key_Region_Anotation/Dotplot_annotation/Domain_Combos.in

# Note: Transmembrane domains are requested so much that they get their own display

# Here is an example:
# Combo_Name      Domains
# Collagen_TNF    IPR001073;IPR050392;IPR008983
# Miosin_Kinesin  IPR001609;IPR027417;IPR036072;IPR036961;IPR010926
# SOCS_box        IPR001496;IPR036036

# 7) folder paths toward the things needed to run run Interproscan through singularity locally
interpro_data=/nfs5/IB/Blouin_Lab/users/javierc/data/Interpro_local_data/interproscan-5.74-105.0/data
sif_image_path=/nfs5/IB/Blouin_Lab/users/javierc/data/Interpro_local_data

################################################################################
################################################################################
################################################################################

CLEAN_BURN=FALSE

# Validation Blocks
VALIDATION_AND_STATS=FALSE # Verify that all the inputs are good and some stats of the fasta files
VALIDATION_REGION=FALSE # Verify that the region table works

# Repeat Masking
SLURM_REPEATMODELER=FALSE # Runs Repeat modeler over the best reference according to its N50
SLURM_REPEATMASKER_Re=FALSE # Runs Repeatmasker to mask and annotate repeats
REFORMAT_REPEATMASKER=FALSE # Reformats the main table of Repeatmasker into something easier to work with

# RagTag Scaffolding and synteny
SLURM_RAGTAG_SCAFFOLDING_Re=FALSE # Runs RAGTAG on the query assemblies vs all the references
SLURM_RAGTAG_SCAFFOLDING_REFERENCES_Re=FALSE # Runs RAGTAG between the references
SLURM_RAGTAG_SUMARY=FALSE # Makes a summary of RAGTAG Results
SLURM_RAGTAG_MINIMAP_SYRI_REGION_Re=FALSE # Takes the scaffold/contig/chr and calculates sunteny regions
REGION_FINDER_Re=FALSE # Track down the synteny blocks overlapping the desired region (the entire chromosome)
FULL_ANCHOR_BLAST=FALSE # Emergency control of quality. Improve to generalize or delete
SLURM_SYRI_PLOTSR_Re=FALSE # Generates the synteny plots for the regions of study

# Protein Annotation Transfer
UNMASKED_SCAFFOLDING_Re=FALSE # Takes Ragstag's output and produces the same scaffolding but with no masking
SLURM_ANNOTATION_TRANSFER_Re=FALSE # Uses Liftoff to transfer annotations from each reference to the RAGTAG Scaffolding.
CONFIRM_LIFTOFF_RUNS=FALSE # Utility block. It check no liftoff command was accidentally skipped
EVALUATE_ANOT_TRANSFER=FALSE # Measures how well each assembly transfer did.

# Syn Block annotation summary
RAW_SYN_BLOCK_TARGET_REGIONS_ANNOTATION=FALSE # Checks all the sources of information and reports on the features found on every sample
RAW_SYN_BLOCK_TARGET_REGIONS_ADDITIONAL_TRACKS=FALSE # Get GFF tracks for the location of each contig location in the scaffolded genome.
SLURM_RAW_SYN_BLOCK_TARGET_REGIONS_ADDITIONAL_TRACKS_SYRI=FALSE # Get GFF tracks for the location of each contig location in the scaffolded genome.
SLURM_RAW_SYN_BLOCK_TARGET_REGIONS_ADDITIONAL_TRACKS_GFF=FALSE # Reformat syri outputs for the reference as GFF file tracks

# Gene annotation
SLURM_GENE_SEQ_INTERPRO_Re=FALSE # Annotates the given protein files for each reference
SLURM_GENE_SEQ_EGGNOG_Re=FALSE # Annotates the given protein files for each reference
SLURM_INTERPRO_PER_GENE_Re=FALSE # This takes the genes in the GFF and adds to interpro's outputs the gene loci code
SLURM_EGGNOG_PER_GENE_Re=FALSE # This takes the genes in the GFF and adds to eggnog's outputs the gene loci code
INTERPRO_PER_GENE_INDEX=FALSE # Generates an index file to guide the retrivery of information
SLURM_SIGNALP6_Re=FALSE # Runs Signalp6 on the annotated genes.
SLURM_TMHMM_Re=FALSE # Identifies transmembrane domains tmhmm
LOCAL_TARGETP=FALSE # I haven't had confirmation from the CQLS team if they plan to do anything with TARGETP. So I'm gonna do a Thanos and do it myself. The online server accepts up to 5,000 sequences so this script devides the reference sequences. The next block will sort them out.
REFORMAT_TARGETP_OUTPUT=FALSE # Combines the outputs run on the web server into a single table and adds the Gene ID on every protein entry

# Final region extraction
REGION_GENE_EXTRACTION_Mk2=TRUE # This block code retrieves the final sequences locations for further study. You need to provide it the table final_region_recovery
REGION_BASIC_SUMMARY=TRUE # Generates a basic summary of the gene counts per gene block
REGION_GENE_ANOTATION_INTERPRO_Re=TRUE # Takes the gene IDs found in the region and summarises the annotation for each gene.
REGION_GENE_ANOTATION_SIGNALP_TARGETP_TMHMM=TRUE # Retrieves the location signals and transmembrane domain counts for the genes within the region
REGION_GENE_ANOTATION_TRACKS=TRUE # Generates filtered gff files for in depth review on IGV viewer.

# Each block assumes that the previopus one ran to completion.
# Do not skip them and if you modify the script make sure the output of each one
# remains consistent. Same place, same file format and same naming convention.
# Only Exception is VALIDATION_AND_STATS. It's only job is to check everything
# is in order. It can be run at any time.

# 4) Other variables
# Here are several varables that are used thought the script but do not affect
# the results.

threads=12 # number of procesors to be used by each tool on this computer
anchor_blast_ident=85 # Identity used to check for the presence/absence of anchors
anchor_blast_qcov=85 # qcovs used to check for the presence/absence of anchors
override_repeatmasker=FALSE # this is an indication to not use the masked genomes for the analysis and use instead the ones provided on genome_table directly. These are assumed to not be masked but they don't have to be. If you want to use masked genomes from other pipeline set this variable to true.
override_repeatmodeler=FALSE # this is an indication to not use a different library for the tandem repeats. To use this feature indicate a path to the new library

################################################################################
################################################################################
################################################################################

# Versions used:
# Minimap2: 2.28-r1209
# RagTag: v2.1.0
# Samtools: 1.20 (using htslib 1.20)
# BLAST: 2.15.0
# Syri: 1.7.0
# Liftoff: 1.6.3

################################################################################
################################################################################
################################################################################

# Functions
# these are little bits of code that are used more than once throught the script

confirm_input_file () {
  # Confirmation that the work table is functional.
  if [ ! -f $genome_table ]
  then
    # The file does not exits...
    echo "Table with the genome data missing. Check this: "
    echo "genome_table="$genome_table
    exit
  else
    # Is there something bad with the genome table?
    check_white_spaces=$(grep -c " " $genome_table)
    check_entries=$(grep -c . $genome_table | awk '{print $1-1}')
    check_N_R=$(awk -F "\t" '{print $2}' $genome_table | sed 1d | grep -w -c R)
    check_N_Q=$(awk -F "\t" '{print $2}' $genome_table | sed 1d | grep -w -c Q)
    check_N_other=$(awk -F "\t" '{print $2}' $genome_table | sed 1d | grep -w -v Q| grep -w -v R | grep -c .)
    check_dup_ID=$(awk -F "\t" '{print $1}' $genome_table | sed 1d | sort | uniq -d | grep -c .)

    if [ $check_white_spaces -gt 0 ] || [ $check_N_R -eq 0 ] || [ $check_N_Q -eq 0 ] || [ $check_dup_ID -gt 0 ] || [ $check_N_other -gt 0 ] || [ $check_entries -le 0 ]
    then
      # There is something bad with the table -.-"
      echo "Something is wrong with: "$genome_table
      echo "N entries: "$check_entries" -- Remember, the script assumes the first line is a header and excludes it."
      echo "N white spaces: "$check_white_spaces" (must be 0)"
      echo "N of reference assemblies: "$check_N_R" (must be greater than 0)"
      echo "N of Query assemblies: "$check_N_Q" (must be greater than 0)"
      echo "N of other assemblies: "$check_N_other" (must be 0)"
      echo "N duplicated IDs: "$check_dup_ID" (must be 0)"
      exit
    fi
  fi

  if [ ! -f $annotation_keystone ]
  then
    echo "No keystone file detected. Consider how do you plan to obtain the loci IDs from the protein IDs on the GFF"
    exit
  fi

  if [ ! -f $interpro_domain_combo ]
  then
    echo "No file with Interpro domains found. These are necesary near the end to facilitate the exploration of the target regions."
  fi
}

confirm_region_file () {
  if [ ! -f $region_table ]
  then
    # The file does not exits...
    echo "Table with the wanted regions are missing. Check this: "
    echo "region_table="$region_table
    exit
  else
    check_entries=$(grep -c . $region_table | awk '{print $1-1}')
    check_white_spaces=$(grep -c " " $region_table)
    check_dup_ID=$(awk -F "\t" '{print $1}' $region_table | sed 1d | sort | uniq -d | grep -c .)
    check_wrong_anchor_points=$(awk -F "\t" '{print $5}' $region_table | sed 1d | grep -w -v Start | grep -w -v Middle | grep -w -v End | grep -c .)

    if [ $check_white_spaces -gt 0 ] || [ $check_dup_ID -gt 0 ] || [ $check_entries -le 0 ] || [ $check_wrong_anchor_points  -gt 0 ]
    then
      echo "Something is wrong with: "$region_table
      echo "N entries: "$check_entries" -- Remember, the script assumes the first line is a header and excludes it."
      echo "N white spaces: "$check_white_spaces" (must be 0)"
      echo "N duplicated IDs: "$check_dup_ID" (must be 0)"
      echo "N Wrong Anchor points: "$check_wrong_anchor_points" (must be 0, every field can be Start, Middle, or End)"
      exit
    fi
  fi
}

start_block_code () {
  # preparation for the run.
  if [ $CLEAN_BURN == "TRUE" ]
  then
    rm -rf $work_folder"/"$storage
  fi

  if [ ! -d $work_folder"/"$storage ]
  then
    mkdir $work_folder"/"$storage
    echo "Already completed..." > $work_folder"/"$storage"/Already_done.txt"
  fi
}

check_anchor_primers () {
  check_duplicated_IDs=$(seqkit seq -n -i $primer_file_location | sort | uniq -d | grep -c .)
  check_primer_IDs=$(seqkit seq -n -i $primer_file_location | sort | grep -c .)
  check_IDs_with_metadata=$(seqkit seq -n $primer_file_location | awk -F " " '{print $2}' | grep -c .)

  if [ $check_duplicated_IDs -gt 0 ] || [ ! $check_primer_IDs == $check_IDs_with_metadata ] || [ $check_primer_IDs -eq 0 ]
  then
    echo "Something is wrong with the primer file"
    echo "Check the following:"
    echo "N Duplicated IDs: "$check_duplicated_IDs" (should be zero)"
    echo "N Primer IDs: "$check_primer_IDs" of which "$check_IDs_with_metadata" have position data. Should be equal and non-zero"
    exit
  fi
}

code_block_completed () {
  # The code block is completed. What to do now.
  echo ""
  echo "Code block completed!!"
  echo "If something needs adjusting but you don't want to repeat a full run remove"
  echo "all files with the specific genome pair and then edit file 'Already_done.txt' for this run"
  echo ""
}

slurm_block_completed () {
  # Useful message about how to submit the commads to SLURM
  echo "This piece of code generated submissions files for SLURM"
  echo "You can find them here: "$cluster_submission_folder
  echo "to run them use this: "
  echo "cat [File to upload] | hpcman queue submit - -q '*' -t array -p [Nº Processors] --watch -r [Name for run] -f [Requested Memory]G"
  echo "..."
  echo "Remember to replace the fraking text between '[ ]' with whatever names or numbers you need!!!"
  echo '(-_-)""'
  echo ""
}

index_fasta () {
  if [ ! -f $work_folder"/"$storage"/Temp_file_sequences/"$indexID"_temp_copy.fasta.fai" ]
  then
    echo "Indexing: "$indexID
    cp $index_genome $work_folder"/"$storage"/Temp_file_sequences/"$indexID"_temp_copy.fasta"
    rename=$(echo ">"$indexID"_")
    sed -i "s/>/$rename/" $work_folder"/"$storage"/Temp_file_sequences/"$indexID"_temp_copy.fasta"
    samtools faidx $work_folder"/"$storage"/Temp_file_sequences/"$indexID"_temp_copy.fasta"
  fi
}

remove_ragtag_previous_run () {
# Ragtag tries to recycle intermediary files by default, but I don't know how it handles it
# therefore while I could set the comand to overwrite output files I'm siding with caution.
# if there is an output folder for an alignment that is not on the control file I want it gone.

if [ -d $work_folder"/"$storage"/Ref_"$refID1"_vs_Ref_"$refID2"_ragtag_output" ]
then
  rm -rf $work_folder"/"$storage"/Ref_"$refID1"_vs_Ref_"$refID2"_ragtag_output"
fi
}

ragtag_summary_start () {
  if [ -f $work_folder"/"$storage"/Run_summary.txt" ]
  then
    rm -f $work_folder"/"$storage"/Run_summary.txt"
  fi

  if [ -d $work_folder"/"$storage"/Temp_file_sequences" ]
  then
    rm -rf $work_folder"/"$storage"/Temp_file_sequences"
  fi

  ragtag_registry_missing=NA
  ragtag_registry_error=NA

  echo "Scaffolding_ID Error Missing N_placed_sequences N_placed_bp N_unplaced_sequences N_unplaced_bp N_gap_bp N_gap_sequences Placed_sequences% Placed_bases%" | tr " " "\t" >> $work_folder"/"$storage"/Run_summary.txt"

  echo "Retrieving all relevant information from the RagTag outputs"
  echo "Group: "$comparison_group
}

recover_information_ragtag_runs () {
  # This is a bunch of data retrieva questions intended to make a table that summarises
  # how the scaffolding runs went.
  if [ ! -d $ragtag_output_path"/"$ragtag_output_folder ]
  then
    echo "Missing output for: "$ragtag_output_path"/"$ragtag_output_folder
    missing=X
    stat_data=$(echo ". . . . . . . .")
    ragtag_registry_missing=$(echo $ragtag_registry_missing" "$print_id)
  else
    check_error=$(grep -c . $ragtag_output_path"/"$ragtag_output_folder"/ragtag.scaffold.err")
    if [ $check_error -gt 0 ]
    then
      error=X
      stat_data=$(echo ". . . . . . . .")
      ragtag_registry_missing=$(echo $ragtag_registry_missing" "$print_id)
    else
      stat_data=$(tail -n 1 $ragtag_output_path"/"$ragtag_output_folder"/ragtag.scaffold.stats")
      per_placed_seq=$(tail -n 1 $ragtag_output_path"/"$ragtag_output_folder"/ragtag.scaffold.stats" | awk -F "\t" '{print ($1/($1+$3))*100"%"}')
      per_placed_bases=$(tail -n 1 $ragtag_output_path"/"$ragtag_output_folder"/ragtag.scaffold.stats" | awk -F "\t" '{print ($2/($2+$4))*100"%"}')
    fi
  fi

  echo $print_id" "$error" "$missing" "$stat_data" "$per_placed_seq" "$per_placed_bases | tr " " "\t" >> $work_folder"/"$storage"/Run_summary.txt"
}

print_issues_ragtag () {
  check_issues=$(echo $registry | tr " " "\n" | grep -v -w -c NA)
  if [ $check_issues -gt 0 ]
  then
    echo ""
    echo "There are errors with "$check_issues" results. Reason: "$source
    echo "Check the following analysis: "
    echo $registry | sed "s/^NA //" | tr " " "\n"
    echo ""
  fi
}

fix_out_of_bounds_coords () {
  chr_lenght=$(awk -F "\t" -v genome=$genome_ID '{if ($1==genome) print}' $work_folder"/Assembly_Stats/Sequence_stats.txt"| grep -w -F $chr_ID | awk -F "\t" '{print $3}')
  fix_start=$(echo $fix_coord | awk -F "," '{if ($1 < 1) print 1"_change"; else print $1}')
  fix_end=$(echo $fix_coord","$chr_lenght | awk -F "," '{if ($2 > $3) print $3"_change"; else print $2}')

  check_adjustment=$(echo $fix_start" "$fix_end | grep -c "_change")
  if [ $check_adjustment -gt 0 ]
  then
    region_coord_adjusted=X
  fi
  fix_coord=$(echo $fix_start","$fix_end | sed "s/_change//g")
}

current_subregion_overlap_check () {
  result_overlap=$(echo $query_coords","$target_coords | awk -F "," '{if ($1 <= $4) print}' | awk -F "," '{if ($3 <= $2) print}' | grep -c . | awk '{if ($1>0) print "X"; else print "."}')
}

do_ragtag_syri () {
  # Exclude sequences from references
  # SLURM Commands

  echo $store_result_ID" "$reference_fasta_minimap" "$target_fasta_minimap | tr " " "\t" >> $local_work_folder"/Run_Registrty.txt"

  check_ref1=$(grep -c . $reference_fasta_minimap)
  check_ref2=$(grep -c . $target_fasta_minimap)

  if [ $check_ref1 -gt 0 ] && [ $check_ref2 -gt 0 ]
  then
    echo "minimap2 --eqx -x asm5 -t "'$NPROCS'" -a "$reference_fasta_minimap" "$target_fasta_minimap" > "$local_work_folder"/Minimap/"$store_result_ID"_scaffold_aln.sam" >> $cluster_submission_folder"/1_"$submit_to_cluster"_minimap.txt"
    echo "samtools sort -@ "'$NPROCS'" -O BAM -o "$local_work_folder"/Minimap/"$store_result_ID"_scaffold_aln.bam" $local_work_folder"/Minimap/"$store_result_ID"_scaffold_aln.sam" >> $cluster_submission_folder"/2_"$submit_to_cluster"_sort_samtools.txt"
    echo "samtools index -@ "'$NPROCS'" "$local_work_folder"/Minimap/"$store_result_ID"_scaffold_aln.bam" >> $cluster_submission_folder"/3_"$submit_to_cluster"_index_samtools.txt"
    echo "syri -F B -c "$local_work_folder"/Minimap/"$store_result_ID"_scaffold_aln.bam -r "$reference_fasta_minimap" -q "$target_fasta_minimap" --dir "$local_work_folder"/Syri --prefix "$store_result_ID"_ --nc "'$NPROCS' >> $cluster_submission_folder"/4_"$submit_to_cluster"_syri.txt"
  fi
}

query_to_ref_correction () {
  if [ $ragtag_strand == "+" ]
  then
    result_start=$(($fix_start + $correction))
    result_end=$(($fix_end + $correction))
  elif [ $ragtag_strand == "-" ]
  then
    full_lenght=$(grep -w -F $genome_ID"_"$chr_ID $ragtag_output"/ragtag.scaffold.agp" | awk -F "\t" '{print $8}')
    result_start=$(($full_lenght - $fix_end + $correction + 1))
    result_end=$(($full_lenght - $fix_start + $correction + 1))
  else
    echo "Am I a joke to you?!"
    echo $ragtag_strand
    exit
  fi
}

extract_syn_block_from_ragtag () {
  if [ $extract_sec_from_ragtag == TRUE ]
  then
    ragtag_start_coord=$(awk -F "\t" -v block=$block '{if ($9==block) print $7}' $block_file)
    ragtag_end_coord=$(awk -F "\t" -v block=$block '{if ($9==block) print $8}' $block_file)

    seqkit grep -p $search_ragtag_id $ragtag_output"/ragtag.scaffold.fasta" | seqkit subseq -r $ragtag_start_coord":"$ragtag_end_coord | sed "s/>.*/$rename_seq/" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Sequences/"$genome_ID"_vs_"$targetID"_synteny_blocks.fasta"
    matching_original_contigs=$(awk -F "\t" -v seq=$search_ragtag_id '{if ($1==seq && $5!="U" && $5!="N") print}' $ragtag_output"/ragtag.scaffold.agp" | awk -F "\t" -v end=$ragtag_end_coord '{if ($2 <= end) print}' |  awk -F "\t" -v start=$ragtag_start_coord '{if (start <= $3) print $6}' | sort -u | sed "s/^$targetID//" | sed "s/^_//" | tr "\n" ";" | sed "s/;$//" )

    N_matching_original_contigs=$(echo $matching_original_contigs | tr ";" "\n" | grep -c .)

    if [ $N_matching_original_contigs -eq 0 ]
    then
      matching_original_contigs=Missing
    fi
  else
    extract_seq=$(echo $search_ragtag_id | sed "s/$targetID//" | sed "s/^_//" | sed "s/_RagTag$//")
    ragtag_start_coord=$(awk -F "\t" -v block=$block '{if ($9==block) print $2}' $block_file)
    ragtag_end_coord=$(awk -F "\t" -v block=$block '{if ($9==block) print $3}' $block_file)

    seqkit grep -p $extract_seq $target_genome_file | seqkit subseq -r $ragtag_start_coord":"$ragtag_end_coord | sed "s/>.*/$rename_seq/" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Sequences/"$genome_ID"_vs_"$targetID"_synteny_blocks.fasta"
    matching_original_contigs=$(awk -F "\t" -v seq=$search_ragtag_id '{if ($1==seq && $5!="U" && $5!="N") print}' $ragtag_output"/ragtag.scaffold.agp" | awk -F "\t" -v end=$ragtag_end_coord '{if ($2 <= end) print}' |  awk -F "\t" -v start=$ragtag_start_coord '{if (start <= $3) print $1}' | sort -u | sed "s/^$targetID//" | sed "s/^_//" | tr "\n" ";" | sed "s/;$//" | sed "s/_RagTag$//")
    N_matching_original_contigs=$(echo $matching_original_contigs | tr ";" "\n" | grep -c .)
  fi
}

process_genes_off_target () {
  count_off_target=1
  walker_off_targe=$(grep -c . $missing_genes_file)

  grep -f $missing_genes_file $annotation_file >> $work_folder"/"$storage"/annotation_reduced.tmp"
  echo "Reference GFF File: "$annotation_file >> $missing_genes_file".info"
  echo "" >> $missing_genes_file".info"
  echo "GenID Expected_Chr//Found_Chr Status" >> $missing_genes_file".info"

  while [ $count_off_target -le $walker_off_targe ]
  do
    off_id=$(sed -n $count_off_target"p" $missing_genes_file)
    awk -F "\t" -v chr=$target_chromosome '{if ($3=="gene") print $1" "$9}' $work_folder"/"$storage"/annotation_reduced.tmp" | grep -F -w "ID="$off_id >> $work_folder"/"$storage"/Off_target_gene.tmp"

    check_found=$(grep -c . $work_folder"/"$storage"/Off_target_gene.tmp")
    if [ $check_found -eq 0 ]
    then
      status=Not_Found
      off_target_chr=NA
    else
      off_target_chr=$(awk -F " " '{print $1}' $work_folder"/"$storage"/Off_target_gene.tmp" | sort -u | tr "\n" ";" | sed "s/;$//")
      if [ $off_target_chr == $target_chromosome ]
      then
        status=off_target
      else
        status=Transposition_candidate
      fi
    fi

    echo $off_id" "$target_chromosome"//"$off_target_chr" "$status >> $missing_genes_file".info"

    rm -rf $work_folder"/"$storage"/Off_target_gene.tmp"
    count_off_target=$(($count_off_target + 1))
  done
  rm -f $work_folder"/"$storage"/annotation_reduced.tmp"
}

get_prefix_info_for_gff () {
  check_info=$(grep -w -F -c $current_ref $annotation_keystone)
  if [ $check_info -eq 0 ]
  then
    echo $current_ref" is not on the instruction file:"
    echo $annotation_keystone
    exit
  else
    prefix_extract=$(grep -w -F $current_ref $annotation_keystone | awk -F "\t" '{print $2}')
    subfix_extract=$(grep -w -F $current_ref $annotation_keystone | awk -F "\t" '{print $3}')
    take_item=$(grep -w -F $current_ref $annotation_keystone | awk -F "\t" '{print $4}')
  fi
}

extract_syntenic_blocks () {
  if [ ! -f $syri_output ]
  then
    echo $targetID" Missing NA NA NA NA NA NA NA NA NA" | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.txt"
  else
    awk -F "\t" '{if ($10=="-") print}' $syri_output | awk -F "\t" -v data_start=$start_col -v search_end=$wanted_end '{if ($data_start <= search_end) print}' |  awk -F "\t" -v data_end=$end_col -v search_start=$wanted_start '{if (search_start <= $data_end) print}' >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks.txt"

    synteny_blocks=$(awk -F "\t" '{print $9}' $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks.txt")
    for block in $synteny_blocks
    do
      # Annotation of structural variants
      check_query_hit=$(awk -F "\t" -v block=$block -v queryid=$Notal_in_query '{if ($9==block && $queryid=="-") print}' $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks.txt" | grep -c .)
      rename_seq=$(echo ">"$targetID"_"$block)

      if [ $check_query_hit -eq 0 ]
      then
        # RAGTAG information
        block_file=$work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks.txt"
        extract_syn_block_from_ragtag
      else
        ragtag_start_coord=NA
        ragtag_end_coord=NA

        matching_original_contigs=NA
        N_matching_original_contigs=0
      fi

      # Overlap with regions
      query_coords=$(awk -F "\t" -v block=$block -v data_start=$start_col -v data_end=$end_col '{if ($9==block) print $data_start","$data_end}' $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks.txt")
      block_coords=$(echo $query_coords | tr "," " ")

      if [ ! $region_lower_boundary == "NA" ]
      then
        target_coords=$region_lower_boundary
        current_subregion_overlap_check
        overlap_lower_boundary=$result_overlap

        target_coords=$region_upper_boundary
        current_subregion_overlap_check
        overlap_upper_boundary=$result_overlap
      else
        current_subregion_overlap_check=NA
        overlap_upper_boundary=NA
      fi

      target_coords=$region_Start_End
      current_subregion_overlap_check
      overlap_target=$result_overlap

      echo $targetID" "$block" "$block_coords" "$ragtag_start_coord" "$ragtag_end_coord" "$N_matching_original_contigs" "$matching_original_contigs" "$overlap_lower_boundary" "$overlap_target" "$overlap_upper_boundary | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.tmp"
    done

    awk -F "\t" '{print $9}' $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks.txt" > $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_asigned_blocks.tmp"

    if [ $extract_sec_from_ragtag == TRUE ]
    then
      check_missing_blocks_start=$(awk -F "\t" '{print $7}' $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks.txt" | tr -d "-" | sort -n |grep . | head -n 1)
      check_missing_blocks_end=$(awk -F "\t" '{print $8}' $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks.txt" | tr -d "-" | sort -n |grep . | tail -n 1)

      awk -F "\t" '{if ($10=="-") print}' $syri_output | awk -F "\t" -v search_start=$check_missing_blocks_start -v search_end=$check_missing_blocks_end '{if ($7 >= search_start && $8 <= search_end ) print}' | grep -v -w -F -f $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_asigned_blocks.tmp" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks_miss_on_ref.txt"
      rm -f $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_asigned_blocks.tmp"
      miss_synteny_blocks=$(awk -F "\t" '{print $9}' $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks_miss_on_ref.txt")

      for block in $miss_synteny_blocks
      do
        rename_seq=$(echo ">"$targetID"_M_"$block)
        block_file=$work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks_miss_on_ref.txt"
        extract_syn_block_from_ragtag

        query_coords=$(awk -F "\t" -v block=$block -v data_start=$start_col -v data_end=$end_col '{if ($9==block) print $data_start","$data_end}' $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks_miss_on_ref.txt")
        block_coords=$(echo $query_coords | tr "," " ")

        overlap_lower_boundary=Skip
        overlap_target=Skip
        overlap_upper_boundary=Skip

        echo $targetID" M_"$block" "$block_coords" "$ragtag_start_coord" "$ragtag_end_coord" "$N_matching_original_contigs" "$matching_original_contigs" "$overlap_lower_boundary" "$overlap_target" "$overlap_upper_boundary | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.tmp"
      done
    else
      extract_seq=$(echo $search_ragtag_id | sed "s/$targetID//" | sed "s/^_//" | sed "s/_RagTag$//")

      check_missing_blocks_start=$(awk -F "\t" '{print $2}' $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks.txt" | tr -d "-" | sort -n |grep . | head -n 1)
      check_missing_blocks_end=$(awk -F "\t" '{print $3}' $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks.txt" | tr -d "-" | sort -n |grep . | tail -n 1)
      awk -F "\t" '{if ($10=="-") print}' $syri_output | awk -F "\t" -v search_start=$check_missing_blocks_start -v search_end=$check_missing_blocks_end '{if ($2 >= search_start && $3 <= search_end ) print}' | grep -v -w -F -f $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_asigned_blocks.tmp" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks_miss_on_ref.txt"
      rm -f $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_asigned_blocks.tmp"

      miss_synteny_blocks=$(awk -F "\t" '{print $9}' $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks_miss_on_ref.txt")

      for block in $miss_synteny_blocks
      do
        rename_seq=$(echo ">"$targetID"_Miss_"$block)
        block_file=$work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks_miss_on_ref.txt"
        extract_syn_block_from_ragtag

        ragtag_start_coord=$(awk -F "\t" -v block=$block '{if ($9==block) print $7}' $block_file)
        ragtag_end_coord=$(awk -F "\t" -v block=$block '{if ($9==block) print $8}' $block_file)

        query_coords=$(awk -F "\t" -v block=$block -v data_start=2 -v data_end=3 '{if ($9==block) print $data_start","$data_end}' $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/"$targetID"_Raw_Syri_blocks_miss_on_ref.txt")
        block_coords=$(echo $query_coords | tr "," " ")

        overlap_lower_boundary=Skip
        overlap_target=Skip
        overlap_upper_boundary=Skip

        echo $targetID" M_"$block" "$ragtag_start_coord" "$ragtag_end_coord" "$block_coords" "$N_matching_original_contigs" "$matching_original_contigs" "$overlap_lower_boundary" "$overlap_target" "$overlap_upper_boundary | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.tmp"
      done
    fi

    if [ -f $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Sequences/"$genome_ID"_vs_"$targetID"_synteny_blocks.fasta" ]
    then
      N_Blocks_On_Target=$(seqkit seq -n $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Sequences/"$genome_ID"_vs_"$targetID"_synteny_blocks.fasta" | grep -v -c "_Miss_" )
      N_Blocks_Off_Target=$(seqkit seq -n $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Sequences/"$genome_ID"_vs_"$targetID"_synteny_blocks.fasta" | grep -c "_Miss_" )
    else
      N_Blocks_On_Target=NA
      N_Blocks_Off_Target=NA
    fi

    echo $targetID" "$check_missing_blocks_start" "$check_missing_blocks_end" "$N_Blocks_On_Target" "$N_Blocks_Off_Target | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Regions.txt"

    ##################################################################

    if [ $Known_Anchor_Fasta == "NA" ]
    then
      cat $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.tmp" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.txt"
    elif [ -f $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Sequences/"$genome_ID"_vs_"$targetID"_synteny_blocks.fasta" ]
    then
      echo "Regions identified and picked. Now for the anchor BLAST"
      echo "Making BLAST Reference"
      makeblastdb -in $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Sequences/"$genome_ID"_vs_"$targetID"_synteny_blocks.fasta" -dbtype nucl -parse_seqids -input_type fasta -out $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/References/"$targetID"_ref"

      echo "" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/"$targetID"_full_anchor_blastn.tmp"
      if [ $check_short_blast -gt 0 ]
      then
        echo "Short BLAST..."
        blastn -task blastn-short -query $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/Temp_File_short.fasta" -db $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/References/"$targetID"_ref" -outfmt "6 std" -out $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/"$targetID"_short_anchor.blastn" -perc_identity $anchor_blast_ident -qcov_hsp_perc $anchor_blast_qcov -num_threads $threads

        cat $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/"$targetID"_short_anchor.blastn" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/"$targetID"_full_anchor_blastn.tmp"
        echo "Done"
      fi

      if [ $check_long_blast -gt 0 ]
      then
        echo "Long BLAST"
        blastn -query $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/Temp_File_long.fasta" -db $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/References/"$targetID"_ref" -outfmt "6 std" -out $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/"$targetID"_long_anchor.blastn" -perc_identity $anchor_blast_ident -qcov_hsp_perc $anchor_blast_qcov -num_threads $threads
        cat $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/"$targetID"_short_anchor.blastn" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/"$targetID"_full_anchor_blastn.tmp"
        echo "Done"
      fi

      anchor_sequences=$(seqkit seq -n -i $Known_Anchor_Fasta)
      traverse_hits=$(grep -c . $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.tmp")
      count_hits=1

      while [ $count_hits -le $traverse_hits ]
      do
        target_info=$(sed -n $count_hits"p" $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.tmp")
        target_seq=$(sed -n $count_hits"p" $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.tmp" | awk -F "\t" -v tar=$targetID '{print tar"_"$2}')

        registry_blast=$(echo "REMOVE")
        for anchor in $anchor_sequences
        do
          check_result=$(grep -w -F "$anchor" $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/"$targetID"_full_anchor_blastn.tmp" | grep -P "$target_seq" | grep -c .)
          registry_blast=$(echo $registry_blast" "$check_result)
        done

        registry_blast=$(echo $registry_blast | sed "s/^REMOVE //")
        echo $target_info" "$registry_blast | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.txt"

        count_hits=$(($count_hits + 1))
      done
      rm -f $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/"$targetID"_full_anchor_blastn.tmp"
    fi
  fi
  rm -f $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.tmp"
}

write_plotsr_files () {
  # 1) Write Genome file
  echo $reference_file" "$original_genome_ID" lw:1.5" | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Genome_Files/"$original_genome_ID"_vs_"$target_id"_genome_file.txt"
  echo $target_file" "$target_id" lw:1.5" | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Genome_Files/"$original_genome_ID"_vs_"$target_id"_genome_file.txt"

  # 2) Write track file
  echo $chr_name" "$region_coords | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Marker_Locations/"$original_genome_ID"_vs_"$target_id"_target_region.bed"
  echo $work_folder"/"$storage"/"$region_ID"/Marker_Locations/"$original_genome_ID"_vs_"$target_id"_target_region.bed "$region_ID" bw:10000;nc:black;ns:6;lc:blue;lw:1;bc:lightblue;ba:0.5"  | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Marker_Locations/"$original_genome_ID"_vs_"$target_id"_track.txt"

  # 3) Write SLURM Commands
  echo "plotsr --sr "$plot_this_syri" --tracks "$work_folder"/"$storage"/"$region_ID"/Marker_Locations/"$original_genome_ID"_vs_"$target_id"_track.txt --genomes "$work_folder"/"$storage"/"$region_ID"/Genome_Files/"$original_genome_ID"_vs_"$target_id"_genome_file.txt -o "$work_folder"/"$storage"/"$region_ID"/"$original_genome_ID"_vs_"$target_id"_full_plotsr.svg" >> $cluster_submission_folder"/"$submit_to_cluster"_syn_plot.txt"
}

non_masked_scaffolding () {
  target_genome_file=$(awk -F "\t" -v genome=$targetID '{if (genome==$1) print $3}' $genome_table)
  # Cleaning any previous interrupted run
  if [ -d $work_folder"/"$storage"/Temp_files" ]
  then
    rm -rf $work_folder"/"$storage"/Temp_files"
  fi

  if [ -f $new_out_file".tmp" ]
  then
    rm -f $new_out_file".tmp"
  fi

  if [ -f $new_out_file".fasta" ]
  then
    rm -f $new_out_file".fasta"
  fi
  # Done cleaning

  mkdir $work_folder"/"$storage"/Temp_files"

  echo "Spliting "$targetID" into neat parts..."
  seqkit split -s 1 --quiet $target_genome_file -O $work_folder"/"$storage"/Temp_files"
  echo "Doing a registry file"
  grep ">" $work_folder"/"$storage"/Temp_files/"* | sed "s/.*\///" | tr -d ">" | tr ":" "\t" >> $work_folder"/"$storage"/Temp_files/Seq_locations.tmp"

  current_ID=SKIP_ME
  count=3
  walker=$(grep -c . $ragtag_output_folder"/ragtag.scaffold.agp")
  while [ $count -le $walker ]
  do
    position_type=$(sed -n $count"p" $ragtag_output_folder"/ragtag.scaffold.agp" | awk -F "\t" '{print $5}')

    if [ $position_type == "U" ] || [ $position_type == "N" ]
    then
      printf %100s | tr " " "N" >> $new_out_file".tmp"
      echo $check_current_ID" --- GAP ("$count")" >> $new_out_file".tracker"
    else
      check_current_ID=$(sed -n $count"p" $ragtag_output_folder"/ragtag.scaffold.agp" | awk -F "\t" '{print $1}')
      get_scaffold=$(sed -n $count"p" $ragtag_output_folder"/ragtag.scaffold.agp" | awk -F "\t" '{print $6}' | sed "s/^$targetID//" | sed "s/^_//")
      get_orientation=$(sed -n $count"p" $ragtag_output_folder"/ragtag.scaffold.agp" | awk -F "\t" '{print $9}')
      get_seq_file=$(grep -w -F $get_scaffold $work_folder"/"$storage"/Temp_files/Seq_locations.tmp" | awk -F "\t" '{print $1}')

      echo $targetID": " $check_current_ID" -- "$get_scaffold" || Entry: "$count" of "$walker
      if [ ! $check_current_ID == $current_ID ]
      then
        rename=$(echo ">"$check_current_ID)
        if [ $get_orientation == "+" ]
        then
          cat $work_folder"/"$storage"/Temp_files/"$get_seq_file | sed "s/>.*/$rename/" >> $new_out_file".tmp"
          echo $check_current_ID" --- "$get_scaffold" "$get_orientation" ("$count") Start" >> $new_out_file".tracker"
        elif [ $get_orientation == "-" ]
        then
          seqkit seq -u -r -p -t DNA -v $work_folder"/"$storage"/Temp_files/"$get_seq_file  | sed "s/>.*/$rename/" >> $new_out_file".tmp"
          echo $check_current_ID" --- "$get_scaffold" "$get_orientation" ("$count") Start" >> $new_out_file".tracker"
        else
          echo "I hate my job... good thing The Sojourn Audiodrama is a thing."
          echo "By the way, somehow the RAGTAG output is compromised"
          exit
        fi
        current_ID=$check_current_ID
      else
        if [ $get_orientation == "+" ]
        then
          cat $work_folder"/"$storage"/Temp_files/"$get_seq_file | grep -v ">" >> $new_out_file".tmp"
          echo $check_current_ID" --- "$get_scaffold" "$get_orientation" ("$count") Middle" >> $new_out_file".tracker"
        elif [ $get_orientation == "-" ]
        then
          seqkit seq -s -u -r -p -t DNA -v $work_folder"/"$storage"/Temp_files/"$get_seq_file  >> $new_out_file".tmp"
          echo $check_current_ID" --- "$get_scaffold" "$get_orientation" ("$count") Middle" >> $new_out_file".tracker"
        else
          echo "I hate my job... good thing The Sojourn Audiodrama is a thing."
          echo "By the way, somehow the RAGTAG output is compromised"
          exit
        fi
      fi
    fi
    count=$(($count + 1))
  done
  seqkit seq -u -w 60 $new_out_file".tmp" >> $new_out_file".fasta"
  rm -f $new_out_file".tmp"
  rm -rf $work_folder"/"$storage"/Temp_files"
  echo "seqkit locate -j "'$NPROCS'" -f "$target_genome_file" "$new_out_file".fasta >> "$new_out_file".quality control" >>  $cluster_submission_folder"/"$submit_to_cluster"_seqkit_location.txt"
}

prepare_for_liftoff () {
  for RefID in $reference_assemblies
  do
    if [ $really_first == TRUE ]
    then
      N_run_fule=1
    fi

    ref_genome=$(awk -F "\t" -v genome=$RefID '{if (genome==$1) print $3}' $genome_table)
    gff_file=$(awk -F "\t" -v genome=$RefID  '{if (genome==$1) print $4}' $genome_table )

    if [ ! -f $work_folder"/"$storage"/Temp_Files/Sequences/"$RefID"_genome.fasta" ]
    then
      echo "Copying and indexing ref genome"
      cp $ref_genome $work_folder"/"$storage"/Temp_Files/Sequences/"$RefID"_genome.fasta"
      samtools faidx $work_folder"/"$storage"/Temp_Files/Sequences/"$RefID"_genome.fasta"
    fi

    if [ ! $gff_file == "NA"  ]
    then
      target_files=$(ls $target_locations"/"*".fasta" | sed "s/.*\///" | sed "s/_unmasked_scaffolds.fasta//")

      if [ ! -f $work_folder"/"$storage"/Temp_Files/Annotations/Round_1/"$RefID"_annotation.gff" ]
      then
        echo "Copying annotation file"
        local_count=1
        while [ $local_count -le 5 ]
        do
          cp $gff_file  $work_folder"/"$storage"/Temp_Files/Annotations/Round_"$local_count"/"$RefID"_annotation.gff"
          local_count=$(($local_count +1))
        done
      fi

      for TargetID in $target_files
      do
        check_done=$(grep -w -F -c $RefID"_transfered_to_"$TargetID $work_folder"/"$storage"/Already_done.txt")
        if [ $check_done -gt 0 ]
        then
          echo $TargetID" allready done. Skipping"
        else
          echo "Transfering: "$RefID" to "$TargetID

          if [ ! -f $target_locations"/"$TargetID"_unmasked_scaffolds.fasta.fai" ]
          then
            echo "Making Query index"
            samtools faidx $target_locations"/"$TargetID"_unmasked_scaffolds.fasta"
          fi

          if [ $N_run_fule == 1 ]
          then
            liftoff_cluster_submissions=$cluster_submission_folder"/Liftoff_runs/Initial_run"
            echo "liftoff -p "'$NPROCS'" -dir "$work_folder"/"$storage"/Temp_Files/"$RefID"transfer_to_"$TargetID"_temp_files -o "$work_folder"/"$storage"/Individual_Transfers/"$RefID"_transfered_to_"$TargetID"_annotation.gff -u "$work_folder"/"$storage"/Individual_Transfers/"$RefID"_transfered_to_"$TargetID"_unnmaped_features.gff -g "$work_folder"/"$storage"/Temp_Files/Annotations/Round_"$round_storage"/"$RefID"_annotation.gff -copies -polish -cds "$target_locations"/"$TargetID"_unmasked_scaffolds.fasta "$work_folder"/"$storage"/Temp_Files/Sequences/"$RefID"_genome.fasta" >> $liftoff_cluster_submissions"/Slurm_Command.sh"
          else
            liftoff_cluster_submissions=$cluster_submission_folder"/Liftoff_runs/Bulk_run"
            echo "liftoff -p "'$NPROCS'" -dir "$work_folder"/"$storage"/Temp_Files/"$RefID"transfer_to_"$TargetID"_temp_files -o "$work_folder"/"$storage"/Individual_Transfers/"$RefID"_transfered_to_"$TargetID"_annotation.gff -u "$work_folder"/"$storage"/Individual_Transfers/"$RefID"_transfered_to_"$TargetID"_unnmaped_features.gff -g "$work_folder"/"$storage"/Temp_Files/Annotations/Round_"$round_storage"/"$RefID"_annotation.gff -copies -polish -cds "$target_locations"/"$TargetID"_unmasked_scaffolds.fasta "$work_folder"/"$storage"/Temp_Files/Sequences/"$RefID"_genome.fasta" >> $liftoff_cluster_submissions"/Slurm_Command_batch_"$cluster_batch_job".sh"

            round_storage=$(($round_storage + 1))
            if [ $round_storage -gt 5 ]
            then
              round_storage=1
              cluster_batch_job=$(($cluster_batch_job + 1))
            fi
          fi

          echo $RefID"_transfered_to_"$TargetID >> $work_folder"/"$storage"/Already_done.txt"
          N_run_fule=2
        fi
      done
    else
      echo $RefID" doesn't has a GFF file. Verify if true. Skipping"
    fi
  done
}

summary_repeat () {
  N_Repeats=$(grep -c . $repeat_storage)

  if [ $N_Repeats -gt 0 ]
  then
    all_repeats=$(awk -F "\t" '{print $11}' $repeat_storage | sort -u)
    Bases_Repeats=$(awk -F "\t" '{print $7-$6+1}' $repeat_storage | awk -F "\t" '{ sum += $1 } END {print sum}')
    for rep in $all_repeats
    do
      if [ ! -f $repeats_to_summary ]
      then
        echo "Repeat Total Coverend_bases" >> $repeats_summary
      fi

      total_repeats=$(awk -F "\t" -v rep=$rep '{if (rep==$11) print}' $repeat_storage | grep -c .)
      covered_bases=$(awk -F "\t" -v rep=$rep '{if (rep==$11) print $7-$6+1}' $repeat_storage | awk -F "\t" '{ sum += $1 } END {print sum}')
      echo $rep" "$total_repeats" "$covered_bases | tr " " "\t" >> $repeats_summary
    done
  else
    Bases_Repeats=0
  fi
}

repeat_recovery () {
  if [ $target_group == "Q" ]
  then
    ragtag_output=$work_folder"/"$query_group_ragtag"/Ref_"$original_genome_ID"/"$TargetID"_ragtag_output/ragtag.scaffold.agp"
  else
    ragtag_output=$(ls -d $work_folder"/"$reference_group_ragtag"/Ref_"* | grep "_"$original_genome_ID"_" | grep "_"$TargetID"_" | awk '{print $0"/ragtag.scaffold.agp"}')
  fi

  if [ $check_origin_region == "Scaffolded" ]
  then
    ragtag_scaffolded_seq=$(echo $original_genome_ID"_"$original_chr"_RagTag") # Depende de quien sea la referencia en la comparación.
    awk -F "\t" -v seq=$ragtag_scaffolded_seq '{if ($1==seq && $5!="U" && $5!="N") print}' $ragtag_output | awk -F "\t" -v start=$Syn_start -v end=$Syn_end '{if ($2<=end && start<=$3 ) print $6" "$2" "$3" "$7" "$8" "$9}' | tr " " "\t"  | sed "s/$TargetID//" | sed "s/^_//" >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_location.temp"
    awk -F "\t" '{print $1}' $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_location.temp" >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_extract.temp"
  else
    ragtag_scaffolded_seq=$(awk -F "\t" -v get_seq=$original_genome_ID"_"$original_chr '{if ($6==get_seq) print $1}' $ragtag_output  | sed "s/^$TargetID//" |  sed "s/^_//" | sed "s/_RagTag//") # Depende de quien sea la referencia en la comparación.
    echo $ragtag_scaffolded_seq" "$Syn_start" "$Syn_end" "$Syn_start" "$Syn_end" Direct" | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_location.temp"
    echo $ragtag_scaffolded_seq>> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_extract.temp"
  fi

  # The repeats files are very big and scan them one contig at the time takes a while. This reduces the computing task.
  repeats_file=$(ls $work_folder"/"$Repeat_information"/"$TargetID"/"*".out.tab")
  grep -w -F -f $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_extract.temp" $repeats_file >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/reduced_repeats.temp"
  repeats_file=$work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/reduced_repeats.temp"

  count_contig_rep=1
  walker_contig_rep=$(grep -c . $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_location.temp")

  repeat_storage=$work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/"$region_ID"_from_"$TargetID"_"$Syn_ID"_repeats.tab"
  repeats_summary=$work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/"$region_ID"_from_"$TargetID"_"$Syn_ID"_repeat_summary.tab"

  echo "TargetID SynID ContigID Start_Ragtag End_Ragtag Start_Contig End_contig Strand_contig Start_SynBlock End_SunBlock Start_Extracted End_Extracted" | tr " " "\t" >> $repeat_storage".controls"
  while [ $count_contig_rep -le $walker_contig_rep ]
  do
    # See this while loop? and the shenanigans saiving time on the repeat file.
    # This is what happens when you do the right thing and mask the repeats before the scaffolding and have to convert the coordinates
    # You are welcomed.
    contigID=$(sed -n $count_contig_rep"p" $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_location.temp" | awk -F "\t" '{print $1}')
    # contig_strand=$(sed -n $count_contig_rep"p" $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_location.temp" | awk -F "\t" '{print $6}')

    # Everything tells me this is wrong. That I need to account for the strand. Or if it is right that the unmasked scaffolding is Wrong
    # I cannot prove it. I've checked contigs in both scaffolds and match. And this recovers the desired area regardless.

    # I will finish the reporting and then work on a map comparison between masked and unmasked scaffolding.
    contig_start=$(sed -n $count_contig_rep"p" $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_location.temp" | awk -F "\t" -v start=$Syn_start '{if ($2>=start) print $4; else print $4+(start-$2)}')
    contig_end=$(sed -n $count_contig_rep"p" $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_location.temp" | awk -F "\t" -v end=$Syn_end '{if ($3<=end) print $5; else print $5-($3-end)}')
    contig_info=$(sed -n $count_contig_rep"p" $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_location.temp" )

    echo $TargetID" "$Syn_ID" "$contig_info" "$Syn_start" "$Syn_end | tr " " "\t" >> $repeat_storage".controls"

    awk -F "\t" -v chr=$contigID -v end=$contig_end '{if ($5==chr && $6 <= end) print}' $repeats_file | awk -F "\t" -v start=$contig_start '{if (start <= $7) print}' >> $repeat_storage

    count_contig_rep=$(($count_contig_rep + 1))
  done
  echo "Repeat Summary: "$TargetID" ("$Syn_ID")"

  summary_repeat
  Per_Repeats=$(echo $Bases_Repeats" "$Syn_lenght | awk -F " " '{printf "%.2f\n", $1/$2*100}'  | awk '{print $1"%"}')

  rm -f $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_location.temp"
  rm -f $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/contig_extract.temp"
  rm -f $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Repeats/reduced_repeats.temp"
}

genes_in_block_counting () {
  N_genes_in_block=$(grep -c . $check_found_genes_file)
  if [ $N_genes_in_block -gt 0 ]
  then
    N_Genes_on_target_zone=$(grep -w -F -c -f $check_found_genes_file $check_genes_file)
  else
    N_Genes_on_target_zone=0
  fi

  check_N_genes=$(grep -c . $check_genes_file)
  if [ $check_N_genes -gt 0 ]
  then
    # I had to fix a division by zero error and I'm too lazy to stream line the code, just for a few nanoseconds.
    Per_found=$(grep -c . $check_genes_file | awk -v found=$N_Genes_on_target_zone '{printf "%.2f\n", found/$1*100}'  | awk '{print $1"%"}')
  else
    Per_found=NA
  fi

  if [ $N_genes_in_block -gt 0 ]
  then
    Per_other=$(echo $N_Genes_on_target_zone" "$N_genes_in_block | awk -F " " '{printf "%.2f\n", ($2-$1)/$2*100}'  | awk '{print $1"%"}')
    grep -v -w -F -f $check_genes_file $check_found_genes_file >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id"_off_target.ID"

    missing_genes_file=$work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id"_off_target.ID"
    target_chromosome=$original_chr

    if [ $transfer_id == $original_genome_ID ]
    then
      annotation_file=$(awk -F "\t" -v genome=$transfer_id '{if (genome==$1) print $4}' $genome_table)
    else
      annotation_file=$work_folder"/"$annotation_transfers"/Individual_Transfers/"$transfer_id"_transfered_to_"$original_genome_ID"_annotation.gff_polished"
    fi

    process_genes_off_target

    if [ -f $missing_genes_file".info" ]
    then
      N_Transposintion_candidates=$(grep -w -F -c "Transposition_candidate" $missing_genes_file".info")
    else
      N_Transposintion_candidates=0
    fi
  else
    Per_other=NA
    N_Transposintion_candidates=0
  fi

  if [ $N_Genes_on_target_zone -gt 0 ]
  then
    Genes_on_target_zone=$(grep -w -F -f $check_found_genes_file $check_genes_file | tr "\n" ";" | sed "s/;$//")
  else
    Genes_on_target_zone=NA
  fi
}

generate_contig_gff () {
  scafolded_files=$(ls $scafolded_dir"/"*"_unmasked_scaffolds.fasta"  | sed "s/.*\///" | sed "s/_unmasked_scaffolds.fasta//")
  for pair in $scafolded_files
  do
    first_assembly=$(echo $pair | sed "s/_vs_/\n/" | head -n 1)
    second_assembly=$(echo $pair | sed "s/_vs_/\n/" | tail -n 1)

    echo $pair

    second_group=$(awk -F "\t" -v second=$second_assembly '{if ($1==second) print $2}' $genome_table)
    if [ $second_group == "Q" ]
    then
      ragtag_output=$work_folder"/"$query_group_ragtag"/Ref_"$first_assembly"/"$second_assembly"_ragtag_output/ragtag.scaffold.agp"
    else
      ragtag_output=$work_folder"/"$reference_group_ragtag"/Ref_"$first_assembly"_vs_Ref_"$second_assembly"_ragtag_output/ragtag.scaffold.agp"
    fi
    awk -F "\t" '{if ($5!="U" && $5!="N") print $1" Manual gene "$2" "$3" . "$9" . ID="$6}' $ragtag_output | grep -v "#" | tr " " "\t" >> $work_folder"/"$storage"/Additional_tracks/"$first_assembly"_vs_"$second_assembly"_contig.gff"
  done
}

confirm_final_region_file () {
  if [ ! -f $final_region_recovery ]
  then
    echo "No final region file found at: "
    echo $final_region_recovery
    exit
  else
    total_regions=$(grep -c . $final_region_recovery)
    total_white_spaces=$(grep -c " " $final_region_recovery)
    if [ $total_regions -lt 2 ]
    then
      echo "The final region file doesn't seems to have anything selected. Please verify and remember that the header is not counted"
      exit
    elif [ $total_white_spaces -gt 0 ]
    then
      echo "There are white spaces on the final region file. For safety none is allowed. They might affect the script's run"
      exit
    fi
  fi

  original_descriptor=$(echo $final_region_recovery | sed "s/.*_from_//" | sed "s/.in$//")
  check_ori=$(awk -F "\t" '{print $1}' $genome_table | grep -w -F -c $original_descriptor)
  if [ ! $check_ori -eq 1 ]
  then
    echo "Unclear original descriptor. Assembly "$check_ori" is not on the genone table provided:"
    echo $genome_table
    exit
  fi

  count=2
  while [ $count -le $total_regions ]
  do
    region_ID=$(sed -n $count"p" $final_region_recovery | awk -F "\t" '{print $1}')
    assembly_ID=$(sed -n $count"p" $final_region_recovery | awk -F "\t" '{print $2}')
    seq_file_path=$(sed -n $count"p" $final_region_recovery | awk -F "\t" '{print $3}')
    location_info=$(sed -n $count"p" $final_region_recovery | awk -F "\t" '{print $4}')
    is_scaffolded=$(sed -n $count"p" $final_region_recovery | awk -F "\t" '{print $5}')
    scaffold_reference=$(sed -n $count"p" $final_region_recovery | awk -F "\t" '{print $6}')

    if [ $scaffold_reference == "NA" ]
    then
      check_scaffold_ref=1
    else
      check_scaffold_ref=$(awk -F "\t" '{print $1}' $genome_table | grep -w -F -c $scaffold_reference)
    fi

    check_assembly=$(awk -F "\t" '{print $1}' $genome_table | grep -w -F -c $assembly_ID)
    check_location=$(echo $location_info | grep -c .)

    if [ ! -f $seq_file_path ]
    then
      echo "No file found at: "
      echo $seq_file_path
      exit
    elif [ ! $check_assembly -eq 1 ] || [ ! $check_location -eq 1 ]
    then
      echo "Something is wrong with the subregion information: "
      echo "Region ID: "$region_ID
      echo "Assembly ID: "$assembly_ID
      echo "Scaffold file path: "$seq_file_path
      echo "Subregion info: "$location_info
      exit
    elif [ ! $check_scaffold_ref -eq 1 ]
    then
      echo "The assembly "$scaffold_reference" is not a valid reference."
      echo "Check the ID used. And remember to indicate non-assembled "
      exit
    elif [ ! $is_scaffolded == "TRUE" ] && [ ! $is_scaffolded == "FALSE" ]
    then
      echo "Unclear indicator for scaffolded file: "$is_scaffolded
      echo "It can only be TRUE or FALSE"
      exit
    fi
    count=$(($count + 1))
  done
}

fraking_gene_things () {
  if [ ! -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"$assembly_ID"_info.txt" ]
  then
    echo "GeneID Gene_On_Target Chr Start End Loci_type Valid_ORFs Total_ORFs Identity Coverage Extra_Copies Description Original_Block Scaffolding_Reference" | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"$assembly_ID"_info.txt"
  fi

  echo $assembly_ID" "$donor" "$original" "$is_scaffolded" "$annot_file | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Region_work_summary.txt"

  # Necesary to ensure we get the description of the gene
  take_item=$(grep -w -F $donor $annotation_keystone | awk -F "\t" '{print $4}')

  block=1
  for target in $location_info
  do
    echo "Current Target: "$target
    chr=$(echo $target | awk -F ":" '{print $1}')
    start=$(echo $target | awk -F ":" '{print $2}' | sed "s/.$//" | awk -F "-" '{print $1}')
    end=$(echo $target | awk -F ":" '{print $2}' | sed "s/.$//" | awk -F "-" '{print $2}')
    strand=$(echo $target | awk -F ":" '{print $2}' | sed "s/./\n&/g" | tail -n 1 ) # if it is stupid but it works...

    ## Extract sequences for transcripts on every gene
    gffread -r $chr":"$start".."$end -w $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Sequences/temp_file.fasta" -g $seq_file_path $annot_file
    cat $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Sequences/temp_file.fasta" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Sequences/"$donor"_genes_on_"$assembly_ID".fasta"
    gffread -r $chr":"$start".."$end $annot_file | awk -F "\t" '{print $9}' >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Sequences/temp_file.annotation"

    echo "All Transcripts"
    all_transcript=$(seqkit seq -n -i $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Sequences/temp_file.fasta")
    echo "Done all transcripts"

    if [ ! -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Sequences/"$donor"_genes_on_"$assembly_ID"_IDs.txt" ]
    then
      echo "GenID TransID" | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Sequences/"$donor"_genes_on_"$assembly_ID"_IDs.txt"
    fi

    for transID in $all_transcript
    do
      gen_ID=$(grep -F "ID="$transID";" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Sequences/temp_file.annotation" | tr ";" "\n" | grep "geneID=" | sed "s/geneID=//")
      echo $gen_ID" "$transID | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Sequences/"$donor"_genes_on_"$assembly_ID"_IDs.txt"
    done

    rm -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Sequences/temp_file.fasta"
    rm -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Sequences/temp_file.annotation"

    if [ $strand == "+" ]
    then
      awk -F "\t" -v chr=$chr -v start=$start -v end=$end -v block=$block '{if ($1==chr && $3=="gene" && $4>=start && $5<=end) print $1"_ooo_"$4"_ooo_"$5"_ooo_"$7"_ooo_"$9"_ooo_Bl-"block}' $annot_file | sed "s/_ooo_/\t/g" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp"
    else
      awk -F "\t" -v chr=$chr -v start=$start -v end=$end -v block=$block '{if ($1==chr && $3=="gene" && $4>=start && $5<=end) print $1"_ooo_"$4"_ooo_"$5"_ooo_"$7"_ooo_"$9"_ooo_Bl-"block}' $annot_file | sed "s/_ooo_/\t/g" | sort -n -r -k 3 | awk -F "\t" -v block=$block '{if ($4=="+") print $1"_ooo_"$2"_ooo_"$3"_ooo_-_ooo_"$5"_ooo_Bl-"block; else print $1"_ooo_"$2"_ooo_"$3"_ooo_+_ooo_"$5"_ooo_Bl-"block}' | sed "s/_ooo_/\t/g" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp"
    fi

    awk -F "\t" -v chr=$chr -v start=$start -v end=$end '{if ($1==chr && $3=="mRNA" && $4>=start && $5<=end) print $9}' $annot_file | sort -u >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_product_annotation.temp"

    block=$(($block +1))
  done

  count_gene_info=1
  walker_gene_info=$(grep -c . $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp")

  while [ $count_gene_info -le $walker_gene_info ]
  do
    gene_ID=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" | awk -F "\t" '{print $5}' | tr ";" "\n" | grep ^"ID=" | sed "s/ID=//")

    if [ ! $original_descriptor == $donor ]
    then
      gene_on_target=$(cat $work_folder"/"$syn_block_inspection"/"$region_ID"/Original_Genome_Fetures/"*"_trans_"$donor"_genes.ID" | grep -w -F -c $gene_ID | awk '{if ($1==1) print "TRUE"; else print "FALSE"}')
    else
      gene_on_target=$(cat $work_folder"/"$syn_block_inspection"/"$region_ID"/Original_Genome_Fetures/"*"_original_genes.ID" | grep -w -F -c $gene_ID | awk '{if ($1==1) print "TRUE"; else print "FALSE"}')
    fi

    chr=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" |awk -F "\t" '{print $1}')
    coord_start=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" |awk -F "\t" '{print $2}')
    coord_end=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" |awk -F "\t" '{print $3}')
    type=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" | awk -F "\t" '{print $5}' | tr ";" "\n" | grep "^gene_biotype=" | sed "s/gene_biotype=//")
    target_block=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" |awk -F "\t" '{print $6}')
    check_description=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" | awk -F "\t" '{print $5}' | tr ";" "\n" | grep -c "^description=")

    if [ $original == FALSE ]
    then
      if [ $type == "protein_coding" ]
      then
        valid_orfs=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" | awk -F "\t" '{print $5}' | tr ";" "\n" | grep -c -w "valid_ORFs=0" | awk -F "\t" '{if ($1==0) print "TRUE"; else print "FALSE"}' )
        total_orfs=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" | awk -F "\t" '{print $5}' | tr ";" "\n" | grep ^"valid_ORFs=" | sed "s/valid_ORFs=//")
      else
        valid_orfs=NA
        total_orfs=NA
      fi

      identity=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" | awk -F "\t" '{print $5}' | tr ";" "\n" | grep ^"sequence_ID=" | sed "s/sequence_ID=//")
      coverage=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" | awk -F "\t" '{print $5}' | tr ";" "\n" | grep ^"coverage=" | sed "s/coverage=//")
      extra_copy_number=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" | awk -F "\t" '{print $5}' | tr ";" "\n" | grep ^"extra_copy_number=" | sed "s/extra_copy_number=//")
    else
      identity=Original
      coverage=Original
      extra_copy_number=Original
      valid_orfs=Original
      total_orfs=Original
    fi

    if [ $check_description -gt 0 ]
    then
      description=$(sed -n $count_gene_info"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp" | awk -F "\t" '{print $5}' | tr ";" "\n" | grep "^description=" | sed "s/description=//")
    else
      if [ ! $extra_copy_number == 0 ] && [ ! $extra_copy_number == Original ]
      then
        search_genid=$(echo $gene_ID | sed "s/^gene-//")
        remove_copy_num=$(echo $search_genid | sed "s/./\n&/g" | grep . | grep -n "_" | tail -n 1 | awk -F ":" '{print $1-1}')
        search_genid=$(echo $search_genid | sed "s/./\n&/g" | grep . | head -n $remove_copy_num | tr -d "\n" | sed "s/$/\n/")
        search_ID_for_product=$(echo $take_item$search_genid";")
      else
        search_ID_for_product=$(echo $take_item"__Delete__"$gene_ID";" | sed "s/__Delete__gene-//")
      fi

      check_product=$(grep $search_ID_for_product $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_product_annotation.temp" | tr ";" "\n" | grep -c "^product=")
      if [ $check_product -gt 0 ]
      then
        description=$(grep $search_ID_for_product $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_product_annotation.temp" | tr ";" "\n" | grep "^product="  | sed "s/product=//" | sort -u | tr "\n" ";" | sed "s/;$/\n/" | sed "s/;/__-__/g")
      else
        description=NA
      fi
    fi

    echo $gene_ID";"$gene_on_target";"$chr";"$coord_start";"$coord_end";"$type";"$valid_orfs";"$total_orfs";"$identity";"$coverage";"$extra_copy_number";"$description";"$target_block";"$scaffold_reference | tr ";" "\t" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"$assembly_ID"_info.txt"
    count_gene_info=$(($count_gene_info + 1))
  done
  rm -rf $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_info.temp"
  rm -rf $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/"$donor"_of_"$assembly_ID"_product_annotation.temp"
}

get_special_contig_genes () {
  mkdir $work_folder"/"$storage"/Additional_tracks/Special_"$contig

  for RefID in $reference_assemblies
  do
    echo "Getting target genes on "$contig": "$RefID
    if [ $RefID == $source ]
    then
      gff_file=$(awk -F "\t" -v genome=$RefID '{if (genome==$1) print $4}' $genome_table)
      awk -F "\t" -v wanted=$contig '{if (wanted==$1 && $3=="gene") print}' $gff_file >> $work_folder"/"$storage"/Additional_tracks/Special_"$contig"/"$RefID"_original.gff"
      awk -F "\t" '{print $9}' $work_folder"/"$storage"/Additional_tracks/Special_"$contig"/"$RefID"_original.gff" | tr ";" "\n" | grep "ID=" | tr -d " " | sed "s/$/;/" >> $work_folder"/"$storage"/Additional_tracks/Special_"$contig"/"$RefID"_search.id"
    else
      gff_file=$work_folder"/"$annotation_transfers"/Individual_Transfers/"$RefID"_transfered_to_"$source"_annotation.gff_polished"
      awk -F "\t" -v wanted=$contig '{if (wanted==$1 && $3=="gene") print}' $gff_file >> $work_folder"/"$storage"/Additional_tracks/Special_"$contig"/"$RefID"_original.gff"
      awk -F "\t" '{print $9}' $work_folder"/"$storage"/Additional_tracks/Special_"$contig"/"$RefID"_original.gff" | tr ";" "\n"  | tr -d " " | grep ^ID= | sed "s/$/;/" >> $work_folder"/"$storage"/Additional_tracks/Special_"$contig"/"$RefID"_search.id"
    fi
  done

  for RefID in $reference_assemblies
  do
    all_transfer_files=$(ls $work_folder"/"$annotation_transfers"/Individual_Transfers/"$RefID"_transfered_to_"*"_annotation.gff_polished" | sed "s/\/$//" | sed "s/.*\///" | grep -v "_transfered_to_"$source"_annotation.gff_polished" | sed "s/_annotation.gff_polished//")
    for transfer in $all_transfer_files
    do
      echo $RefID" "$contig": "$transfer
      awk -F "\t" '{if ($3=="gene") print}' $work_folder"/"$annotation_transfers"/Individual_Transfers/"$transfer"_annotation.gff_polished" | grep -F -f $work_folder"/"$storage"/Additional_tracks/Special_"$contig"/"$RefID"_search.id" >> $work_folder"/"$storage"/Additional_tracks/Special_"$contig"/"$transfer"_custom_track.gff"
    done
  done
}

do_genes_on_target () {
  for transfer in $transfer_gff_files
  do
    transfer_id=$(echo $transfer | sed "s/.*\///" | sed "s/_transfered_to_/\n/" | head -n 1 )
    if [ $transfer_id == $original_genome_ID ]
    then
      check_genes_file=$work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_original_genes.ID"
    else
      check_genes_file=$work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_trans_"$transfer_id"_genes.ID"
    fi

    if [ ! -f $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID"/Target_genes_from_"$transfer_id".gff" ]
    then
      grep -w -F -f $check_genes_file $transfer  | awk -F "\t" '{if ($3=="gene") print}' >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID"/Target_genes_from_"$transfer_id".gff"
    fi
  done
}

detailed_gene_track () {
  echo "Getting genes... "$print_current
  if [ ! -f $out_track ]
  then
    echo "##gff-version 3" > $out_track
  fi

  if [ $run_duplicated == TRUE ]
  then
    count_dup=1
    all_dup=$(grep -c . $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/duplicated.tmp")

    while [ $count_dup -le $all_dup ]
    do
      dup_id=$(sed -n $count_dup"p" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/duplicated.tmp" | awk -F "\t" '{print $1}')
      search_id=$(sed -n $count_dup"p"  $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/duplicated.tmp" | awk -F "\t" '{print $2}')

      check_dup_gene=$(grep -c $search_id $wanted_genes)
      if [ $check_dup_gene -gt 0 ]
      then
        echo $take_item$search_id >> $wanted_genes
      fi

      count_dup=$(($count_dup + 1))
    done
  fi

  grep -F -f $wanted_genes $current_annotation_file >> $out_track
  check_abemus_data=$(grep -c . $out_track)
  if [ $check_abemus_data -eq 0 ]
  then
    rm -f $out_track
  fi
  rm -f $wanted_genes
}

################################################################################
################################################################################
################################################################################

# Block of code
# Now is the script proper doing its thing

if [ $VALIDATION_AND_STATS == TRUE ]
then
  # This block of codes verifies that the genome assemblies are properly indicated
  # and basic folders exist. Be warned: The validation process is not exahustive.

  # 1) Check if working folder exists
  if [ ! -d $work_folder ]
  then
    # The work folder does not exits...
    echo "Working folder does not exist. Check this: "
    echo "work_folder="$work_folder
    exit
  else
    echo "Work Folder is real, good..."
  fi

  if [ ! -d $cluster_submission_folder ]
  then
    # The work folder does not exits...
    echo "The folder for the cluster submissions not exist. Check this: "
    echo "cluster_submission_folder="$cluster_submission_folder
    exit
  else
    echo "...There is a place where to store batch jobs..."
  fi

  # 2) Check if Genome files exists and has the right overall data
  confirm_input_file

  echo "Basic checks completed, checking if every assembly is actually real."
  echo "Stand by while seqkit gets stats for each one"

  if [ -d $work_folder"/Assembly_Stats" ]
  then
    rm -rf $work_folder"/Assembly_Stats"
  fi

  mkdir $work_folder"/Assembly_Stats"

  echo "ID Group File Format Type Num_seqs Sum_len Min_len Avg_len Max_len Q1 Q2 Q3 Sum_gap N50 N50_num Q20(%) Q30(%) AvgQual GC(%)" | tr " " "\t" >> $work_folder"/Assembly_Stats/Assemblies_stats.txt"
  echo "ID contig/scaffold Lenght GC% GC_Skew" | tr " " "\t" > $work_folder"/Assembly_Stats/Sequence_stats.txt"

  # Checking each genome assembly and GFF are fine
  count=2
  walker=$(grep -c . $genome_table)
  while [ $count -le $walker ]
  do
    path_to_genome=$(sed -n $count"p" $genome_table | awk -F "\t" '{print $3}')
    path_to_prot=$(sed -n $count"p" $genome_table | awk -F "\t" '{print $5}')

    check_good_data=$(echo $path_to_genome | grep -c .)
    if [ $check_good_data -eq 0 ]
    then
      echo "Cannot retrieve location of for genome at line "$count
      echo "Abort"
      exit
    elif [ ! -f $path_to_genome ]
    then
      echo "Cannot retrieve location of for genome at line "$count
      echo "Abort"
      exit
    else
      assembly_ID=$(sed -n $count"p" $genome_table | awk -F "\t" '{print $1}')
      group_ID=$(sed -n $count"p" $genome_table | awk -F "\t" '{print $2}')

      echo "Checking "$assembly_ID" - "$group_ID" - "$path_to_genome

      seqkit stats -T -a $path_to_genome > $work_folder"/Assembly_Stats/temp_stats.temp"
      check_stats=$(grep -c . $work_folder"/Assembly_Stats/temp_stats.temp")

      if [ $check_stats -eq 2 ]
      then
        data=$(tail -n 1 $work_folder"/Assembly_Stats/temp_stats.temp")
        echo $assembly_ID" "$group_ID" "$data | tr " " "\t" >> $work_folder"/Assembly_Stats/Assemblies_stats.txt"
        seqkit fx2tab -n -l -g -G $path_to_genome | awk -F "\t" -v name=$assembly_ID '{print name"\t"$0}' >> $work_folder"/Assembly_Stats/Sequence_stats.txt"
        rm -f $work_folder"/Assembly_Stats/temp_stats.temp"
      else
        rm -f $work_folder"/Assembly_Stats/temp_stats.temp"
        echo "Something is bad with this file. Abort"
        exit
      fi
    fi

  # Checking GFF file.
    gff_file=$(sed -n $count"p" $genome_table | awk -F "\t" '{print $4}')
    if [ ! $gff_file == "NA" ]
    then
      if [ ! -f $gff_file ]
      then
        echo "Can't find annotation the file: " $gff_file
        exit
      else
        check_first_line=$(grep -v "#" $gff_file | head -n 1 | tr "\t" "\n" | grep -c .)
        check_last_line=$(grep -v "#" $gff_file | tail -n 1 | tr "\t" "\n" | grep -c .)

        if [ $check_first_line == $check_last_line ] && [ $check_first_line -eq 9 ]
        then
          echo "Annotation file seems good. Moving on"
        else
          echo "Something is wrong with. Check its integrity with GFFread: " $gff_file
          exit
        fi
      fi
    fi

    # Checking protein file
    check_good_data=$(echo $path_to_prot | grep -c .)
    if [ $check_good_data -eq 0 ]
    then
      echo "Cannot retrieve location of for proteins at line "$count
      echo "Abort"
      exit
    elif [ ! -f $path_to_genome ] && [ $path_to_genome != "NA" ]
    then
      echo "Cannot retrieve location of for proteins at line "$count
      echo "Abort"
      exit
    else
      if [ $path_to_prot != "NA" ]
      then
        seqkit stats -T -a $path_to_prot > $work_folder"/Assembly_Stats/temp_stats.temp"
        check_stats=$(grep -c . $work_folder"/Assembly_Stats/temp_stats.temp")

        if [ $check_stats -eq 2 ]
        then
          rm -f $work_folder"/Assembly_Stats/temp_stats.temp"
        else
          rm -f $work_folder"/Assembly_Stats/temp_stats.temp"
          echo "Something is bad with file: "$path_to_prot
          exit
        fi
      fi
    fi

    count=$(($count + 1))
  done
  echo ""
  echo "All seems in order..."
  echo "Just REMEMBER!!! The script assumes there is a header in the file and ignores the first line."
  echo "I will change that beheavior when I have confirmation humans have telepathy. Until then keep a header explaining what everything is."

  echo "Disclaimer 1: to validate Fasta files I use seqkit, but recently I found out that empty files are considered valid for some reason."
  echo "Please don't use empty files."
  echo ""
  echo "Disclaimer 2: the GFF files are not trully validated. GFFread can do it but it takes too long"
  echo "This script only checks that the file exists and that that the first and last entry have 9 columns"
fi

################################################################################

if [ $VALIDATION_REGION == TRUE ]
then
  echo "Checking desired Region file..."
  confirm_region_file
  echo "File seems in order... just remember the script won't check if you put letters where coordinates or lenght should be"
  echo ""
  echo "Checking regions..."

  N_region_entries_with_primers=$(awk -F "\t" '{if ($9 != "NA") print $9}' $region_table | sed 1d | grep -c .)

  if [ $N_region_entries_with_primers -eq 0 ]
  then
    echo "No region seems to have a primer file associated, guess that's all."
  else
    echo $N_region_entries_with_primers" have a primer file asigned. Time to check it"

    all_regions=$(awk -F "\t" '{print $1}' $region_table | sed 1d)
    for regionID in $all_regions
    do
      genome_ID=$(awk -F "\t" -v region=$regionID '{if ($1==region) print $2}' $region_table)
      chr_ID=$(awk -F "\t" -v region=$regionID '{if ($1==region) print $3}' $region_table)

      check_genome_in_analysis=$(awk -F "\t" -v genome=$genome_ID '{if ($1==genome) print}' $genome_table | grep -c .)
      check_chr_in_analysis=$(awk -F "\t" -v genome=$genome_ID '{if ($1==genome) print}' $work_folder"/Assembly_Stats/Sequence_stats.txt"| grep -w -F -c $chr_ID)

      if [ $check_genome_in_analysis -eq 0 ]
      then
        echo $regionID" is not on the genome table! Aborting"
        echo "check your speLing on: "$genome_table
        exit
      elif [ $check_chr_in_analysis -eq 0 ]
      then
        echo $chr_ID" is not part of the assembly "$genome_ID
      fi
      echo "Check Region: "$regionID
    done

    regions_with_primers=$(awk -F "\t" '{if ($9 != "NA") print $1}' $region_table | sed 1d)
    for regionID in $regions_with_primers
    do
      primer_file_location=$(awk -F " " -v region=$regionID '{if ($1==region) print $9}' $region_table)
      if [ ! -f $primer_file_location ]
      then
        echo "Primer File for "$regionID" is missing."
        echo "Check: "$primer_file_location
        exit
      else
        check_file_is_fasta=$(seqkit stats $primer_file_location | grep -c .)
        if [ ! $check_file_is_fasta -eq 2 ]
        then
          echo "Something is wrong with file: "$primer_file_location
          echo "doesn't seems to be a fasta"
          exit
        else
          check_anchor_primers
        fi
      fi
    done
  fi

  echo ""
  echo "All seems in order..."
  echo "Just REMEMBER!!! The script assumes there is a header in the file and ignores the first line."
  echo "I will change that beheavior when I have confirmation humans have telepathy. Until then keep a header explaining what everything is."
fi

################################################################################

if [ $SLURM_REPEATMODELER == TRUE ]
then
  storage=REPEAT_MODELER_RUN

  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Repeats_"$current_date)

  confirm_input_file
  start_block_code

  rm -f $work_folder"/"$storage"/Already_done.txt"

  best_assembly_refID=$(awk -F "\t" '{if ($2=="R") print $1"\t"$15}' $work_folder"/Assembly_Stats/Assemblies_stats.txt" | sort -r -n -k2 | awk -F "\t" '{print $1}' | head -n 1)
  genome_file=$(awk -F "\t" -v genome=$best_assembly_refID '{if (genome==$1) print $3}' $genome_table)

  echo "Selected assembly: "$best_assembly_refID
  echo $genome_file
  echo ""
  echo "BuildDatabase "$genome_file" -name "$work_folder"/"$storage"/"$best_assembly_refID"_reference" >>  $cluster_submission_folder"/1_"$submit_to_cluster"_BuildDatabase.txt"
  echo "RepeatModeler -threads "'$NPROCS'" -database "$work_folder"/"$storage"/"$best_assembly_refID"_reference" >> $cluster_submission_folder"/2_"$submit_to_cluster"_RepeatModeler.txt"

  code_block_completed
  slurm_block_completed

  echo "Special Note: "
  echo "Repeat Modeler is... annoying to say the least. Many intermediary files dumped on the working folder. So wherever you launched the upload to SLURM"
  echo "Go back where you launched the run, search for a file named: RM_[Number_ID].[Date] and compreses it or delete it"
  echo ""
fi

################################################################################

if [ $SLURM_REPEATMASKER_Re == TRUE ]
then
  storage=REPEAT_MASKER_RUN
  repeat_modeler_storage=REPEAT_MODELER_RUN

  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Repeats_"$current_date)

  confirm_input_file
  start_block_code

  if [ ! $override_repeatmodeler == "FALSE" ]
  then
    repeat_modeler_ref=$override_repeatmodeler
  else
    repeat_modeler_ref=$(ls $work_folder"/"$repeat_modeler_storage"/"*"_reference-families.fa")
  fi

  check_references=$(echo $repeat_modeler_ref | tr " " "\n" | grep -c .)
  if [ ! $check_references -eq 1 ]
  then
    echo "Found "$check_references" potential files for repeatmasker. There should be one alone."
    echo $repeat_modeler_ref
    exit
  fi

  genomes_IDs=$(awk -F "\t" '{print $1}' $genome_table | sed 1d | sort -u)
  for genome_id in $genomes_IDs
  do
    check_done=$(grep -w -F -c $genome_id $work_folder"/"$storage"/Already_done.txt")

    if [ $check_done -gt 0 ]
    then
      echo $genome_id" is already done. Skipping!"
    else
      if [ -d $work_folder"/"$storage"/"$genome_id ]
      then
        rm -rf $work_folder"/"$storage"/"$genome_id
      fi

      genome_file=$(awk -F "\t" -v genome=$genome_id '{if (genome==$1) print $3}' $genome_table)
      echo "RepeatMasker -e rmblast -pa 2 -lib "$repeat_modeler_ref" "$genome_file" -dir "$work_folder"/"$storage"/"$genome_id >> $cluster_submission_folder"/"$submit_to_cluster"_RepeatMasker.txt"
      echo $genome_id >> $work_folder"/"$storage"/Already_done.txt"
    fi
  done
  code_block_completed
  slurm_block_completed

  echo ""
  echo "Special note: Remember that RepeatMasker has special rules for  its procesors. Here is the explanation from the -h:"
  echo "----------"
  echo "-pa(rallel) [number]"
  echo "The number of sequence batch jobs [50kb minimum] to run in parallel. RepeatMasker will fork off this number of parallel jobs, each running the search engine specified. For each search engine invocation ( where applicable ) a fixed the number of cores/threads is used:"
  echo "RMBlast: 4 | ABBlast: 4 | nhmmer: 2 | crossmatch: 1"
  echo "To estimate the number of cores a RepeatMasker run will use simply multiply the -pa value by the number of cores the particular search engine will use."
  echo "----------"
  echo "This script sets -pa 2, therefore you need to ask at least 8 procesors."

fi

################################################################################

# Repeatmasker is annoying and evil.its output table is formated with arbitrary numbers of white spaces and I don't like that
# This if will make a civilized copy.

if [ $REFORMAT_REPEATMASKER == TRUE ]
then
  storage=REPEAT_MASKER_RUN
  assemblies=$(ls -d $work_folder"/"$storage"/"* | grep -v .txt$ | sed "s/\/$//" | sed "s/.*\///")

  for assembly in $assemblies
  do
    echo $assembly
    repeats_file=$(ls  $work_folder"/"$storage"/"$assembly"/"*".out")
    if [ -f $repeats_file".tab" ]
    then
      rm -f $repeats_file".tab"
    fi
    sed 1,2d $repeats_file | tr -s " " | sed "s/^ //" | tr " " "\t" >> $repeats_file".tab"
  done
fi

################################################################################

if [ $SLURM_RAGTAG_SCAFFOLDING_Re == TRUE ]
then
  # This block of code aligns all posible pairs between query assemblies (Q) against reference ones (R)
  # using minimap2. Then it compresses the SAM file and generates an index for further use.

  storage=RAGTAG_Scaffolding
  masked_genomes=REPEAT_MASKER_RUN

  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Scaffolding_"$current_date)

  confirm_input_file
  start_block_code

  query_assemblies=$(awk -F "\t" '{if ($2=="Q") print $1}' $genome_table | sort)
  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $genome_table | sort)

  if [ ! -d $work_folder"/"$storage"/Temp_file_sequences" ]
  then
    mkdir $work_folder"/"$storage"/Temp_file_sequences"
  fi

  for refID in $reference_assemblies
  do
    echo "Working relative to Reference: "$refID

    if [ ! -d $work_folder"/"$storage"/Ref_"$refID ]
    then
      mkdir $work_folder"/"$storage"/Ref_"$refID
    fi

    if [ $override_repeatmasker == "TRUE" ]
    then
      ref_genome=$(awk -F "\t" -v wanted=$refID '{if ($1==wanted) print $3}' $genome_table )
    else
      ref_genome=$(ls $work_folder"/"$masked_genomes"/"$refID"/"*".masked")
    fi

    index_genome=$ref_genome
    indexID=$refID
    index_fasta

    for queryID in $query_assemblies
    do
      echo "Scaffolding "$refID" vs "$queryID
      check_done=$(grep -w -c "This_"$refID"_vs_"$queryID"_Done" $work_folder"/"$storage"/Already_done.txt")

      if [ $check_done -eq 0 ]
      then
        remove_ragtag_previous_run

        if [ $override_repeatmasker == "TRUE" ]
        then
          target_genome=$(awk -F "\t" -v wanted=$queryID '{if ($1==wanted) print $3}' $genome_table )
        else
          target_genome=$(ls $work_folder"/"$masked_genomes"/"$queryID"/"*".masked")
        fi

        index_genome=$target_genome
        indexID=$queryID
        index_fasta

        echo "ragtag.py scaffold -u -o "$work_folder"/"$storage"/Ref_"$refID"/"$queryID"_ragtag_output -t "'$NPROCS'" "$work_folder"/"$storage"/Temp_file_sequences/"$refID"_temp_copy.fasta "$work_folder"/"$storage"/Temp_file_sequences/"$queryID"_temp_copy.fasta" >> $cluster_submission_folder"/"$submit_to_cluster"_ragtag.txt"
        echo "This_"$refID"_vs_"$queryID"_Done" >> $work_folder"/"$storage"/Already_done.txt"
      fi
    done
  done
  code_block_completed
  slurm_block_completed

  echo ""
  echo "---------------------------------------------------------------"
  echo ""
  echo "Special Note 1: After you are done with the SLURM go back to:"
  echo $work_folder"/"$storage
  echo "And delete the folder Temp_file_sequences, or just run SLURM_RAGTAG_SUMARY to remove it."
  echo "RagTag needed access to a fasta.fai file and nor every genome had it."
  echo "So for ease of use the script copies the genome to this folder and creates one"
  echo "For the most part this is just redundant information and re-creating doesn't takes long."
  echo ""
  echo "---------------------------------------------------------------"
  echo ""
  echo "Special Note 2:"
  echo "Since I was copying the genomes, I took the liberty of adding the sample ID at the start of"
  echo "every sequence. This might make the ragtag.scaffold.agp file easier to read,"
  echo "but i'm going to reverse it on the next block of code for file ragtag.scaffold.fasta"
fi

################################################################################

if [ $SLURM_RAGTAG_SCAFFOLDING_REFERENCES_Re == TRUE ]
then
  # This block of code aligns all posible pairs between query assemblies (Q) against reference ones (R)
  # using minimap2. Then it compresses the SAM file and generates an index for further use.

  storage=RAGTAG_Scaffolding_References
  masked_genomes=REPEAT_MASKER_RUN
  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Scaffolding_References_"$current_date)

  confirm_input_file
  start_block_code

  if [ ! -f $work_folder"/Assembly_Stats/Assemblies_stats.txt" ]
  then
    echo "Code block SLURM_RAGTAG_SCAFFOLDING_REFERENCES_Re requires the file Assemblies_stats.txt to get and order the references."
    echo "Please run VALIDATION_AND_STATS first"

    echo "The file should be here: "$work_folder"/Assembly_Stats/Assemblies_stats.txt"
    exit
  fi

  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1"\t"$15}' $work_folder"/Assembly_Stats/Assemblies_stats.txt" | sort -r -n -k2 | awk -F "\t" '{print $1}')
  N_references=$(echo $reference_assemblies | tr " " "\n" | grep -c .)

  if [ $N_references -le 1 ]
  then
    echo "Not enough references to compare. Only "$N_references" in work table"
    echo $reference_assemblies
  else
    if [ ! -d $work_folder"/"$storage"/Temp_file_sequences" ]
    then
      mkdir $work_folder"/"$storage"/Temp_file_sequences"
    fi

    for refID1 in $reference_assemblies
    do
      echo "Working relative to Reference: "$refID1

      if [ $override_repeatmasker == "TRUE" ]
      then
        ref_genome1=$(awk -F "\t" -v wanted=$refID1 '{if ($1==wanted) print $3}' $genome_table )
      else
        ref_genome1=$(ls $work_folder"/"$masked_genomes"/"$refID1"/"*".masked")
      fi

      index_genome=$ref_genome1
      indexID=$refID1
      index_fasta

      for refID2 in $reference_assemblies
      do
        if [ ! $refID1 == $refID2 ]
        then
          echo "Scaffolding "$refID1" vs "$refID2
          check_done=$(grep "_"$refID1"_" $work_folder"/"$storage"/Already_done.txt" | grep -c "_"$refID2"_")
          if [ $check_done -eq 0 ]
          then
            remove_ragtag_previous_run
            ref_genome2=$(awk -F "\t" -v wanted=$refID2 '{if ($1==wanted) print $3}' $genome_table )

            if [ $override_repeatmasker == "TRUE" ]
            then
              ref_genome2=$(awk -F "\t" -v wanted=$refID2 '{if ($1==wanted) print $3}' $genome_table )
            else
              ref_genome2=$(ls $work_folder"/"$masked_genomes"/"$refID2"/"*".masked")
            fi

            index_genome=$ref_genome2
            indexID=$refID2
            index_fasta

            echo "ragtag.py scaffold -u -o "$work_folder"/"$storage"/Ref_"$refID1"_vs_Ref_"$refID2"_ragtag_output -t "'$NPROCS'" "$work_folder"/"$storage"/Temp_file_sequences/"$refID1"_temp_copy.fasta "$work_folder"/"$storage"/Temp_file_sequences/"$refID2"_temp_copy.fasta" >> $cluster_submission_folder"/"$submit_to_cluster"_ragtag.txt"
            echo "This_"$refID1"_vs_"$refID2"_Done" >> $work_folder"/"$storage"/Already_done.txt"
          fi
        fi
      done
    done

    code_block_completed
    slurm_block_completed

    echo ""
    echo "---------------------------------------------------------------"
    echo ""
    echo "Special Note 1: After you are done with the SLURM go back to:"
    echo $work_folder"/"$storage
    echo "And delete the folder Temp_file_sequences, or just run SLURM_RAGTAG_SUMARY to remove it."
    echo "RagTag needed access to a fasta.fai file and nor every genome had it."
    echo "So for ease of use the script copies the genome to this folder and creates one"
    echo "For the most part this is just redundant information and re-creating doesn't takes long."
    echo ""
    echo "---------------------------------------------------------------"
    echo ""
    echo "Special Note 2:"
    echo "Since I was copying the genomes, I took the liberty of adding the sample ID at the start of"
    echo "every sequence. This might make the ragtag.scaffold.agp file easier to read,"
    echo "but i'm going to reverse it on the next block of code for file ragtag.scaffold.fasta"
  fi
fi

################################################################################

if [ $SLURM_RAGTAG_SUMARY == TRUE ]
then
  # Thise code block summarizes all the RagTag results and checks if everything is in order to continue
  # working on this data. Or if some scaffolding attempts failed / are missing.

  # Because my coding skills have limits, if for whatever reason tomorrow there are more groups of samples
  # than B and R, this block will need to be updated. I will try it to make easier by setting most of it into
  # a single variabe.

  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Scaffold_cords_"$current_date)

  confirm_input_file

  storage=RAGTAG_Scaffolding
  query_assemblies=$(awk -F "\t" '{if ($2=="Q") print $1}' $genome_table | sort)
  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $genome_table | sort)

  start_block_code

  comparison_group=Q_vs_R # Assemblies from the lab vs a Reference

  ragtag_summary_start

  for refID in $reference_assemblies
  do
    echo "Working relative to Reference: "$refID
    ragtag_output_path=$work_folder"/"$storage"/Ref_"$refID

    for queryID in $query_assemblies
    do
      missing=$(echo .)
      error=$(echo .)

      ragtag_output_folder=$queryID"_ragtag_output"
      target_genome=$(awk -F "\t" -v assembly_id=$queryID '{if ($1==assembly_id) print $3}' $genome_table)
      print_id=$(echo $ragtag_output_folder | sed "s/_ragtag_output$//" | awk -v add=$refID '{print add"_vs_"$0}')
      recover_information_ragtag_runs

    done
  done

  registry=$ragtag_registry_missing
  source=$(echo "Missing Files")

  ragtag_registry_missing=NA
  ragtag_registry_error=NA
  print_issues_ragtag

  registry=$ragtag_registry_error
  source=$(echo "RagTag Error")
  print_issues_ragtag

  ##############################################
  ##############################################
  ##############################################

  storage=RAGTAG_Scaffolding_References
  comparison_group=R_vs_R # Assemblies from the lab vs a Reference

  ragtag_summary_start

  ragtag_output_path=$work_folder"/"$storage
  for refID1 in $reference_assemblies
  do
    for refID2 in $reference_assemblies
    do
      if [ ! $refID1 == $refID2 ]
      then
        missing=$(echo .)
        error=$(echo .)
        ragtag_output_folder="Ref_"$refID1"_vs_Ref_"$refID2"_ragtag_output"


        if [ -d $ragtag_output_path"/"$ragtag_output_folder ]
        then
          # This might be an issue eventually: the script doesn't checks if a Ref vs Ref scaffolding is missing or not
          print_id=$(echo $refID1"_vs_"$refID2)
          recover_information_ragtag_runs
        fi
      fi
    done
  done

  echo ""
  echo "RAGTAG Summary complete"
  echo "Note: in its current form, the script doesn't checks if all the reference assemblies comparisons are present or not"
fi

################################################################################

if [ $SLURM_RAGTAG_MINIMAP_SYRI_REGION_Re == TRUE ]
then
  # This block of code will identify rougly where everything is relative to each other
  storage=RagTag_Plus_Syri_Region
  query_group_ragtag=RAGTAG_Scaffolding
  reference_group_ragtag=RAGTAG_Scaffolding_References
  masked_genomes=REPEAT_MASKER_RUN

  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Region_Focused_Syri_"$current_date)

  confirm_input_file
  confirm_region_file
  start_block_code

  count_regions=2
  recorrer_regions=$(grep -c . $region_table)

  while [ $count_regions -le $recorrer_regions ]
  do
    region_ID=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $1}')

    rename_fasta_seq=$(echo ">"$region_ID"_sequence")

    original_genome_ID=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $2}')
    original_chr_ID=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $3}')
    original_group=$(awk -F "\t" -v genome=$original_genome_ID '{if (genome==$1) print $2}' $genome_table)
    original_genome_file=$(awk -F "\t" -v genome=$original_genome_ID '{if (genome==$1) print $3}' $genome_table)

    check_region_done=$(grep -w -F -c $original_chr_ID $work_folder"/"$storage"/Already_done.txt")

    echo "Doing: "$region_ID" || "$original_chr_ID
    echo $count_regions" / "$recorrer_regions

    if [ $check_region_done -gt 0 ]
    then
      echo $original_chr_ID" Already done. Skipping."
    else
      echo $original_genome_ID" "$original_group" "$original_chr_ID
      echo $original_genome_file
      echo ""

      local_work_folder=$work_folder"/"$storage"/"$original_chr_ID
      if [ -d $local_work_folder ]
      then
        rm -rf $local_work_folder
      fi

      mkdir $local_work_folder
      mkdir $local_work_folder"/Minimap"
      mkdir $local_work_folder"/Syri"

      if [ ! -f $local_work_folder"/Run_Registrty.txt" ]
      then
        echo "Genome_v_Genome Reference_Fasta Target_Fasta" | tr " " "\t" >> $local_work_folder"/Run_Registrty.txt"
      fi

      reference_fasta_minimap=$local_work_folder"/"$original_genome_ID"_reference.fasta"
      seqkit grep -j $threads -p $original_chr_ID $original_genome_file | seqkit seq -u | sed "s/>.*/$rename_fasta_seq/" > $reference_fasta_minimap

      target_genome_IDs=$(sed 1d $genome_table | awk -F "\t" '{print $1}' | grep -v -F -w $original_genome_ID | sort -u)
      for targetID in $target_genome_IDs
      do
        echo "Current Target ID: "$targetID

        target_group=$(awk -F "\t" -v genome=$targetID '{if (genome==$1) print $2}' $genome_table)
        # target_genome_file=$(awk -F "\t" -v genome=$targetID '{if (genome==$1) print $3}' $genome_table)

        store_result_ID=$(echo $original_genome_ID"_vs_"$targetID"_"$original_chr_ID)

        if [ $original_group == "R" ] && [ $target_group == "Q" ]
        then
          reference_fasta_minimap=$local_work_folder"/"$original_genome_ID"_reference.fasta" # This is awkward. But I need to update the file name to avoid a bug every time there is a R vs R comparison
          ragtag_output_folder=$work_folder"/"$query_group_ragtag"/Ref_"$original_genome_ID"/"$targetID"_ragtag_output"
          matching_sequence=$(echo $original_genome_ID"_"$original_chr_ID"_RagTag")
          target_fasta_minimap=$local_work_folder"/"$store_result_ID"_target.fasta"

          seqkit grep -j $threads -p $matching_sequence $ragtag_output_folder"/ragtag.scaffold.fasta" | sed "s/>.*/$rename_fasta_seq/" > $target_fasta_minimap
          do_ragtag_syri
        elif [ $original_group == "R" ] && [ $target_group == "R" ]
        then
          confirm_reference=$(ls -d $work_folder"/"$reference_group_ragtag"/"* | grep -v .txt$ | sed "s/\/$//" | sed "s/.*\///" | sed "s/ragtag_output//" | grep "_"$original_genome_ID"_" | grep "_"$targetID"_" | sed "s/_vs_/_\n_/" | grep -n "_"$original_genome_ID"_" | awk -F ":" '{print $1}')

          if [ $confirm_reference -eq 1 ]
          then
            ragtag_output_folder=$work_folder"/"$reference_group_ragtag"/Ref_"$original_genome_ID"_vs_Ref_"$targetID"_ragtag_output"
            reference_fasta_minimap=$local_work_folder"/"$original_genome_ID"_reference.fasta"

            matching_sequence=$(echo $original_genome_ID"_"$original_chr_ID"_RagTag")

            target_fasta_minimap=$local_work_folder"/"$store_result_ID"_target.fasta"
            seqkit grep -j $threads -p $matching_sequence $ragtag_output_folder"/ragtag.scaffold.fasta" | sed "s/>.*/$rename_fasta_seq/" > $target_fasta_minimap
            do_ragtag_syri
          else
            # Just for this section, some of the roles are reversed. The alternative would be to run RAGTAG for all possible combinations of the reference assemblies.
            # Not just using use the best assembly as a reference and the other as query
            ragtag_output_folder=$work_folder"/"$reference_group_ragtag"/Ref_"$targetID"_vs_Ref_"$original_genome_ID"_ragtag_output"
            # target_genome_file=$(awk -F "\t" -v genome=$original_genome_ID '{if (genome==$1) print $3}' $genome_table)

            matching_sequence=$(grep -w -F $original_genome_ID"_"$original_chr_ID $ragtag_output_folder"/ragtag.scaffold.agp" | awk -F "\t" '{print $1}' | sort -u)
            confirm_match=$(echo $matching_sequence | tr "\t" "\n" | grep -c .)
            store_result_ID=$(echo $targetID"_vs_"$original_genome_ID"_"$original_chr_ID)

            ## HERE BE DRAGONS!!! FIX REQUIRED!!
            other_ref_genome_file=$(ls $work_folder"/"$masked_genomes"/"$targetID"/"*".masked")

            reference_fasta_minimap=$local_work_folder"/"$targetID"_reference.fasta"
            if [ ! -f $reference_fasta_minimap ]
            then
              extract_sequence=$(echo $matching_sequence | sed "s/$targetID//" | sed "s/_RagTag$//" | sed "s/^_//" )
              seqkit grep -j $threads -p $extract_sequence $other_ref_genome_file | seqkit seq -u | sed "s/>.*/$rename_fasta_seq/" > $reference_fasta_minimap
            fi

            target_fasta_minimap=$local_work_folder"/"$store_result_ID"_target.fasta"
            seqkit grep -j $threads -p $matching_sequence $ragtag_output_folder"/ragtag.scaffold.fasta" | sed "s/>.*/$rename_fasta_seq/" > $target_fasta_minimap
            do_ragtag_syri
          fi
        fi
      done
      echo $original_chr_ID >> $work_folder"/"$storage"/Already_done.txt"
    fi

    count_regions=$(($count_regions +1))
  done

  code_block_completed
  slurm_block_completed
  echo ""
  echo "Special Note"
  echo "This block of code generates a lot of sam files that are redundant. The summary block will delete them."
  echo "Known issues:"
  echo "1) Samtools can be unreliable sorting big files. If the computer runs out of memory the run will be aborted"
  echo "Solution: brute force"
  echo "Try again any sample that failed to sort itself with more memory requested and without the -@ "'$NPROCS'" parameter"
  echo "A failed run can be identifued by the presence of extra files with '.bam.tmp' in the name, or by a failure to produce an index"

  echo "2) Minimap might also run out of memory. This block compares ch/romosomes one to one for that reason. But there is no guarantee that a really big one will not break it."
  echo "You can know it failed because the samtools sort will fail."
  echo "Solution: more brute force"
  echo "You will need to contact the people in charge of the cluster and request for more memory."
fi

################################################################################

if [ $REGION_FINDER_Re == "TRUE" ]
then
  confirm_input_file
  confirm_region_file
  storage=Extracted_Regions_of_Interest_Syri
  syri_synteny=RagTag_Plus_Syri_Region
  masked_genomes=REPEAT_MASKER_RUN
  override_repeatmasker=FALSE

  query_group_ragtag=RAGTAG_Scaffolding
  reference_group_ragtag=RAGTAG_Scaffolding_References

  start_block_code

  query_assemblies=$(awk -F "\t" '{if ($2=="Q") print $1}' $genome_table | sort)
  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $genome_table | sort)

  count_regions=2
  traverse_regions_of_interest=$(grep -c . $region_table)
  while [ $count_regions -le $traverse_regions_of_interest ]
  do
    # First things first, lets get the data and be sure this regions is not already processed.
    region_ID=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $1}')
    check_already_done=$(grep -w -F -c $region_ID $work_folder"/"$storage"/Already_done.txt")
    region_coord_adjusted=$(echo .)

    if [ $check_already_done -gt 0 ]
    then
      echo "Information for "$region_ID" already in the folder"
      echo "Skipping..."
    else
      # Now to retrieve the data for the run and tell the user
      genome_ID=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $2}')
      Group_ID=$(awk -F "\t" -v genome=$genome_ID '{if ($1==genome) print $2}' $genome_table)

      if [ $override_repeatmasker == "TRUE" ]
      then
        genome_fasta=$(awk -F "\t" -v wanted=$genome_ID '{if ($1==wanted) print $3}' $genome_table )
      else
        genome_fasta=$(ls $work_folder"/"$masked_genomes"/"$genome_ID"/"*".masked")
      fi

      chr_ID=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $3}')
      Ref_Anchor_point=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $4}')
      Anchor_location=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $5}')
      Extension=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $6}')
      Ref_Anchor_Tolerance=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $7}')
      Error_Margin=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $8}')
      Known_Anchor_Fasta=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $9}')

      if [ -d $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID ]
      then
        rm -rf $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID
      fi

      mkdir $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID
      mkdir $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Sequences"
      mkdir $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS"
      mkdir $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Syri_BLOCKS/Feature_Details"
      mkdir $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST"
      mkdir $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/References"

      echo "Retrieving information for: "$region_ID
      echo "Desired Region: "$region_ID > $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "Starting point on Genome: "$genome_ID" ("$Group_ID") -- Chr: "$chr_ID" -- Anchor Point: "$Ref_Anchor_point >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "Relative Ancrhor location: "$Anchor_location" -- Aproximate missmatch: "$Ref_Anchor_Tolerance >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "Expeted region extension: "$Extension" -- Set Error Margin: "$Error_Margin >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "Fasta File with reference Anchors: "$Known_Anchor_Fasta  >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"

      if [ $Anchor_location == "Start" ]
      then
        region_Start_End=$(echo $Ref_Anchor_point" "$Extension" "$Ref_Anchor_Tolerance | awk -F " " '{print $1-$3","$1+$2}')
      elif [ $Anchor_location == "End" ]
      then
        region_Start_End=$(echo $Ref_Anchor_point" "$Extension" "$Ref_Anchor_Tolerance | awk -F " " '{print $1-$2","$1+$3}')
      elif [ $Anchor_location == "Middle" ]
      then
        echo "Note: Anchor is in the Middle. Remember that column 7 does not apply, and that the Extension is devided by 2 (rounded to the lower integer)" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
        echo "" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"

        region_range=$(echo $Extension | awk -F " " '{printf "%.0f\n", $1/2 }')
        region_Start_End=$(echo $Ref_Anchor_point" "$region_range | awk -F " " '{print $1-$2","$1+$2}')
      else
        echo "......"
        echo "How the hell did we even got here?!"
        echo "Somehow Col 5 on the region file is broken and skipped the controls -____-"
        exit
      fi

      echo "Target_Assembly Putative_Start Putative_End N_Blocks_On_Target N_Blocks_Off_Target" | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Regions.txt"

      if [ $Known_Anchor_Fasta == "NA" ]
      then
        echo "Target_Assembly Synteny_block Origin_Coord_Start Origin_Coord_End Target_Coord_Start Target_Coord_End N_Match_Contigs Match_Contigs Overlap_Lower_Boundary Overlap_Target Overlap_Upper_Boundary" | tr " " "\t" > $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.txt"
      else
        if [ -f $Known_Anchor_Fasta ]
        then
          anchor_information=$(seqkit seq -n $Known_Anchor_Fasta | awk -F " " '{print $1"_("$2")"}' | tr "\n" " ")
          echo "Target_Assembly Synteny_block Origin_Coord_Start Origin_Coord_End Target_Coord_Start Target_Coord_End N_Match_Contigs Match_Contigs Overlap_Lower_Boundary Overlap_Target Overlap_Upper_Boundary "$anchor_information | tr " " "\t" > $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Matching_Blocks.txt"
        else
          echo "No anchor file found -.-''"
          echo "Check: "$Known_Anchor_Fasta
          echo "Region: "$region_ID
          exit
        fi
      fi

      echo "Original_Genome Target_Genome RagTag_Strand Correction Study_Region Lower_Boundary Target_Region Upper_Boundary"  | tr " " "\t" >>  $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Coordinate_corrections.txt"

      fix_coord=$region_Start_End
      fix_out_of_bounds_coords
      region_Start_End=$fix_coord

      if [ $Error_Margin -eq 0 ]
      then
        region_lower_boundary=NA
        region_upper_boundary=NA

        search_start=$(echo $region_Start_End | awk -F "," '{print $1}')
        search_end=$(echo $region_Start_End | awk -F "," '{print $2}')
      else
        region_lower_boundary=$(echo $region_Start_End","$Error_Margin | awk -F "," '{print $1-$3","$1-1}')

        fix_coord=$region_lower_boundary
        fix_out_of_bounds_coords
        region_lower_boundary=$fix_coord
        region_upper_boundary=$(echo $region_Start_End","$Error_Margin | awk -F "," '{print $2+1","$2+$3}')
        fix_coord=$region_upper_boundary
        fix_out_of_bounds_coords
        region_upper_boundary=$fix_coord

        # These are going to be used to filter the matching regions. The ranges for upper and lower boundaries are meant for reporting
        search_start=$(echo $region_lower_boundary | awk -F "," '{print $1}')
        search_end=$(echo $region_upper_boundary | awk -F "," '{print $2}')
      fi

      echo "---------" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "Work Coordinates: "$search_start":"$search_end >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "Lower margin boundary: "$region_lower_boundary | tr "," ":" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "Actual query area: "$region_Start_End  | tr "," ":" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "Upper margin boundary: "$region_upper_boundary  | tr "," ":" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "Were coords out of bounds? ("$region_coord_adjusted")" | sed "s/(\.)/No/" | sed "s/(X)/Yes/" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "---------" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
      echo "" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"

      echo "Extracting control sequences..."
      rename_seq=$(echo ">"$get_sequence" Coords:"$search_start":"$search_end)

      seqkit grep -j $threads -p "$chr_ID" $genome_fasta | seqkit subseq -j $threads -r $search_start":"$search_end | seqkit seq -u | sed "s/>.*/$rename_seq/" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Sequences/Control_Target_sequences.fasta"
      seqkit grep -j $threads -p "$chr_ID" $genome_fasta | seqkit seq -u >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Sequences/Full_sequence_control.fasta"

      if [ ! $Known_Anchor_Fasta == "NA" ]
      then
        echo "Anchors provided. Verifying positions in the original region"
        makeblastdb -in $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Sequences/Control_Target_sequences.fasta" -dbtype nucl -parse_seqids -input_type fasta -out $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/References/"$genome_ID"_Control_ref"
        echo "Done with reference for BLAST."
        echo ""
        seqkit seq -g -M 30 $Known_Anchor_Fasta >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/Temp_File_short.fasta"
        seqkit seq -g -m 31 $Known_Anchor_Fasta >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/Temp_File_long.fasta"

        check_short_blast=$(grep -c ">" $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/Temp_File_short.fasta")
        check_long_blast=$(grep -c ">" $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/Temp_File_long.fasta")
        echo "We have "$check_short_blast" short anchors and "$check_long_blast" long ones. Now for the BLAST"

        if [ $check_short_blast -gt 0 ]
        then
          echo "Short BLAST..."
          blastn -task blastn-short -query $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/Temp_File_short.fasta" -db $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/References/"$genome_ID"_Control_ref" -outfmt "6 std" -out $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/"$genome_ID"_Control_short_anchor.blastn" -perc_identity $anchor_blast_ident -qcov_hsp_perc $anchor_blast_qcov -num_threads $threads
        fi

        if [ $check_long_blast -gt 0 ]
        then
          blastn -query $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/Temp_File_long.fasta" -db $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/References/"$genome_ID"_Control_ref" -outfmt "6 std" -out $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Anchor_BLAST/"$genome_ID"_Control_long_anchor.blastn" -perc_identity $anchor_blast_ident -qcov_hsp_perc $anchor_blast_qcov -num_threads $threads
        fi
      fi

      echo "Initial adjustments completed. Looking for matches"
      echo "Matches for "$genome_ID"("$Group_ID")"

      if [ $Group_ID == "R" ]
      then
        search_ragtag_id=$(echo $genome_ID"_"$chr_ID"_RagTag")

        wanted_start=$search_start
        wanted_end=$search_end

        start_col=2
        end_col=3
        Notal_in_query=6

        extract_sec_from_ragtag=TRUE
        for targetID in $query_assemblies
        do
          echo $targetID
          ragtag_output=$work_folder"/"$query_group_ragtag"/Ref_"$genome_ID"/"$targetID"_ragtag_output"
          syri_output=$work_folder"/"$syri_synteny"/"$chr_ID"/Syri/"$genome_ID"_vs_"$targetID"_"$chr_ID"_syri.out"

          echo $targetID" Scaffolded" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
          extract_syntenic_blocks
        done

        ragtag_output_path=$work_folder"/"$reference_group_ragtag
        ref_to_ref_Comparison=$(ls -d $ragtag_output_path"/"*"_"$genome_ID"_"* | grep -v .txt$ | sed "s/\/$//" | sed "s/.*\///" | sed "s/_ragtag_output$//" )

        for ref_to_ref in $ref_to_ref_Comparison
        do
          siri_out_id=$(echo $ref_to_ref | sed "s/Ref_//g")
          check_first=$(echo $siri_out_id | sed "s/_vs_/\n/" | grep -w -F -n $genome_ID | awk -F ":" '{print $1}')

          ragtag_output=$ragtag_output_path"/"$ref_to_ref"_ragtag_output"
          syri_output=$work_folder"/"$syri_synteny"/"$chr_ID"/Syri/"$siri_out_id"_"$chr_ID"_syri.out"

          # Because in some comparisons the target boundaries are going to be updated, I need to re-set them every time.
          region_Start_End=$(sed -n 11p $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"  | sed "s/.* //" | tr ":" ",")

          if [ $region_upper_boundary != "NA" ]
          then
            region_lower_boundary=$(sed -n 10p $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt" | sed "s/.* //" | tr ":" ",")
            region_upper_boundary=$(sed -n 12p $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt" | sed "s/.* //" | tr ":" ",")
          fi

          if [ $check_first -eq 1 ]
          then
            extract_sec_from_ragtag=TRUE
            targetID=$(echo $siri_out_id | sed "s/_vs_/\n/" | tail -n 1)
            search_ragtag_id=$(echo $genome_ID"_"$chr_ID"_RagTag")

            wanted_start=$search_start
            wanted_end=$search_end

            start_col=2
            end_col=3
            Notal_in_query=6

            echo $targetID" Scaffolded" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
            extract_syntenic_blocks
          else
            targetID=$(echo $siri_out_id | sed "s/_vs_/\n/" | head -n 1)

            echo $targetID" Reference" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
            extract_sec_from_ragtag=FALSE

            if [ $override_repeatmasker == "TRUE" ]
            then
              target_genome_file=$(awk -F "\t" -v wanted=$targetID '{if ($1==wanted) print $3}' $genome_table )
            else
              target_genome_file=$(ls $work_folder"/"$masked_genomes"/"$targetID"/"*".masked")
            fi

            search_ragtag_id=$(grep -w -F $genome_ID"_"$chr_ID $ragtag_output"/ragtag.scaffold.agp" | awk -F "\t" '{print $1}' | sort -u)
            correction=$(grep -w -F $genome_ID"_"$chr_ID $ragtag_output"/ragtag.scaffold.agp" | awk -F "\t" '{print $2-$7}')
            ragtag_strand=$(grep -w -F $genome_ID"_"$chr_ID $ragtag_output"/ragtag.scaffold.agp" | awk -F "\t" '{print $9}')

            # Fix search area
            fix_start=$search_start
            fix_end=$search_end
            query_to_ref_correction
            wanted_start=$result_start
            wanted_end=$result_end

            # Fix original target area
            fix_start=$(echo $region_Start_End | awk -F "," '{print $1}')
            fix_end=$(echo $region_Start_End | awk -F "," '{print $2}')
            query_to_ref_correction
            region_Start_End=$(echo $result_start","$result_end)

            # Fix boundaries
            if [ $region_upper_boundary != "NA" ]
            then
              fix_start=$(echo $region_lower_boundary | awk -F "," '{print $1}')
              fix_end=$(echo $region_lower_boundary | awk -F "," '{print $2}')
              query_to_ref_correction
              region_lower_boundary=$(echo $result_start","$result_end)

              fix_start=$(echo $region_upper_boundary | awk -F "," '{print $1}')
              fix_end=$(echo $region_upper_boundary | awk -F "," '{print $2}')
              query_to_ref_correction
              region_upper_boundary=$(echo $result_start","$result_end)
            fi

            echo $genome_ID" "$targetID" "$ragtag_strand" "$correction" "$wanted_start","$wanted_end" "$region_lower_boundary" "$region_Start_End" "$region_upper_boundary | tr " " "\t" | tr "," ":" >>  $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Coordinate_corrections.txt"

            start_col=7
            end_col=8
            Notal_in_query=1
            # Find the issue here
            extract_syntenic_blocks
          fi
        done
      else
        echo ""
        echo "--------------"
        echo "The original assembly with the region is designated as: Query"
        echo "Be advised: this portion of the script recieved less debugging than its counter part."
        echo "While it should work as a Ref_1 vs Ref_2 where Ref_2 has the region, I didn't expend nearly as much time working on it."
        echo 'Therefore, a healthy dose of "Trust but verify" is advised.'
        echo "--------------"
        echo ""
        ref_ids=$(ls -d $work_folder"/"$query_group_ragtag"/"* | grep -v .txt$ | sed "s/\/$//" | sed "s/.*\///" | sed "s/^Ref_//" )

        for targetID in $ref_ids
        do
          echo $targetID" Original" >> $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt"
          ragtag_output=$work_folder"/"$query_group_ragtag"/Ref_"$targetID"/"$genome_ID"_ragtag_output"
          syri_output=$work_folder"/"$syri_synteny"/"$chr_ID"/Syri/"$targetID"_vs_"$genome_ID"_"$chr_ID"_syri.out"

          if [ -f $syri_output ]
          then
            extract_sec_from_ragtag=FALSE

            if [ $override_repeatmasker == "TRUE" ]
            then
              target_genome_file=$(awk -F "\t" -v wanted=$targetID '{if ($1==wanted) print $3}' $genome_table )
            else
              target_genome_file=$(ls $work_folder"/"$masked_genomes"/"$targetID"/"*".masked")
            fi

            search_ragtag_id=$(grep -w -F $genome_ID"_"$chr_ID $ragtag_output"/ragtag.scaffold.agp" | awk -F "\t" '{print $1}' | sort -u)
            correction=$(grep -w -F $genome_ID"_"$chr_ID $ragtag_output"/ragtag.scaffold.agp" | awk -F "\t" '{print $2-$7}')
            ragtag_strand=$(grep -w -F $genome_ID"_"$chr_ID $ragtag_output"/ragtag.scaffold.agp" | awk -F "\t" '{print $9}')

            if [ $region_upper_boundary != "NA" ]
            then
              region_lower_boundary=$(sed -n 10p $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt" | sed "s/.* //" | tr ":" ",")
              region_upper_boundary=$(sed -n 12p $work_folder"/"$storage"/"$region_ID"_from_"$genome_ID"_"$Group_ID"/Summary_region.txt" | sed "s/.* //" | tr ":" ",")
            fi

            # Fix search area
            fix_start=$search_start
            fix_end=$search_end
            query_to_ref_correction
            wanted_start=$result_start
            wanted_end=$result_end

            # Fix original target area
            fix_start=$(echo $region_Start_End | awk -F "," '{print $1}')
            fix_end=$(echo $region_Start_End | awk -F "," '{print $2}')
            query_to_ref_correction
            region_Start_End=$(echo $result_start","$result_end)

            # Fix boundaries
            if [ $region_upper_boundary != "NA" ]
            then
              fix_start=$(echo $region_lower_boundary | awk -F "," '{print $1}')
              fix_end=$(echo $region_lower_boundary | awk -F "," '{print $2}')
              query_to_ref_correction
              region_lower_boundary=$(echo $result_start","$result_end)

              fix_start=$(echo $region_upper_boundary | awk -F "," '{print $1}')
              fix_end=$(echo $region_upper_boundary | awk -F "," '{print $2}')
              query_to_ref_correction
              region_upper_boundary=$(echo $result_start","$result_end)
            fi

            start_col=7
            end_col=8
            Notal_in_query=1
            extract_syntenic_blocks
          fi
        done
      fi
      echo $region_ID >> $work_folder"/"$storage"/Already_done.txt"
    fi

    count_regions=$(($count_regions + 1))
  done

  code_block_completed

  echo "Special Notes:"
  echo "1) Given the complexity observed in some of these regions I'm unable to automate the final refinement of the region of interest."
  echo "At least not in a way where I can ensure success."
  echo "The best I can do is to give you information on what is on every synteny block and where. And with that information at hand allow you to refinse the search."
  echo ""
  echo "2) Going from a region defined on a query assembly to the references is not supported beyond this point due to lack of real need. Plus some of them have proven to be affected by missasemblies on key regions."
  echo "If you have to use this feature, use the results here to guide an inspection of your reference of choice."
  echo ""
  echo "On an ideal scenario, the desired region will be located on a single synteny block on both genomes. An after identifying it it will be a simple matter of translate coordinates to one region to the other."
  echo "Regions like PTC2 however where there is a huge fragmentation and potential re-arragement of segments, they either hate to all be considered or a few of them handpicked."
  echo ""
  echo "Known issues:"
  echo "1) Because I'm lazy, at one moment in the reporting this code counts synteny blocks out of the expected mark by counting the sequence headers with '_Miss_' generated."
  echo "If one of your assemblies IDs includes the string _Miss_ every block will be reported as out of its intended place. "
  echo "Solution: Don't name an assembly _Miss_"
  echo ""
  echo "2) Due to decision latter down the line, the same Region ID described from two differnt assemblies are not compatible."
  echo "Solution: Ensure it is not the case."
fi

################################################################################

if [ $FULL_ANCHOR_BLAST == "TRUE" ]
then
  confirm_input_file
  confirm_region_file
  storage=Full_Anchor_Blasts
  masked_genomes=REPEAT_MASKER_RUN

  clean_run_references=FALSE # Remaking the references all the time is tiresome. But it is not worth the time to make this also recyclable. If you add a new assembly sewt to true.
  if [ -d $work_folder"/"$storage ]
  then
    rm -rf $work_folder"/"$storage
  fi

  mkdir $work_folder"/"$storage

  echo "This block of code was added out of despair"
  echo "Tests with PTC2 gave inconsistent results and I plan to confirm what is going on."
  echo "The idea is to first confirm if the marker sequences are in the genome or not to begin with."
  echo "At least for PTC2 they are"

  genomes_IDs=$(awk -F "\t" '{print $1}' $genome_table | sed 1d)
  regions_with_anchor=$(sed 1d $region_table | awk -F "\t" '{if ($9!="NA") print $1}' )
  if [ $clean_run_references == "TRUE" ]
  then
    if [ -d $work_folder"/"$storage"/REFS" ]
    then
      rm -rf $work_folder"/"$storage"/REFS"
    fi
    mkdir $work_folder"/"$storage"/REFS"
    for genome in $genomes_IDs
    do
      genome_file=$(awk -F "\t" -v genome=$genome '{if (genome==$1) print $3}' $genome_table)
      makeblastdb -in $genome_file -dbtype nucl -parse_seqids -input_type fasta -out $work_folder"/"$storage"/REFS/"$genome"_unmasked_ref"

      masked_genome_file=$(ls $work_folder"/"$masked_genomes"/"$genome"/"*".masked")
      makeblastdb -in $masked_genome_file -dbtype nucl -parse_seqids -input_type fasta -out $work_folder"/"$storage"/REFS/"$genome"_masked_ref"
    done
  fi

  for region_id in $regions_with_anchor
  do
    mkdir $work_folder"/"$storage"/"$region_id"_Unmasked"
    mkdir $work_folder"/"$storage"/"$region_id"_Masked"

    anchor_file=$(awk -F "\t" -v region=$region_id '{if ($1==region) print $9}' $region_table)

    for genome in $genomes_IDs
    do
      echo $region_id" --- "$genome" ||| "$anchor_file

      blastn -task blastn-short -query $anchor_file -db $work_folder"/"$storage"/REFS/"$genome"_unmasked_ref" -outfmt "6 std" -out $work_folder"/"$storage"/"$region_id"_Unmasked/"$genome".blastn" -perc_identity $anchor_blast_ident -qcov_hsp_perc $anchor_blast_qcov -num_threads $threads
      blastn -task blastn-short -query $anchor_file -db $work_folder"/"$storage"/REFS/"$genome"_masked_ref" -outfmt "6 std" -out $work_folder"/"$storage"/"$region_id"_Masked/"$genome".blastn" -perc_identity $anchor_blast_ident -qcov_hsp_perc $anchor_blast_qcov -num_threads $threads
    done
  done
  echo "Note: This miniblock currently only does blastn-short"
fi

################################################################################

if [ $SLURM_SYRI_PLOTSR_Re == TRUE ]
then
  # CURRENT ISSUE: the genome ID and the allocation of not aligned blocks in the bed files.
  confirm_input_file
  confirm_region_file
  storage=REGION_PLOTSR
  summary_information=Extracted_Regions_of_Interest_Syri
  syri_synteny=RagTag_Plus_Syri_Region

  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "PLOTSR_"$current_date)

  start_block_code

  count_regions=2
  traverse_regions_of_interest=$(grep -c . $region_table)
  while [ $count_regions -le $traverse_regions_of_interest ]
  do
    # First things first, lets get the data and be sure this regions is not already processed.
    region_ID=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $1}')
    check_already_done=$(grep -w -F -c $region_ID $work_folder"/"$storage"/Already_done.txt")
    if [ $check_already_done -gt 0 ]
    then
      echo "Already done. Skipping: "$region_ID" ($count_regions -le $traverse_regions_of_interest)"
    else
      echo "Plotsr: "$region_ID
      if [ -d $work_folder"/"$storage"/"$region_ID ]
      then
        rm -rf $work_folder"/"$storage"/"$region_ID
      fi

      mkdir $work_folder"/"$storage"/"$region_ID
      mkdir $work_folder"/"$storage"/"$region_ID"/Genome_Files"
      mkdir $work_folder"/"$storage"/"$region_ID"/Marker_Locations"

      original_genome_ID=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $2}')
      original_chr=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $3}')
      original_genome_group_ID=$(awk -F "\t" -v genome=$original_genome_ID '{if (genome==$1) print $2}' $genome_table)

      work_genome_ids=$(sed 1d $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_genome_group_ID"/Matching_Regions.txt" | awk -F "\t" '{print $1}')
      for target_id in $work_genome_ids
      do
        check_first=$(awk -F "\t" '{print $1}' $work_folder"/"$syri_synteny"/"$original_chr"/Run_Registrty.txt" | grep -F $original_genome_ID | grep -F $target_id"_" | grep -c ^$original_genome_ID)
        check_second=$(awk -F "\t" '{print $1}' $work_folder"/"$syri_synteny"/"$original_chr"/Run_Registrty.txt" | grep -F $original_genome_ID |  grep -F $target_id"_" | grep -c "_"$original_genome_ID)

        echo $region_ID" "$original_genome_ID"_vs_"$target_id" ("$check_first"/"$check_second")"

        if [ $check_first -gt 0 ]
        then
          plot_this_syri=$work_folder"/"$syri_synteny"/"$original_chr"/Syri/"$original_genome_ID"_vs_"$target_id"_"$original_chr"_syri.out"
          reference_file=$(awk -F "\t" -v run_id=$original_genome_ID"_vs_"$target_id"_"$original_chr '{if ($1==run_id) print $2}' $work_folder"/"$syri_synteny"/"$original_chr"/Run_Registrty.txt")
          target_file=$(awk -F "\t" -v run_id=$original_genome_ID"_vs_"$target_id"_"$original_chr '{if ($1==run_id) print $3}' $work_folder"/"$syri_synteny"/"$original_chr"/Run_Registrty.txt")
          chr_name=$(seqkit seq -n $reference_file)

          region_coords=$(grep -F "Work Coordinates:" $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_genome_group_ID"/Summary_region.txt" | tr -d " " | awk -F ":" '{print $2" "$3}' )
          write_plotsr_files
        elif [ $check_second -gt 0 ]
        then
          plot_this_syri=$work_folder"/"$syri_synteny"/"$original_chr"/Syri/"$target_id"_vs_"$original_genome_ID"_"$original_chr"_syri.out"
          reference_file=$(awk -F "\t" -v run_id=$target_id"_vs_"$original_genome_ID"_"$original_chr '{if ($1==run_id) print $2}' $work_folder"/"$syri_synteny"/"$original_chr"/Run_Registrty.txt")
          target_file=$(awk -F "\t" -v run_id=$target_id"_vs_"$original_genome_ID"_"$original_chr '{if ($1==run_id) print $3}' $work_folder"/"$syri_synteny"/"$original_chr"/Run_Registrty.txt")
          chr_name=$(seqkit seq -n $reference_file)

          region_coords=$(awk -F "\t" -v genome=$target_id '{if ($1==genome) print $2" "$3}' $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_genome_group_ID"/Matching_Regions.txt")
          write_plotsr_files
        else
          echo "Cannot find syri blocks between "$original_genome_ID" and "$TargetID
          echo "Reconsider your life choices, listen to The Sojourn Audio Drama, and check if you have all the data needed for this"
          exit
        fi
      done
      echo $region_ID >> $work_folder"/"$storage"/Already_done.txt"
    fi
    count_regions=$(($count_regions + 1))
  done
  code_block_completed
  slurm_block_completed
fi

################################################################################

if [ $UNMASKED_SCAFFOLDING_Re == TRUE ]
then
  echo "WARNING!!!"
  echo "This block can take up forever and I have not yet found a way to reliably run it at OSU computers in the background."

  confirm_input_file
  storage=UNMASKED_SCAFFOLDING
  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Scaffolding_References_"$current_date)

  start_block_code
  if [ ! -d $work_folder"/"$storage"/Ref_vs_Query" ]
  then
    mkdir $work_folder"/"$storage"/Ref_vs_Query"
    mkdir $work_folder"/"$storage"/Ref_vs_Ref"
  fi

  query_group_ragtag=RAGTAG_Scaffolding
  reference_group_ragtag=RAGTAG_Scaffolding_References

  query_assemblies=$(awk -F "\t" '{if ($2=="Q") print $1}' $genome_table | sort)
  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $genome_table | sort)

  for ref_ID in $reference_assemblies
  do
    echo "Re-scaffolding Queries relative to: "$ref_ID
    for query_ID in $query_assemblies
    do
      check_done=$(grep -c -w -F $ref_ID"_vs_"$query_ID $work_folder"/"$storage"/Already_done.txt")

      if [ $check_done -gt 0 ]
      then
        echo $ref_ID"_vs_"$query_ID" already present. Skipping"
      else
        ragtag_output_folder=$work_folder"/"$query_group_ragtag"/Ref_"$ref_ID"/"$query_ID"_ragtag_output"
        new_out_file=$(echo $work_folder"/"$storage"/Ref_vs_Query/"$ref_ID"_vs_"$query_ID"_unmasked_scaffolds")
        targetID=$query_ID

        non_masked_scaffolding
        echo $ref_ID"_vs_"$query_ID >> $work_folder"/"$storage"/Already_done.txt"
      fi
    done
  done

  ragtag_output_path=$work_folder"/"$reference_group_ragtag
  ref_to_ref_Comparison=$(ls -d $ragtag_output_path"/"*"_ragtag_output" | sed "s/\/$//" | sed "s/.*\///" | sed "s/_ragtag_output$//" )

  for ref_to_ref in $ref_to_ref_Comparison
  do
    check_done=$(grep -c -w -F $ref_to_ref $work_folder"/"$storage"/Already_done.txt")

    if [ $check_done -gt 0 ]
    then
      echo $ref_to_ref" already present. Skipping"
    else
      print_ref=$(echo $ref_to_ref | sed "s/Ref_//g")
      targetID=$(echo $print_ref | sed "s/_vs_/\n/" | tail -n 1 )
      ragtag_output_folder=$ragtag_output_path"/"$ref_to_ref"_ragtag_output"

      new_out_file=$(echo $work_folder"/"$storage"/Ref_vs_Ref/"$print_ref"_unmasked_scaffolds")

      non_masked_scaffolding

      echo $ref_to_ref >> $work_folder"/"$storage"/Already_done.txt"
    fi
  done

  code_block_completed
  echo "Known issues:"
  echo "1) I did a dirty trick to simplify the code, the first sequence the script tries to check is named SKIP_ME. Please do not name any conting like that or the output file might break."
  echo "2) For simplicity, this block doesn't cleans a partial run where the scaffolding was interrupted half way. If that is the case you will need to fix it manually by removing the files for that assembly, and the Temp folder"
fi

################################################################################

if [ $SLURM_ANNOTATION_TRANSFER_Re == TRUE ]
then
  echo "Prepare for take off!!"
  echo "Except we wont. Liftoff is a pain, it cannot run in paralell from the same folder. So I either make over 101 folders for submission to SLURM and make you launch them one by one, or make you run it from a single computer "
  echo "Sadly for you (and me as I write this) the latter should take over 2 weeks to complete with the current datasets. Unless I disable the polishing step and convince everyone to work with suboptimal annotations"
  echo "So prepare for a ride... each individual command needs to be sent to the cluster manually, unless I get express permission to use a script."

  confirm_input_file
  storage=Annotation_Transfer_Liftoff
  unmasked_scaffolds=UNMASKED_SCAFFOLDING

  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Annotation_Transfer_"$current_date)

  start_block_code

  if [ ! -d $work_folder"/"$storage"/Individual_Transfers" ]
  then
    mkdir $work_folder"/"$storage"/Individual_Transfers"
  fi

  if [ ! -d $work_folder"/"$storage"/Temp_Files" ]
  then
    mkdir $work_folder"/"$storage"/Temp_Files"
    mkdir $work_folder"/"$storage"/Temp_Files/Sequences"
    mkdir $work_folder"/"$storage"/Temp_Files/Annotations"

    local_count=1
    while [ $local_count -le 5 ]
    do
      mkdir $work_folder"/"$storage"/Temp_Files/Annotations/Round_"$local_count
      local_count=$(($local_count +1))
    done
  fi

  if [ -d $cluster_submission_folder"/Liftoff_runs" ]
  then
    rm -rf $cluster_submission_folder"/Liftoff_runs"
  fi

  mkdir $cluster_submission_folder"/Liftoff_runs"
  mkdir $cluster_submission_folder"/Liftoff_runs/Initial_run"
  mkdir $cluster_submission_folder"/Liftoff_runs/Bulk_run"

  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $genome_table | sort)
  #command_number=1
  cluster_batch_job=1
  round_storage=1

  echo "Going after Ref vs Queries"
  really_first=TRUE
  target_locations=$work_folder"/"$unmasked_scaffolds"/Ref_vs_Query"
  prepare_for_liftoff

  echo "Going after Ref vs Ref"
  really_first=FALSE
  N_run_fule=2
  target_locations=$work_folder"/"$unmasked_scaffolds"/Ref_vs_Ref"
  prepare_for_liftoff

  # This is a bit redundant, but this part transfers the annotation from the worst references to the best one
  # It might be worth checking eventually

  echo "Final transfers"
  N_assemblies_ref=$(awk -F "\t" '{if ($2=="R" && $4!=NA) print $1}' $genome_table | grep -c .)
  if [ $N_assemblies_ref -gt 1 ]
  then

    final_ref_assemblies=$(awk -F "\t" '{if ($2=="R") print $1"\t"$15}' $work_folder"/Assembly_Stats/Assemblies_stats.txt" | sort -r -n -k2 | awk -F "\t" '{print $1}')
    for RefID1 in $final_ref_assemblies
    do
      ref_genome_1=$(awk -F "\t" -v genome=$RefID1 '{if (genome==$1) print $3}' $genome_table)
      round_storage=1

      if [ ! -f $work_folder"/"$storage"/Temp_Files/Sequences/"$RefID1"_genome.fasta" ]
      then
        echo "Copying and indexing ref genome"
        cp $target_genome $work_folder"/"$storage"/Temp_Files/Sequences/"$RefID1"_genome.fasta"
        samtools faidx $work_folder"/"$storage"/Temp_Files/Sequences/"$RefID1"_genome.fasta"
      fi

      gff_file_1=$(awk -F "\t" -v genome=$RefID1  '{if (genome==$1) print $4}' $genome_table )
      if [ ! -f $work_folder"/"$storage"/Temp_Files/Annotations/"$RefID1"_annotation.gff" ]
      then
        echo "Copying annotation file"
        cp $gff_file_1 $work_folder"/"$storage"/Temp_Files/Annotations/"$RefID1"_annotation.gff"
      fi

      for RefID2 in $final_ref_assemblies
      do
        if [ ! $RefID1 == $RefID2 ]
        then
          echo $RefID1" to "$RefID2
          check_final_run=$(grep -w -F Final_Ref_to_Ref $work_folder"/"$storage"/Already_done.txt" | grep -w -F $RefID1 | grep -w -F -c $RefID2)

          if [ $check_final_run -gt 0 ]
          then
            echo $RefID1" vs "$RefID2" Already Done. Skip!!"
          else
            ref_genome_2=$(awk -F "\t" -v genome=$RefID2 '{if (genome==$1) print $3}' $genome_table)

            liftoff_cluster_submissions=$cluster_submission_folder"/Liftoff_runs/Bulk_run"
            echo "liftoff -p "'$NPROCS'" -dir "$work_folder"/"$storage"/Temp_Files/"$RefID1"transfer_to_"$RefID2"_temp_files -o "$work_folder"/"$storage"/Individual_Transfers/"$RefID1"_transfered_to_"$RefID2"_annotation.gff -u "$work_folder"/"$storage"/Individual_Transfers/"$RefID1"_transfered_to_"$RefID2"_unnmaped_features.gff -g "$work_folder"/"$storage"/Temp_Files/Annotations/Round_"$round_storage"/"$RefID1"_annotation.gff -copies -polish -cds "$work_folder"/"$storage"/Temp_Files/Sequences/"$RefID2"_genome.fasta" $work_folder"/"$storage"/Temp_Files/Sequences/"$RefID1"_genome.fasta" >> $liftoff_cluster_submissions"/Slurm_Command_batch_Final.sh"

            round_storage=$(($round_storage + 1))
            echo "Final_Ref_to_Ref "$RefID1" "$RefID2 >> $work_folder"/"$storage"/Already_done.txt"
          fi
        fi
      done
    done
  fi

  code_block_completed
  slurm_block_completed

  echo "Special notes: "
  echo "1) The Liftoff runs are going to share some input files generated by the program: the GFF, one per reference, and the first thing it is going to do is to create an index for it."
  echo "I refuse to investigate what would happen if multiple runs try to create the same file. If they are going to crash, overwrite each other or what not."
  echo "Therefore the submission is devided in 2. The first file has 1 work per reference. It should create all the needed indexes without conflicts with one another."
  echo "The second file has the rest of them. To avoid headaches make sure the first one completes before launching any other."
  echo ""
  echo "And again, due to poor implementation of Liftoff you will have to go into almost a hundred of folders and launch them into the cluster. Manually"
  echo "Try to do only 20 at the time so you can check if you need to remove some temporary files at "$work_folder"/"$storage"/Temp_Files"
  echo "DO NOT REMOVE Annotations and Sequences!!!! Some other code blocks will check on them if they use software with the same requirements."
  echo ""
  echo "2) After you are done with Liftoff check the folders where you stored the sequences. Particularly the one with the one with the unmasked scaffolds. If I'm not mistaken the minimap alignments generated during Liftoff are stored wherever the sequences happen to be."
  echo ""

  echo "Known issues:"
  echo "1) I was in a rush with this one, the script will break if you have more than 5 references. Call Javier Calvelo if that becomes a problem, the solution is somewhat trivial, just increase the number of copies of the annotation. A REAL solution is not trivial at all (at least for me)"
  echo "2) I didn't managed to solve the conflict between multiple runs of liftoff completly. The remaining error is a silent bug where some of the annotation transfers generate a GFF with only comments."
  echo "The solution is to check and re-run those lines alone."
  echo "3) The separation between Initial and Bulk runs were a mistake that is not worth addressing at this moment. It adds nothing and under some conditions it might contain runs that are not compatible to run toghether. If it happens just do it one at the time"
fi

################################################################################

if [ $CONFIRM_LIFTOFF_RUNS == TRUE ]
then
  echo "Because lifftof is so annoying to run, a code block is needed to ensure no transfer is missing. And if so, which one"
  storage=Annotation_Transfer_Liftoff

  liftoff_run_dir=$work_folder"/Cluster_Submissions/Liftoff_runs"
  if [ -d $work_folder"/"$storage"/Control_Runs" ]
  then
    rm -rf $work_folder"/"$storage"/Control_Runs"
  fi
  mkdir $work_folder"/"$storage"/Control_Runs"

  cat $liftoff_run_dir"/Initial_run/"*"sh" | grep . >> $work_folder"/"$storage"/Control_Runs/Run_commands.sh"
  cat $liftoff_run_dir"/Bulk_run/"*"sh" | grep . >> $work_folder"/"$storage"/Control_Runs/Run_commands.sh"

  count=1
  walker=$(grep -c . $work_folder"/"$storage"/Control_Runs/Run_commands.sh")

  while [ $count -le  $walker ]
  do
    echo "Checking file Nº "$count" of "$walker
    out_file=$(sed -n $count"p" $work_folder"/"$storage"/Control_Runs/Run_commands.sh" | tr " " "\n" | grep -A 1 \\-o | tail -n 1)
    if [ ! -f $out_file ]
    then
      echo "File Missing: "
      echo $out_file
    else
      check_data=$(head $out_file | grep -c .)
      if [ $check_data -lt 10 ]
      then
        echo "File too small, might need a re-run: "
        echo $out_file
      fi
    fi
    count=$(($count + 1))
  done
  echo "Done!!"
  echo "If all you saw was a count down then everything is ok"
  echo "I think..."
fi

################################################################################

if [ $EVALUATE_ANOT_TRANSFER == TRUE ]
then
  confirm_input_file
  confirm_region_file
  storage=Annotation_Transfer_Liftoff

  if [ ! -d $work_folder"/"$storage"/Sumary_stats" ]
  then
    mkdir $work_folder"/"$storage"/Sumary_stats"
  fi

  reference_assemblies=$(awk -F "\t"  '{if ($2=="R") print $1}' $genome_table)

  # there is no need for complex manipulation. Let's just count how many genes there were and how many were transfered and under what conditions
  echo "RefID Num_genes" | tr " " "\t" >  $work_folder"/"$storage"/Sumary_stats/Ref_gene_counts.txt"
  for refID in $reference_assemblies
  do
    if [ ! -f $work_folder"/"$storage"/Temp_Files/Annotations/"$refID"_annotation.gff" ]
    then
      echo "Missing gff at: "$work_folder"/"$storage"/Temp_Files/Annotations/"$refID"_annotation.gff"
      echo "Remember that this script relies on this temporary files to work. Check for errors while listenning to Faithless. A sojourn story at Nebula."
      exit
    fi

    N_genes=$(awk -F "\t" '{if ($3=="gene") print $9}' $work_folder"/"$storage"/Temp_Files/Annotations/"$refID"_annotation.gff" | awk -F ";" '{print $1}' | sort -u | grep -c .)
    echo $refID" "$N_genes | tr " " "\t" >> $work_folder"/"$storage"/Sumary_stats/Ref_gene_counts.txt"
  done

  echo "TransferID Reference N_unique_genes N_genes_good_cov N_genes_cluster" | tr " " "\t" > $work_folder"/"$storage"/Sumary_stats/Transfers_gene_counts.txt"
  for refID in $reference_assemblies
  do
    all_transfers=$(ls $work_folder"/"$storage"/Individual_Transfers/"*"_polished" | sed "s/.*_transfered_to_//" | sed "s/_annotation.gff_polished//" | sort -u)
    for transfer in $all_transfers
    do
      if [ -f $work_folder"/"$storage"/Individual_Transfers/"$refID"_transfered_to_"$transfer"_annotation.gff_polished" ]
      then
        echo "Working on: "$refID" to "$transfer
        awk -F "\t" '{if ($3=="gene") print $9}' $work_folder"/"$storage"/Individual_Transfers/"$refID"_transfered_to_"$transfer"_annotation.gff_polished" > $work_folder"/"$storage"/Sumary_stats/current_genes.temp"

        N_Genes_transfered_non_cluster=$(grep -c ";extra_copy_number=0" $work_folder"/"$storage"/Sumary_stats/current_genes.temp" )
        N_Genes_cluster=$(grep -v ";extra_copy_number=0" $work_folder"/"$storage"/Sumary_stats/current_genes.temp" | tr ";" "\n" |grep "Name=" | sort -u | grep -c .)
        N_genes_cov_good=$(grep ";extra_copy_number=0" $work_folder"/"$storage"/Sumary_stats/current_genes.temp"  | tr ";" "\n" | grep "coverage=" | awk -F "=" '{if ($2>=0.95) print}' | grep -c .)

        echo $transfer" "$refID" "$N_Genes_transfered_non_cluster" "$N_genes_cov_good" "$N_Genes_cluster | tr " " "\t" >> $work_folder"/"$storage"/Sumary_stats/Transfers_gene_counts.txt"
      fi
    done
  done
fi

################################################################################

if [ $RAW_SYN_BLOCK_TARGET_REGIONS_ANNOTATION == TRUE ]
then
  confirm_input_file
  confirm_region_file

  storage=Raw_Syn_Block_Inspection_Mk2
  summary_information=Extracted_Regions_of_Interest_Syri
  Repeat_information=REPEAT_MASKER_RUN
  query_group_ragtag=RAGTAG_Scaffolding
  reference_group_ragtag=RAGTAG_Scaffolding_References
  annotation_transfers=Annotation_Transfer_Liftoff

  start_block_code

  echo "Extracting gathered information from each block"
  echo "Note: this will only work if the annotation files are in GFF3 format"

  count_regions=2
  traverse_regions_of_interest=$(grep -c . $region_table)

  ref_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $work_folder"/Assembly_Stats/Assemblies_stats.txt"  | sort)
  columns_on_header=NA
  for Ref in $ref_assemblies
  do
    columns_on_header=$(echo $columns_on_header" N_Genes_in_Block_"$Ref" N_Genes_on_Target_"$Ref" Per_genes_found_"$Ref" Per_genes_Other_"$Ref" Target_Gene_IDs_"$Ref" N_transposition_Candidates"$Ref)
    columns_on_missing=$(echo $columns_on_missing" N_Genes_in_Block_"$Ref)
  done
  columns_on_header=$(echo $columns_on_header | sed "s/NA //")
  columns_on_missing=$(echo $columns_on_missing | sed "s/NA //")

  while [ $count_regions -le $traverse_regions_of_interest ]
  do
    # First things first, lets get the data and be sure this regions is not already processed.
    region_ID=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $1}')
    check_already_done=$(grep -w -F -c $region_ID $work_folder"/"$storage"/Already_done.txt")

    if [ $check_already_done -gt 0 ]
    then
      echo "Information for "$region_ID" already in the folder"
      echo "Skipping..."
    else
      echo "Working on "$region_ID
      if [ -d $work_folder"/"$storage"/"$region_ID ]
      then
        rm -rf $work_folder"/"$storage"/"$region_ID
      fi

      original_genome_ID=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $2}' )
      catch_up_ref_assemblies=$(echo $ref_assemblies | tr " " "\n" | grep -v $original_genome_ID)

      original_chr=$(sed -n $count_regions"p" $region_table | awk -F "\t" '{print $3}' )

      original_group=$(awk -F "\t" -v genome=$original_genome_ID '{if (genome==$1) print $2}' $genome_table)
      original_annotation=$(awk -F "\t" -v genome=$original_genome_ID '{if (genome==$1) print $4}' $genome_table)

      original_work_coords=$(grep -F "Work Coordinates:" $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_group"/Summary_region.txt" | awk -F " " '{print $3}')

      start_search_area=$(echo $original_work_coords | awk -F ":" '{print $1}')
      end_search_area=$(echo $original_work_coords | awk -F ":" '{print $2}')

      mkdir $work_folder"/"$storage"/"$region_ID
      mkdir $work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures"
      mkdir $work_folder"/"$storage"/"$region_ID"/Target_Assemblies"
      mkdir $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes"
      mkdir $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Missing_Blocks"
      mkdir $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks"

      echo "AssemblyID SynBlockID Start End Block_Lenght Region_Boundaries "$columns_on_header | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Synteny_Block_content.tab"
      echo "AssemblyID SynBlockID Start End Block_Lenght "$columns_on_missing | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Missing_Blocks_on_targets.tab"

      if [ ! $original_annotation == NA ]
      then
        echo "Getting original annotation"
        awk -F "\t" -v chr=$original_chr -v end=$end_search_area '{if ($1==chr && $3=="gene" && $4 <= end) print}' $original_annotation | awk -F "\t" -v start=$start_search_area '{if (start <= $5) print}' >> $work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_original_genes.gff"
        awk -F "\t" '{print $9}' $work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_original_genes.gff" | tr ";" "\n" | grep ^ID= | sed "s/ID=//" >> $work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_original_genes.ID"

        ref_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $work_folder"/Assembly_Stats/Assemblies_stats.txt" | grep -w -F -v $original_genome_ID)
        for refID in $ref_assemblies
        do
          awk -F "\t" -v chr=$original_chr -v end=$end_search_area  '{if ($1==chr && $3=="gene" && $4 <= end) print}' $work_folder"/"$annotation_transfers"/Individual_Transfers/"$refID"_transfered_to_"$original_genome_ID"_annotation.gff_polished" | awk -F "\t" -v start=$start_search_area '{if (start <= $5) print}' >> $work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_trans_"$refID"_genes.gff"
          awk -F "\t" '{print $9}' $work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_trans_"$refID"_genes.gff" | tr ";" "\n" | grep ^ID= | sed "s/ID=//" >> $work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_trans_"$refID"_genes.ID"
        done
      else
        ref_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $work_folder"/Assembly_Stats/Assemblies_stats.txt"  | sort)
        for refID in $ref_assemblies
        do
          awk -F "\t" -v chr=$original_chr -v end=$end_search_area  '{if ($1==chr && $3=="gene" && $4 <= end) print}' $work_folder"/"$annotation_transfers"/Individual_Transfers/"$refID"_transfered_to_"$original_genome_ID"_annotation.gff_polished" | awk -F "\t" -v start=$start_search_area '{if (start <= $5) print}' >> $work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_trans_"$refID"_genes.gff"
          awk -F "\t" '{print $9}' $work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_trans_"$refID"_genes.gff" | tr ";" "\n" | grep ^ID= | sed "s/ID=//" >> $work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_trans_"$refID"_genes.ID"
        done
      fi

      ref_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $work_folder"/Assembly_Stats/Assemblies_stats.txt" | sort)

      count_synth_block=2
      walker_syn_block=$(grep -c . $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_group"/Matching_Blocks.txt")

      while [ $count_synth_block -le $walker_syn_block ]
      do
        TargetID=$(sed -n $count_synth_block"p" $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_group"/Matching_Blocks.txt" | awk -F "\t" '{print $1}')
        target_group=$(awk -F "\t" -v genome=$TargetID '{if (genome==$1) print $2}' $genome_table)

        Syn_ID=$(sed -n $count_synth_block"p" $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_group"/Matching_Blocks.txt" | awk -F "\t" '{print $2}')
        check_origin_region=$(grep -F $TargetID" " $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_group"/Summary_region.txt" | awk -F " " '{print $2}')
        Boundary_information=$(sed -n $count_synth_block"p" $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_group"/Matching_Blocks.txt" | awk -F "\t" '{print $9"-"$10"-"$11}')

        Syn_start=$(sed -n $count_synth_block"p" $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_group"/Matching_Blocks.txt" | awk -F "\t" '{print $5}')
        Syn_end=$(sed -n $count_synth_block"p" $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_group"/Matching_Blocks.txt" | awk -F "\t" '{print $6}')
        Syn_lenght=$(($Syn_end - $Syn_start + 1))

        if [ ! -d $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID ]
        then
          mkdir $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID
        fi

        # Add oc  movement of checking genes on target. I needed to fix a special case so I moved it here
        if [ $check_origin_region == "Scaffolded" ]
        then
          transfer_gff_files=$(ls $work_folder"/"$annotation_transfers"/Individual_Transfers/"*"_to_"$original_genome_ID"_vs_"$TargetID"_annotation.gff_polished" | sort)
          do_genes_on_target
        else
          target_annotation=$(awk -F "\t" -v genome=$TargetID '{if (genome==$1) print $4}' $genome_table)
          transfer_gff_files=$(ls $work_folder"/"$annotation_transfers"/Individual_Transfers/"*"_transfered_to_"$TargetID"_annotation.gff_polished")
          do_genes_on_target
        fi

        if [ $check_origin_region == "Scaffolded" ]
        then
          ragtag_scaffolded_seq=$(echo $original_genome_ID"_"$original_chr"_RagTag")
          transfer_gff_files=$(ls $work_folder"/"$annotation_transfers"/Individual_Transfers/"*"_to_"$original_genome_ID"_vs_"$TargetID"_annotation.gff_polished" | sort)

          if [ ! $Syn_start == "NA" ]
          then
            # repeat_recovery
            echo $ragtag_scaffolded_seq" Syn_Block gene "$Syn_start" "$Syn_end" . + . ID="$Syn_ID | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID"/Syn_blocks.gff"

            transfer_genes_registry=NA
            for transfer in $transfer_gff_files
            do
              transfer_id=$(echo $transfer | sed "s/.*\///" | sed "s/_transfered_to_/\n/" | head -n 1 )
              echo "Checking Annotation from "$transfer_id" transfered to "$TargetID" ("$Syn_ID")"

              awk -F "\t" -v chr=$ragtag_scaffolded_seq -v end=$Syn_end '{if ($1==chr && $3=="gene" && $4 <= end) print}' $transfer | awk -F "\t" -v start=$Syn_start '{if (start <= $5) print}' >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".gff"
              awk -F "\t" '{print $9}' $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".gff" | tr ";" "\n" | grep ^ID= | sed "s/ID=//" >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".ID"
              cat $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".gff" >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID"/Genes_on_synblock_"$Syn_ID"_from_"$transfer_id".gff"

              check_found_genes_file=$work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".ID"

              genes_in_block_counting
              transfer_genes_registry=$(echo $transfer_genes_registry" "$N_genes_in_block" "$N_Genes_on_target_zone" "$Per_found" "$Per_other" "$Genes_on_target_zone" "$N_Transposintion_candidates)
            done

            transfer_genes_registry=$(echo $transfer_genes_registry | sed "s/NA //")

            echo $TargetID" "$Syn_ID" "$Syn_start" "$Syn_end" "$Syn_lenght" "$Boundary_information" "$transfer_genes_registry | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Synteny_Block_content.tab"
          else
            # Summary of blocks not found in the target
            missing_block_start=$(sed -n $count_synth_block"p" $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_group"/Matching_Blocks.txt" | awk -F "\t" '{print $3}')
            missing_block_end=$(sed -n $count_synth_block"p" $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_group"/Matching_Blocks.txt" | awk -F "\t" '{print $4}')
            missing_lenght=$(($missing_block_end - $missing_block_start + 1))

            registro_missing_genes=NA
            for transfer_id in $ref_assemblies
            do
              if [ $transfer_id == $original_genome_ID ]
              then
                transfer=$work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_original_genes.gff"
                check_genes_file=$work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_original_genes.ID"
              else
                transfer=$work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_trans_"$refID"_genes.gff"
                check_genes_file=$work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_trans_"$transfer_id"_genes.ID"
              fi

              awk -F "\t" -v chr=$original_chr -v end=$missing_block_end '{if ($1==chr && $3=="gene" && $4 <= end) print}' $transfer | awk -F "\t" -v start=$missing_block_start '{if (start <= $5) print}' >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Missing_Blocks/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".gff"
              awk -F "\t" '{print $9}' $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Missing_Blocks/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".gff" | tr ";" "\n" | grep ^ID= | sed "s/ID=//" >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Missing_Blocks/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".ID"

              N_genes_in_block=$(grep -c . $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Missing_Blocks/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".ID")
              registro_missing_genes=$(echo $registro_missing_genes" "$N_genes_in_block)
            done

            registro_missing_genes=$(echo $registro_missing_genes | sed "s/NA //")
            echo $TargetID" "$Syn_ID" "$missing_block_start" "$missing_block_end" "$missing_lenght" "$registro_missing_genes | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Missing_Blocks_on_targets.tab"
          fi
        else

          echo "Ref_assembly: "$TargetID
          ragtag_output=$(ls -d $work_folder"/"$reference_group_ragtag"/Ref_"* | grep "_"$original_genome_ID"_" | grep "_"$TargetID"_" | awk '{print $0"/ragtag.scaffold.agp"}')

          ragtag_scaffolded_seq=$(awk -F "\t" -v get_seq=$original_genome_ID"_"$original_chr '{if ($6==get_seq) print $1}' $ragtag_output  | sed "s/^$TargetID//" |  sed "s/^_//" | sed "s/_RagTag//") # Depende de quien sea la referencia en la comparación.
          target_annotation=$(awk -F "\t" -v genome=$TargetID '{if (genome==$1) print $4}' $genome_table)
          transfer_gff_files=$(ls $work_folder"/"$annotation_transfers"/Individual_Transfers/"*"_"$TargetID"_annotation.gff_polished")

          if [ ! $Syn_start == "NA" ]
          then
            # repeat_recovery
            echo $ragtag_scaffolded_seq" Syn_Block gene "$Syn_start" "$Syn_end" . + . ID="$Syn_ID | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID"/Syn_blocks.gff"
            transfer_genes_registry=NA
            for transfer_id in $ref_assemblies
            do
              if [ $transfer_id == $TargetID ]
              then
                if [ ! $target_annotation == "NA" ] && [ $target_group == "R" ]
                then
                  # Recover info from Original
                  awk -F "\t" -v chr=$ragtag_scaffolded_seq -v end=$Syn_end '{if ($1==chr && $3=="gene" && $4 <= end) print}' $target_annotation | awk -F "\t" -v start=$Syn_start '{if (start <= $5) print}' >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$TargetID".gff"
                  awk -F "\t" '{print $9}' $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$TargetID".gff" | tr ";" "\n" | grep ^ID= | sed "s/ID=//" >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$TargetID".ID"
                  cat $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$TargetID".gff" >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID"/Genes_on_synblock_"$Syn_ID"_from_"$TargetID".gff"

                  check_found_genes_file=$work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$TargetID".ID"
                  check_genes_file=$work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_trans_"$TargetID"_genes.ID"

                  if [ ! -f $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID"/Target_genes_from_"$TargetID".gff" ]
                  then
                    grep -w -F -f $check_genes_file $target_annotation  | awk -F "\t" '{if ($3=="gene") print}' >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID"/Target_genes_from_"$TargetID".gff"
                  fi
                  genes_in_block_counting
                fi
              else
                transfer=$(ls $work_folder"/"$annotation_transfers"/Individual_Transfers/"$transfer_id"_transfered_to_"$TargetID"_annotation.gff_polished")

                awk -F "\t" -v chr=$ragtag_scaffolded_seq -v end=$Syn_end '{if ($1==chr && $3=="gene" && $4 <= end) print}' $transfer | awk -F "\t" -v start=$Syn_start '{if (start <= $5) print}' >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".gff"
                awk -F "\t" '{print $9}' $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".gff" | tr ";" "\n" | grep ^ID= | sed "s/ID=//" >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".ID"
                cat $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".gff" >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID"/Genes_on_synblock_"$Syn_ID"_from_"$transfer_id".gff"
                check_found_genes_file=$work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Syn_Blocks_Genes/"$TargetID"_"$region_ID"_"$Syn_ID"_genes_from_"$transfer_id".ID"

                if [ $transfer_id == $original_genome_ID ]
                then
                  check_genes_file=$work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_original_genes.ID"
                else
                  check_genes_file=$work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_trans_"$transfer_id"_genes.ID"
                fi

                if [ ! -f $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID"/Target_genes_from_"$transfer_id".gff" ]
                then
                  grep -w -F -f $check_genes_file $transfer | awk -F "\t" '{if ($3=="gene") print}' >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$TargetID"/Target_genes_from_"$transfer_id".gff"
                fi
                genes_in_block_counting
              fi
              transfer_genes_registry=$(echo $transfer_genes_registry" "$N_genes_in_block" "$N_Genes_on_target_zone" "$Per_found" "$Per_other" "$Genes_on_target_zone" "$N_Transposintion_candidates)
            done
            transfer_genes_registry=$(echo $transfer_genes_registry | sed "s/NA //")
            echo $TargetID" "$Syn_ID" "$Syn_start" "$Syn_end" "$Syn_lenght" "$Boundary_information" "$transfer_genes_registry | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Synteny_Block_content.tab"
          fi
        fi
        count_synth_block=$(($count_synth_block + 1))
      done

      echo "Done with Blocks. Now for alternative tracks on other  "

      # Yes dear code reader. This is an add on feature that I forgot to include in the first version of the script, and it is easier to code after all else is done.
      all_targets=$(awk -F "\t" '{print $1}' $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_group"/Matching_Blocks.txt" | sed 1d | sort -u )

      for Target in $all_targets
      do
        check_origin_region=$(grep -F $Target" " $work_folder"/"$summary_information"/"$region_ID"_from_"$original_genome_ID"_"$original_group"/Summary_region.txt" | awk -F " " '{print $2}')
        if [ $check_origin_region == "Scaffolded" ]
        then
          if [ ! -d $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$Target"/Additional_Scaffoldings" ]
          then
            mkdir $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$Target"/Additional_Scaffoldings"
          fi

          for catch in $catch_up_ref_assemblies
          do
            if [ ! $catch == $Target ]
            then
              for ref_data in $ref_assemblies
              do
                gff_thing=$(echo $work_folder"/"$annotation_transfers"/Individual_Transfers/"$ref_data"_transfered_to_"$catch"_vs_"$Target"_annotation.gff_polished")
                if [ -f $gff_thing ]
                then
                  echo "Chechiking scaffolding: "$catch"_vs_"$Target" || Anot: "$ref_data
                  if [ $ref_data == $original_genome_ID ]
                  then
                    region_genes_ID=$work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_original_genes.ID"
                  else
                    region_genes_ID=$work_folder"/"$storage"/"$region_ID"/Original_Genome_Fetures/"$region_ID"_from_"$original_genome_ID"_trans_"$ref_data"_genes.ID"
                  fi
                  cat $region_genes_ID | awk '{print "ID="$1";"}' >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$Target"/Additional_Scaffoldings/Matching_Genes_"$catch"_vs_"$Target"_Anot_"$ref_data".get"

                  grep -F -f $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$Target"/Additional_Scaffoldings/Matching_Genes_"$catch"_vs_"$Target"_Anot_"$ref_data".get" $gff_thing >> $work_folder"/"$storage"/"$region_ID"/Target_Assemblies/Display_tracks/"$Target"/Additional_Scaffoldings/Matching_Genes_"$catch"_vs_"$Target"_Anot_"$ref_data".gff"
                fi
              done
            fi
          done
        fi
      done
      echo $region_ID >> $work_folder"/"$storage"/Already_done.txt"
    fi
    count_regions=$(($count_regions + 1 ))
  done

  code_block_completed

  echo "Remember that the results of this block are meant to guide the manual delimitation of the target areas. For most cases everything up to this point will be an overkill."
  echo "You will get huge synteny blocks filled with genes and you will basically pick a few markers and work with the surrounding area. With the RAGTAG results helping by putting all the fragments in order."
  echo "The fun begins with more complex regions where significant re-arragements took place. In those cases the synteny block information is meant to help guide what happened with every gene."
  echo "Where did they go? Do you need to increase the studdy area or can you narrow it?"
  echo ""
  echo "The output tables Synteny_Block_content.tab and Missing_Blocks_on_targets.tab should give you some idea of what you are working on for each assembly. Which blocks have the most genes, which ones seem to be out of target by some re-arragement and wich ones cover far more than the target sequence."
  echo "Then, within the folder Display_tracks you will find many *.gff files with information: Location of the synteny blocks, genes found in the synteny blocks, and location of the target genes found in the reference assembly where the region was defined."
  echo "With these GFFs you can go to any genome visualization software and visually inspect what is going on at the gene level. What you want to see is the annotation for the same gene overlapping in the last two. If there are genes missing in one you need to see if they are entirely missing or located in another place."
  echo ""
  echo "And remember: Paranoia is the key to success in genomics. If something is off check if it might not be a bug in the pipeline that was left unadressed."
  echo "Only then can we start talking about changing the world"
  echo ""
  echo "Known issues:"
  echo "Ragtag is an all or nothing process when it assings contigs. This becomes an issue when the reference is split over the desired region and contigs get asigned to one end or the other. You will find such cases by looking at the Matching_Blocks.txt file and see that an assembly has only NOTAL blocks asinged and their display track folder will only have genes on target."
  echo "on this section of the code, such samples will be missing. If this happens the plan is for you to identify them and define a new region based on the genes."
fi

################################################################################

if [ $RAW_SYN_BLOCK_TARGET_REGIONS_ADDITIONAL_TRACKS == "TRUE" ]
then
  confirm_input_file
  confirm_region_file
  storage=Raw_Syn_Block_Inspection_Mk2
  query_group_ragtag=RAGTAG_Scaffolding
  reference_group_ragtag=RAGTAG_Scaffolding_References
  scaffolded=UNMASKED_SCAFFOLDING
  annotation_transfers=Annotation_Transfer_Liftoff

  if [ -d $work_folder"/"$storage"/Additional_tracks" ]
  then
    rm -rf $work_folder"/"$storage"/Additional_tracks"
  fi

  # Original contig location
  mkdir $work_folder"/"$storage"/Additional_tracks"
  scafolded_dir=$work_folder"/"$scaffolded"/Ref_vs_Query"
  generate_contig_gff

  scafolded_dir=$work_folder"/"$scaffolded"/Ref_vs_Ref"
  generate_contig_gff

  # Genes on special contigs
  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $genome_table | sort)
  contig=JAKZJK010000215.1
  source=BS90_bu_et_al_qtl
  get_special_contig_genes
fi

################################################################################

if [ $SLURM_RAW_SYN_BLOCK_TARGET_REGIONS_ADDITIONAL_TRACKS_SYRI == TRUE ]
then
  echo "You are trying to run: RAW_SYN_BLOCK_TARGET_REGIONS_ADDITIONAL_TRACKS_SYRI"
  echo "This block requires a complete 'confirm_final_region_file' input file to run"
  echo "The coordinates on the queries do not matter and they can be a place-holder, but this script needs:"
  echo "1) Scaffolded file for each query"
  echo "2) ID of the reference assembly (it will be compared with the one in the name of the file, the script will run only if they do not match)"
  echo "3) Accuarate scaffold ID on the query"
  echo ""
  echo "The region location is taken from the 'confirm_final_region_file'"
  echo "just to protect my sanity every block asigned will be included in the output gffs. There original script proved to be an overkill."

  confirm_input_file
  confirm_region_file
  confirm_final_region_file

  storage=Additional_Syri_tracks

  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Additional_syri_"$current_date)
  query_group_ragtag=RAGTAG_Scaffolding
  reference_group_ragtag=RAGTAG_Scaffolding_References
  masked_genomes=REPEAT_MASKER_RUN
  local_work_folder=$work_folder"/"$storage

  start_block_code

  if [ ! -d $local_work_folder"/Minimap" ]
  then
    mkdir $local_work_folder"/Inputs"
    mkdir $local_work_folder"/Minimap"
    mkdir $local_work_folder"/Syri"
  fi

  original_genome_file=$(ls $work_folder"/"$masked_genomes"/"$original_descriptor"/"*".masked")
  count_final_regions=2
  while [ $count_final_regions -le $total_regions ]
  do
    is_scaffold=$(sed -n $count_final_regions"p" $final_region_recovery | awk -F "\t" '{print $5}')
    scaffold_reference=$(sed -n $count_final_regions"p" $final_region_recovery | awk -F "\t" '{print $6}')

    if [ ! $scaffold_reference == $original_descriptor ] && [ $is_scaffold == "TRUE" ]
    then
      # Reference info
      region_ID=$(sed -n $count_final_regions"p" $final_region_recovery | awk -F "\t" '{print $1}')
      rename_fasta_seq=$(echo ">"$region_ID"_sequence")
      original_chr=$(awk -F "\t" -v region=$region_ID -v ref=$original_descriptor '{if ($1==region && $2==ref) print $3}' $region_table)

      reference_fasta_minimap=$local_work_folder"/Inputs/"$region_ID"_ref.fasta"
      target_fasta_minimap=$local_work_folder"/Inputs/"$region_ID"_"$assembly_ID"_ref.fasta"

      # Preparing input files
      if [ ! -f $reference_fasta_minimap ]
      then
        seqkit grep -p $original_chr $original_genome_file | sed "s/>.*/$rename_fasta_seq/" > $reference_fasta_minimap
      fi

      # Query info
      assembly_ID=$(sed -n $count_final_regions"p" $final_region_recovery | awk -F "\t" '{print $2}')
      assembly_group=$(awk -F "\t" -v genome=$assembly_ID '{if (genome==$1) print $2}' $genome_table)
      query_genome_sequence_ID=$(sed -n $count_final_regions"p" $final_region_recovery | awk -F "\t" '{print $4}' | awk -F ":" '{print $1}')

      if [ $assembly_group == "Q" ]
      then
        ragtag_output_folder=$work_folder"/"$query_group_ragtag"/Ref_"$scaffold_reference"/"$assembly_ID"_ragtag_output"
      else
        ragtag_output_folder=$work_folder"/"$reference_group_ragtag"/Ref_"$scaffold_reference"_vs_Ref_"$assembly_ID"_ragtag_output"
      fi

      store_result_ID=$(echo $region_ID"__Scaffolding_"$scaffold_reference"_vs_"$assembly_ID)
      seqkit grep -p $query_genome_sequence_ID $ragtag_output_folder"/ragtag.scaffold.fasta" | sed "s/>.*/$rename_fasta_seq/" > $target_fasta_minimap
      do_ragtag_syri

      echo $region_ID" "$scaffold_reference" "$assembly_ID" "$original_chr" "$query_genome_sequence_ID | tr " " "\t" >> $local_work_folder"/Syri/Registry.txt"
    fi

    count_final_regions=$(($count_final_regions + 1))
  done
fi

################################################################################

if [ $SLURM_RAW_SYN_BLOCK_TARGET_REGIONS_ADDITIONAL_TRACKS_GFF == TRUE ]
then
  confirm_input_file
  confirm_region_file
  confirm_final_region_file

  storage=Additional_Syri_tracks/GFF_Format
  syri_out=Additional_Syri_tracks/Syri

  if [ -d $work_folder"/"$storage ]
  then
    rm -rf $work_folder"/"$storage
  fi

  start_block_code

  if [ -f $work_folder"/"$syri_out"/Registry.txt" ]
  then
    count=1
    walker=$(grep -c . $work_folder"/"$syri_out"/Registry.txt")

    while [ $count -le $walker ]
    do
      echo $count" // "$walker

      syri_information=$(sed -n $count"p" $work_folder"/"$syri_out"/Registry.txt")

      region_ID=$(echo $syri_information | awk -F " " '{print $1}')
      Reference_ID=$(echo $syri_information | awk -F " " '{print $2}')
      Query_ID=$(echo $syri_information | awk -F " " '{print $3}')
      Target_seq=$(echo $syri_information | awk -F " " '{print $5}')

      awk -F "\t" -v chr=$Target_seq '{if ($10=="-" && $7!="-") print chr" Syn_Block gene "$7" "$8" . + . ID="$9}' $work_folder"/"$syri_out"/"$region_ID"__Scaffolding_"$Reference_ID"_vs_"$Query_ID"_syri.out" | sed "s/ID=NOTAL/ID=M_NOTAL/" | tr " " "\t" >> $work_folder"/"$storage"/Syri_Blocks__"$region_ID"_"$Reference_ID"_vs_"$Query_ID"_syri.out"
      count=$(($count + 1))
    done

  else
    echo "Cannot find the Registry file at: "$work_folder"/"$syri_out"/Registry.txt"
    exit
  fi
fi

################################################################################

if [ $SLURM_GENE_SEQ_INTERPRO_Re == TRUE ]
then
  echo "NOTE: This block will be deprecated because the local version is not up to date. Instead I have to use singularity to run the analysis and that beheaves somewhat differently and might prove hard to just automate"
  echo "Inestead of running this on Slurm, modify the script Utility_up_to_date_Interpro.sh"

  storage=Ref_Genes_Annotation

  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Gene_Annotation_"$current_date)

  confirm_input_file
  start_block_code

  if [ ! -d $work_folder"/"$storage"/Interpro_scan" ]
  then
    mkdir $work_folder"/"$storage"/Interpro_scan"
  fi

  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $genome_table | sort)
  for RefID in $reference_assemblies
  do
    echo "Working on: "$RefID
    check_already_done=$(grep -w -F -c $RefID"_interproscan" $work_folder"/"$storage"/Already_done.txt")

    if [ $check_already_done -gt 0 ]
    then
      echo "Already done"
      echo "Skipping..."
    else
      if [ -d $work_folder"/"$storage"/Interpro_scan/"$RefID ]
      then
        rm -rf $work_folder"/"$storage"/Interpro_scan/"$RefID
      fi

      mkdir $work_folder"/"$storage"/Interpro_scan/"$RefID

      prot_file=$(awk -F "\t" -v genome=$RefID '{if (genome==$1) print $5}' $genome_table)
      if [ $prot_file != "NA" ]
      then
        echo "interproscan.sh -cpu "'$NPROCS'" -f TSV -pa -b "$work_folder"/"$storage"/Interpro_scan/"$RefID" -exclappl SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE -dra -i "$prot_file >> $cluster_submission_folder"/"$submit_to_cluster"_Interpro.txt"
      else
        echo "No Protein file for this reference. Skipping"
      fi
      echo $RefID"_interproscan" >> $work_folder"/"$storage"/Already_done.txt"
    fi
  done
  code_block_completed
  slurm_block_completed
fi

################################################################################

if [ $SLURM_INTERPRO_PER_GENE_Re == TRUE ]
then
  # I didn't wanted to make this. I just wanted to add a simple loci ID in front of the transcript ID on Interpro...
  # But... 6.1G for one assembly? And I need to potentially check it multiple times? Unaceptable.
  # So I'm gonna push my limited bash skills to the limits and abuse SLURM just a little bit.

  storage=Ref_Genes_Annotation
  confirm_input_file

  if [ ! -d  $work_folder"/"$storage"/Work_version" ]
  then
    mkdir $work_folder"/"$storage"/Work_version"
    mkdir $work_folder"/"$storage"/Work_version/Interpro"
    mkdir $work_folder"/"$storage"/Work_version/Interpro/split_temp"
    echo "Already done" >> $work_folder"/"$storage"/Work_version/Already_done_interpro.txt"
  fi

  if [ -d $cluster_submission_folder"/Add_small_name_big_stuff" ]
  then
    rm -rf $cluster_submission_folder"/Add_small_name_big_stuff"
  fi

  mkdir $cluster_submission_folder"/Add_small_name_big_stuff"
  mkdir $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage"

  annotated_assemblies=$(ls $work_folder"/"$storage"/Interpro_scan/"*".tsv" | sed "s/\/$//" | sed "s/.*\///" | sed "s/.tsv//")
  for RefID in $annotated_assemblies
  do
    check_done=$(grep -w -F -c $RefID $work_folder"/"$storage"/Work_version/Already_done_interpro.txt")
    if [ $check_done -gt 0 ]
    then
      echo $RefID" already re-formated. Skipping"
    else
      echo "Working on: "$RefID

      current_ref=$RefID
      get_prefix_info_for_gff

      gff_file=$(awk -F "\t" -v genome=$RefID '{if (genome==$1) print $4}' $genome_table )

      # Oh? You dirty boy has over a million lines to work with? with each itteration becoming slower and more memory draining than the last?
      # That's cute
      # The idea here is to split the original file into chunks of 1000 lines and work through them one by one. This way there are never more than 1000 lines in memory at any time

      echo "Spliting the original file into workable chunks..."
      split -l 50000 $work_folder"/"$storage"/Interpro_scan/"$RefID".tsv" $work_folder"/"$storage"/Work_version/Interpro/split_temp/"$RefID"_"
      echo "Done. Writing SLURM scripts..."

      all_splits=$(ls $work_folder"/"$storage"/Work_version/Interpro/split_temp/"$RefID"_"*)
      N_splits=$(echo $all_splits | tr " " "\n" | grep -c .)
      current_split=1
      for splat in $all_splits
      do
        echo $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh" >> $cluster_submission_folder"/Add_small_name_big_stuff/Run_me.sh"
        #######################################################
        ################## Script for SLURM ###################
        #######################################################
        echo "#!/usr/bin/env bash" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo ""  >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo "echo "'"'"LociID ProteinID Sequence_MD5_digest Sequence_Lenght Analysis SignatureID Signature_Description Start End Score Status Date InterproID Interpro_Description GO_TERMS Pathways"'"'" | tr "'"'" "'"'' "\t"'" >> "$work_folder"/"$storage"/Work_version/Interpro/"$RefID"_work_table_"$current_split".tsv" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo "current_protID=PLACE_HOLDER_DONT_USE_ME" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo "count=1"  >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo "walker="'$'"(grep -c . "$splat")" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo "while [ "'$'"count -le "'$'"walker ]" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo "do" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo " protID="'$'"(sed -n "'$'"count"'"'"p"'"' $splat" | awk -F "'"'"\\t"'"' "'"'{print $1}'"')" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo " if [ ! "'$'"current_protID == "'$'"protID ]" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo " then" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo "  gene_loci="'$''(grep -F '$prefix_extract'$'"protID"'"'$subfix_extract'"'" "$gff_file" | awk -F "'"'"\\t"'"' "'"'{print $9}'"' | head -n 1 | tr "'"'";"'"' '"'"\\n"'"'" | grep -F "$take_item' | sed "s/'$take_item'//")' >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo "  current_protID="'$'"protID" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo " fi" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo "  sed -n "'$'"count"'"'"p"'"' $splat" | awk -F "'"'"\\t"'"'" -v loci="'$'"gene_loci ""'{print loci"'"'"\\t"'"''$0'"}' >> "$work_folder"/"$storage"/Work_version/Interpro/"$RefID"_work_table_"$current_split".tsv" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo "  count="'$'"(("'$'"count + 1))" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        echo "done" >> $cluster_submission_folder"/Add_small_name_big_stuff/Script_storage/Add_gene_ID_"$RefID"_"$current_split".sh"
        ###################################################
        ###################################################

        current_split=$(($current_split + 1))
      done
      echo $RefID >> $work_folder"/"$storage"/Work_version/Already_done_interpro.txt"
    fi
  done
  code_block_completed
  slurm_block_completed
  echo ""
  echo "Note: This block is a bit difference. It generates the scripts that are going to be run in SLURM. Before running them you need to give them execution permits"
fi

################################################################################

if [ $SLURM_EGGNOG_PER_GENE_Re == TRUE ]
then
  storage=Ref_Genes_Annotation
  confirm_input_file

  if [ ! -d  $work_folder"/"$storage"/Work_version/Eggnog" ]
  then
    mkdir $work_folder"/"$storage"/Work_version/Eggnog"
    echo "Already done" >> $work_folder"/"$storage"/Work_version/Already_done_eggnog.txt"
  fi

  annotated_assemblies=$(ls $work_folder"/"$storage"/Eggnog_mapper/"*".emapper.annotations" | sed "s/\/$//" | sed "s/.*\///" | sed "s/.emapper.annotations//")
  for RefID in $annotated_assemblies
  do
    check_done=$(grep -w -F -c $RefID $work_folder"/"$storage"/Work_version/Already_done_eggnog.txt")
    if [ $check_done -gt 0 ]
    then
      echo $RefID" already Done. Skipping"
    else
      echo "Working on: "$RefID

      check_info=$(grep -w -F -c $RefID $annotation_keystone)
      if [ $check_info -eq 0 ]
      then
        echo $RefID" is not on the instruction file:"
        echo $annotation_keystone
        exit
      else
        prefix_extract=$(grep -w -F $RefID $annotation_keystone | awk -F "\t" '{print $2}')
        subfix_extract=$(grep -w -F $RefID $annotation_keystone | awk -F "\t" '{print $3}')
        take_item=$(grep -w -F $RefID $annotation_keystone | awk -F "\t" '{print $4}')
      fi

      gff_file=$(awk -F "\t" -v genome=$RefID '{if (genome==$1) print $4}' $genome_table )

      echo "GeneID ProteinID GO_TERM" | tr " " "\t" >> $work_folder"/"$storage"/Work_version/Eggnog/"$RefID"_goterms_per_isoform.tsv"
      count=6
      walker=$(grep -c . $work_folder"/"$storage"/Eggnog_mapper/"$RefID".emapper.annotations")
      while [ $count -le $walker ]
      do
        protID=$(sed -n $count"p" $work_folder"/"$storage"/Eggnog_mapper/"$RefID".emapper.annotations" | awk -F "\t" '{print $1}')
        check_end_of_file=$(echo $protID | grep -c -F "##")
        if [ $check_end_of_file -eq 0 ]
        then
          goterms=$(sed -n $count"p" $work_folder"/"$storage"/Eggnog_mapper/"$RefID".emapper.annotations" | awk -F "\t" '{print $10}')
          gene_loci=$(grep -F $prefix_extract$protID$subfix_extract $gff_file | awk -F "\t" '{print $9}' | head -n 1 | tr ";" "\n" | grep -F $take_item | sed "s/$take_item//")

          echo $gene_loci" "$protID" "$goterms | tr " " "\t" >> $work_folder"/"$storage"/Work_version/Eggnog/"$RefID"_goterms_per_isoform.tsv"
          count=$(($count + 1))
        else
          count=$(($walker + $walker))
        fi
      done

      echo $RefID >> $work_folder"/"$storage"/Work_version/Already_done_eggnog.txt"
    fi
  done
  code_block_completed
  slurm_block_completed
  echo ""
  echo "Note: This block is a bit difference. It generates the scripts that are going to be run in SLURM. Before running them you need to give them execution permits"
fi

################################################################################

if [ $SLURM_GENE_SEQ_EGGNOG_Re == TRUE ]
then
  storage=Ref_Genes_Annotation

  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Gene_Annotation_"$current_date)

  confirm_input_file
  start_block_code

  if [ ! -d $work_folder"/"$storage"/Eggnog_mapper" ]
  then
    mkdir $work_folder"/"$storage"/Eggnog_mapper"
  fi

  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $genome_table | sort)
  for RefID in $reference_assemblies
  do
    echo "Working on: "$RefID
    check_already_done=$(grep -w -F -c $RefID"_emapper" $work_folder"/"$storage"/Already_done.txt")

    if [ $check_already_done -gt 0 ]
    then
      echo "Already done"
      echo "Skipping..."
    else
      prot_file=$(awk -F "\t" -v genome=$RefID '{if (genome==$1) print $5}' $genome_table)
      if [ $prot_file != "NA" ]
      then
        echo "emapper.py --cpu "'$NPROCS'" -i "$prot_file" -o "$work_folder"/"$storage"/Eggnog_mapper/"$RefID >> $cluster_submission_folder"/"$submit_to_cluster"_Eggnog.txt"
      else
        echo "No Protein file for this reference. Skipping"
      fi
      echo $RefID"_emapper" >> $work_folder"/"$storage"/Already_done.txt"
    fi
  done
  code_block_completed
  slurm_block_completed
fi

################################################################################

if [ $INTERPRO_PER_GENE_INDEX == TRUE ]
then
  storage=Ref_Genes_Annotation
  confirm_input_file

  if [ -d $work_folder"/"$storage"/Work_version/Interpro/"$RefID".index" ]
  then
    rm -rf $work_folder"/"$storage"/Work_version/Interpro/split_temp"
  fi

  annotated_assemblies=$(ls $work_folder"/"$storage"/Interpro_scan/"*".tsv" | sed "s/\/$//" | sed "s/.*\///" | sed "s/.tsv//")
  for RefID in $annotated_assemblies
  do
    if [ -f $work_folder"/"$storage"/Work_version/Interpro/"$RefID".index" ]
    then
      rm -f $work_folder"/"$storage"/Work_version/Interpro/"$RefID".index"
    fi

    split_files=$(ls $work_folder"/"$storage"/Work_version/Interpro/"$RefID"_work_table_"*".tsv" | sed "s/\/$//" | sed "s/.*\///")
    for split in $split_files
    do
      echo "Indexing: "$RefID" -- "$split
      awk -F "\t" '{print $1}' $work_folder"/"$storage"/Work_version/Interpro/"$split | sed 1d | sort -u | awk -F "\t" -v file=$split '{print $1"\t"file}' >> $work_folder"/"$storage"/Work_version/Interpro/"$RefID".index"
    done
  done
  echo "Done with the index"
fi

################################################################################


if [ $SLURM_SIGNALP6_Re == TRUE ]
then
  confirm_input_file
  storage=Ref_Genes_Annotation

  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Gene_Annotation_"$current_date)

  if [ ! -d $work_folder"/"$storage"/SignalP6" ]
  then
    mkdir $work_folder"/"$storage"/SignalP6"
    mkdir $work_folder"/"$storage"/SignalP6/Temp_sequences"
    echo "Already Done" >> $work_folder"/"$storage"/SignalP6/Already_done.txt"
  fi

  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $genome_table | sort)
  for RefID in $reference_assemblies
  do
    echo "Working on: "$RefID
    check_already_done=$(grep -w -F -c $RefID $work_folder"/"$storage"/SignalP6/Already_done.txt")

    if [ $check_already_done -gt 0 ]
    then
      echo "Already done"
      echo "Skipping..."
    else
      if [ -d $work_folder"/"$storage"/SignalP6/"$RefID ]
      then
        rm -rf $work_folder"/"$storage"/SignalP6/"$RefID
      fi

      prot_file=$(awk -F "\t" -v genome=$RefID '{if (genome==$1) print $5}' $genome_table)
      seqkit seq -i $prot_file >> $work_folder"/"$storage"/SignalP6/Temp_sequences/"$RefID".aa"

      if [ $prot_file != "NA" ]
      then
        echo "signalp6-cpu --fastafile "$work_folder"/"$storage"/SignalP6/Temp_sequences/"$RefID".aa --organism eukarya --output_dir "$work_folder"/"$storage"/SignalP6/"$RefID >> $cluster_submission_folder"/"$submit_to_cluster"_signalp6.txt"
      else
        echo "No Protein file for this reference. Skipping"
      fi
      echo $RefID >> $work_folder"/"$storage"/SignalP6/Already_done.txt"
    fi
  done
  code_block_completed
  slurm_block_completed
  echo ""
  echo "Known issues: "
  echo "For unknown reasons I cannot run multple threads on SLURM with signalp6-cpu."
  echo "Temporary solution: Just don't. It will take longer but it will complete"
fi

################################################################################

if [ $SLURM_TMHMM_Re == TRUE ]
then
  confirm_input_file
  storage=Ref_Genes_Annotation

  current_date=$(date | tr " " "_" | tr -s "_")
  submit_to_cluster=$(echo "Gene_Annotation_"$current_date)

  if [ ! -d $work_folder"/"$storage"/TMHMM" ]
  then
    mkdir $work_folder"/"$storage"/TMHMM"
    echo "Already Done" >> $work_folder"/"$storage"/TMHMM/Already_done.txt"
  fi

  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $genome_table | sort)
  for RefID in $reference_assemblies
  do
    echo "Working on: "$RefID
    check_already_done=$(grep -w -F -c $RefID $work_folder"/"$storage"/TMHMM/Already_done.txt")
    prot_file=$(awk -F "\t" -v genome=$RefID '{if (genome==$1) print $5}' $genome_table)

    if [ $prot_file == "NA" ]
    then
      echo "No Protein file for this reference. Skipping"
    elif [ ! -f $work_folder"/"$storage"/SignalP6/Temp_sequences/"$RefID".aa" ]
    then
      echo "Can't find the sequence file generated for signalP run: "$work_folder"/"$storage"/SignalP6/Temp_sequences/"$RefID".aa"
      echo "run that module first"
      exit
    elif [ $check_already_done -gt 0 ]
    then
      echo "Already done"
      echo "Skipping..."
    else
      if [ -d $work_folder"/"$storage"/TMHMM/"$RefID ]
      then
        rm -rf $work_folder"/"$storage"/TMHMM/"$RefID
      fi

      echo "tmhmm "$work_folder"/"$storage"/SignalP6/Temp_sequences/"$RefID".aa >> "$work_folder"/"$storage"/TMHMM/"$RefID"_tmhmm.txt" >> $cluster_submission_folder"/"$submit_to_cluster"_tmhmm.txt"
      echo $RefID >> $work_folder"/"$storage"/TMHMM/Already_done.txt"
    fi
  done

  code_block_completed
  slurm_block_completed
  echo ""
  echo "Known issues: "
  echo "tmhmm will generate additional outputs in the working directory and there is not a way to configure where. At least not a simple one"
  echo "after the script is done you will have to work with them as you see fit"
fi

################################################################################

if [ $LOCAL_TARGETP == TRUE ]
then
  confirm_input_file
  storage=Ref_Genes_Annotation

  if [ -d $work_folder"/"$storage"/TargetP" ]
  then
    rm -rf -d $work_folder"/"$storage"/TargetP"
  fi

  mkdir $work_folder"/"$storage"/TargetP"
  mkdir $work_folder"/"$storage"/TargetP/Sequnces"
  mkdir $work_folder"/"$storage"/TargetP/Output_External"
  echo "Already Done" >> $work_folder"/"$storage"/TargetP/Already_done.txt"
  echo "Please downladu this folder locally and upload each part to this server https://services.healthtech.dtu.dk/services/TargetP-2.0" >> $work_folder"/"$storage"/TargetP/Instructions.txt"
  echo "Store the summary table as [RefID].part_[num].txt" >> $work_folder"/"$storage"/TargetP/Instructions.txt"
  echo "The script will only use that table for results. Anything else will require to re-run particular sequences on the server." >> $work_folder"/"$storage"/TargetP/Instructions.txt"
  echo "Reason: outputs are messy and the 'Download everything button' does not work at time of writing this." >> $work_folder"/"$storage"/TargetP/Instructions.txt"

  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $genome_table | sort)
  for RefID in $reference_assemblies
  do
    echo "Working on: "$RefID
    mkdir $work_folder"/"$storage"/TargetP/Sequnces/"$RefID

    prot_file=$(awk -F "\t" -v genome=$RefID '{if (genome==$1) print $5}' $genome_table)
    if [ $prot_file != "NA" ]
    then
      cp $prot_file $work_folder"/"$storage"/TargetP/Sequnces/"$RefID".aa"
      seqkit split -s 5000 $work_folder"/"$storage"/TargetP/Sequnces/"$RefID".aa" -O $work_folder"/"$storage"/TargetP/Sequnces/"$RefID
      rm -f $work_folder"/"$storage"/TargetP/Sequnces/"$RefID".aa"
    else
      echo "No Protein file for this reference. Skipping"
    fi
    echo $RefID >> $work_folder"/"$storage"/SignalP6/Already_done.txt"
  done
fi

################################################################################

if [ $REFORMAT_TARGETP_OUTPUT == TRUE ]
then
  confirm_input_file
  storage=Ref_Genes_Annotation

  reference_assemblies=$(awk -F "\t" '{if ($2=="R") print $1}' $genome_table | sort)
  for RefID in $reference_assemblies
  do
    echo "Working on: "$RefID
    if [ -f $work_folder"/"$storage"/TargetP/"$RefID"_TargetP.txt" ]
    then
      rm -f $work_folder"/"$storage"/TargetP/"$RefID"_TargetP.txt"
    fi

    echo "GeneID guidestone"
    gff_file=$(awk -F "\t" -v genome=$RefID '{if (genome==$1) print $4}' $genome_table)
    awk -F "\t" '{if ($3!="exon") print $9}' $gff_file | sort -u > $work_folder"/"$storage"/TargetP/current_genes.tmp"
    echo "done"

    # extract gene info
    current_ref=$RefID
    get_prefix_info_for_gff

    targetp_parts=$(ls $work_folder"/"$storage"/TargetP/Output_External/"$RefID"."*"txt")
    echo "GeneID PortID Prediction OTHER SP mTP CS_Position" | tr " " "\t"  > $work_folder"/"$storage"/TargetP/"$RefID"_TargetP.txt"

    for part in $targetp_parts
    do
      echo "current: "$part
      echo ""
      count=3
      walker=$(grep -c . $part)
      while [ $count -le $walker ]
      do
        protID=$(sed -n $count"p" $part | awk -F "\\t" '{print $1}')
        gene_loci=$(grep -F $prefix_extract$protID$subfix_extract $work_folder"/"$storage"/TargetP/current_genes.tmp" | tr ";" "\n" | grep -F $take_item | sed "s/$take_item//" | sort -u)
        echo $gene_loci": "$protID
        sed -n $count"p" $part | awk -F "\t" -v loci=$gene_loci '{if ($2=="OTHER") print loci"\t"$0"NA";else print loci"\t"$0}' >> $work_folder"/"$storage"/TargetP/"$RefID"_TargetP.txt"

        count=$(($count + 1))
      done
    done
    rm -f $work_folder"/"$storage"/TargetP/current_genes.tmp"
  done
fi

################################################################################

if [ $REGION_GENE_EXTRACTION_Mk2 == TRUE ]
then
  confirm_input_file
  confirm_region_file
  confirm_final_region_file

  storage=Definitive_Regions
  scaffolded_assemblies=UNMASKED_SCAFFOLDING
  liftoff_annotation=Annotation_Transfer_Liftoff
  syn_block_inspection=Raw_Syn_Block_Inspection_Mk2
  query_group_ragtag=RAGTAG_Scaffolding
  reference_group_ragtag=RAGTAG_Scaffolding_References

  start_block_code
  rm -f $work_folder"/"$storage"/Already_done.txt"

  all_regions=$(sed 1d  $final_region_recovery | awk -F "\t" '{print $1}' | sort -u)
  for region_ID in $all_regions
  do
    if [ -d $work_folder"/"$storage"/"$region_ID ]
    then
      rm -rf $work_folder"/"$storage"/"$region_ID
    fi

    mkdir $work_folder"/"$storage"/"$region_ID
    mkdir $work_folder"/"$storage"/"$region_ID"/Other_Scaffolded_checks"
    mkdir $work_folder"/"$storage"/"$region_ID"/Other_Scaffolded_checks/Sequences"
    mkdir $work_folder"/"$storage"/"$region_ID"/Other_Scaffolded_checks/Contig_Order"
    mkdir $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions"
    mkdir $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Sequences"
  done

  echo "AssemblyID Annotiation_Donor Is_Original Is_scaffolded Annotation_file" | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Region_work_summary.txt"

  count_final_regions=2
  while [ $count_final_regions -le $total_regions ]
  do
    region_ID=$(sed -n $count_final_regions"p" $final_region_recovery | awk -F "\t" '{print $1}')

    assembly_ID=$(sed -n $count_final_regions"p" $final_region_recovery | awk -F "\t" '{print $2}')
    assembly_group=$(awk -F "\t" -v genome=$assembly_ID '{if (genome==$1) print $2}' $genome_table)

    is_scaffolded=$(sed -n $count_final_regions"p" $final_region_recovery | awk -F "\t" '{print $5}')
    scaffold_reference=$(sed -n $count_final_regions"p" $final_region_recovery | awk -F "\t" '{print $6}')

    seq_file_path=$(sed -n $count_final_regions"p" $final_region_recovery | awk -F "\t" '{print $3}')
    scaffold_file=$(echo $seq_file_path | sed "s/.*\///" )

    location_info=$(sed -n $count_final_regions"p" $final_region_recovery | awk -F "\t" '{print $4}' | tr ";" " ")

    if [ $is_scaffolded == "TRUE" ]
    then
      if [ $assembly_group == "Q" ]
      then
        ragtag_output=$work_folder"/"$query_group_ragtag"/Ref_"$original_descriptor"/"$assembly_ID"_ragtag_output/ragtag.scaffold.agp"
      else
        ragtag_folder=$(echo $seq_file_path| sed "s/.*\///" | sed "s/_unmasked_scaffolds.fasta//" | sed "s/_vs_/_vs_Ref_/" | awk '{print "Ref_"$1"_ragtag_output"}')
        ragtag_output=$work_folder"/"$reference_group_ragtag"/"$ragtag_folder"/ragtag.scaffold.agp"
      fi
    fi

    if [ ! -f $work_folder"/"$storage"/"$region_ID"/Other_Scaffolded_checks/Contig_Order/"$assembly_ID"_original_contig_order.txt" ]
    then
      echo "Block_Num Scaffolded_ID Contig_ID Coods Strand" | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Other_Scaffolded_checks/Contig_Order/"$assembly_ID"_original_contig_order.txt"
    fi

    echo "Working on: "$region_ID" -- "$assembly_ID" ("$is_scaffolded")"

    block=1
    echo $assembly_ID" Extracting sequences..."
    for target in $location_info
    do
      chr=$(echo $target | awk -F ":" '{print $1}')
      coords=$(echo $target | awk -F ":" '{print $2}' | sed "s/.$//" | tr "-" ":")
      strand=$(echo $target | awk -F ":" '{print $2}' | sed "s/./\n&/g" | tail -n 1 ) # if it is stupid but it works...
      header=$(echo ">"$region_ID"_"$assembly_ID"_block-"$block)

      if [ -f $seq_file_path".seqkit.fai" ]
      then
        rm -f $seq_file_path".seqkit.fai"
      fi

      echo "Block "$block": "$chr" | "$coords
      seqkit subseq --quiet -j $threads --chr $chr -r $coords $seq_file_path | sed "s/>.*/$header/" >> $work_folder"/"$storage"/"$region_ID"/Other_Scaffolded_checks/Sequences/"$assembly_ID"_region_blocks.fasta"
      echo "Done Block Extraction"

      if [ $is_scaffolded == "FALSE" ]
      then
        echo $block" NA "$chr" "$coords" "$strand | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Other_Scaffolded_checks/Contig_Order/"$assembly_ID"_original_contig_order.txt"
      else
        start_coord_ragtag=$(echo $coords | awk -F ":" '{print $1}')
        end_coord_ragtag=$(echo $coords | awk -F ":" '{print $2}')

        if [ $strand == "+" ]
        then
          awk -F "\t" -v chr=$chr -v end=$end_coord_ragtag '{if ($1==chr && $5!="U" && $5!="N" && $2 <= end) print}' $ragtag_output | awk -F "\t" -v start=$start_coord_ragtag -v block=$block '{if (start <= $3) print block"\t"$1"\t"$6"\t"$2":"$3"\t"$9}' >> $work_folder"/"$storage"/"$region_ID"/Other_Scaffolded_checks/Contig_Order/"$assembly_ID"_original_contig_order.txt"
        else
          awk -F "\t" -v chr=$chr -v end=$end_coord_ragtag '{if ($1==chr && $5!="U" && $5!="N" && $2 <= end) print}' $ragtag_output | sort -r -n -k 2 | awk -F "\t" -v start=$start_coord_ragtag -v block=$block '{if (start <= $3) print block"\t"$1"\t"$6"\t"$2":"$3"\t"$9}' | awk -F "\t" '{if ($5=="+") print $1"\t"$2"\t"$3"\t"$4"\t-"; else print $1"\t"$2"\t"$3"\t"$4"\t+"}'  >> $work_folder"/"$storage"/"$region_ID"/Other_Scaffolded_checks/Contig_Order/"$assembly_ID"_original_contig_order.txt"
        fi
      fi

      block=$(($block + 1))
    done

    echo "Working on annotations..."
    if [ $is_scaffolded == "TRUE" ]
    then
      check_liftoff=$(echo $scaffold_reference"_vs_"$assembly_ID)

      liftoff_annotation_transfers=$(ls $work_folder"/"$liftoff_annotation"/Individual_Transfers/"*"_transfered_to_"$check_liftoff"_annotation.gff_polished" | sed "s/.*\///")
      for liftoff in $liftoff_annotation_transfers
      do
        annot_file=$work_folder"/"$liftoff_annotation"/Individual_Transfers/"$liftoff
        donor=$(echo $liftoff | sed "s/_transfered_to_/\n/" | head -n 1)
        original=FALSE
        fraking_gene_things
      done
    else
      liftoff_annotation_transfers=$(ls $work_folder"/"$liftoff_annotation"/Individual_Transfers/"*"_transfered_to_"$assembly_ID"_annotation.gff_polished" | sed "s/.*\///")
      for liftoff in $liftoff_annotation_transfers
      do
        annot_file=$work_folder"/"$liftoff_annotation"/Individual_Transfers/"$liftoff
        donor=$(echo $liftoff | sed "s/_transfered_to_/\n/" | head -n 1)
        original=FALSE
        fraking_gene_things
      done
      annot_file=$(awk -F "\t" -v genome=$assembly_ID '{if (genome==$1) print $4}' $genome_table )
      donor=$assembly_ID
      original=TRUE
      fraking_gene_things
    fi
    echo "Next!!"
    echo ""
    count_final_regions=$(($count_final_regions + 1))
  done
  code_block_completed
  echo ""
  echo "Known issues:"
  echo "Unfortunately due to time and the inner approach I took for this block, I wasn't able to make it aware of previous runs"
  echo "You can still work with new regions without removing previos work by changing the file given on the variable final_region_recovery"
fi

################################################################################

if [ $REGION_BASIC_SUMMARY == TRUE ]
then
  confirm_input_file
  confirm_region_file
  confirm_final_region_file
  storage=Definitive_Regions

  all_regions=$(sed "1d" $final_region_recovery | awk -F "\t" '{print $1}' | sort -u)
  for region_ID in $all_regions
  do
    all_samples=$(awk -v region=$region_ID '{if ($1==region) print $2}' $final_region_recovery | sort -u)
    all_annotation_donors=$(ls $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"* | sed "s/.*Source_//" | sed "s/_transfered_to_.*//" | sort -u)
    if [ -f $work_folder"/"$storage"/"$region_ID"/Sequence_Block_summary.txt" ]
    then
      rm -f $work_folder"/"$storage"/"$region_ID"/Sequence_Block_summary.txt"
    fi

    header=$(echo $all_annotation_donors | tr " " "\n" | awk '{print "N_Genes_"$1}')
    echo "AssemblyID Block Chr Start End Lenght "$header" Scaffolding_Reference" | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Sequence_Block_summary.txt"
    for sample in $all_samples
    do
      blocks=$(grep -v -w Gene_On_Target $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"*"_transfered_to_"$sample"_info.txt" | awk -F "\t" '{print $13}' | uniq -c | awk '{print $2}' | sort -u)
      for bl in $blocks
      do
        chr=$(cat $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"*"_transfered_to_"$sample"_info.txt" | awk -F "\t" -v block=$bl '{if ($13==block) print $3}' | sort -u)
        start=$(cat $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"*"_transfered_to_"$sample"_info.txt" | awk -F "\t" -v block=$bl '{if ($13==block) print $4}' | sort -n | head -n 1 )
        end=$(cat $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"*"_transfered_to_"$sample"_info.txt" | awk -F "\t" -v block=$bl '{if ($13==block) print $5}' | sort -n | tail -n 1 )
        lenght=$(($end - $start))
        scaffold_reference=$(cat $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"*"_transfered_to_"$sample"_info.txt" | awk -F "\t" -v block=$bl '{if ($13==block) print $14}' | sort -u)

        print_information=$(echo $sample" "$bl" "$chr" "$start" "$end" "$lenght )
        for donor in $all_annotation_donors
        do
          N_genes=$(awk -F "\t" -v block=$bl '{if ($13==block) print}' $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"$sample"_info.txt" | grep -c .)
          echo $donor": "$N_genes
          print_information=$(echo $print_information" "$N_genes)
        done

        echo $print_information" "$scaffold_reference | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Sequence_Block_summary.txt"
      done
    done
  done
fi

################################################################################

if [ $REGION_GENE_ANOTATION_INTERPRO_Re == TRUE ]
then
  confirm_input_file
  confirm_region_file
  confirm_final_region_file

  storage=Definitive_Regions
  interpro_annotation=Ref_Genes_Annotation/Work_version/Interpro

  all_regions=$(sed "1d" $final_region_recovery | awk -F "\t" '{print $1}' | sort -u)
  for region_ID in $all_regions
  do
    all_samples=$(awk -v region=$region_ID '{if ($1==region) print $2}' $final_region_recovery | sort -u)
    all_annotation_donors=$(ls $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"* | sed "s/.*Source_//" | sed "s/_transfered_to_.*//" | sort -u)

    if [ ! -d $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation" ]
    then
      mkdir $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation"
    fi

    if [ ! -d $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro" ]
    then
      mkdir $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro"
      echo "Already done" >>  $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/Already_done.txt"
    fi

    for donor in $all_annotation_donors
    do
      echo "LociID N_Assemblies Analysis SignatureID Signature_Description InterproID Interpro_Description" | tr " " "\t" > $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/Gene_annotation_"$donor".txt"

      if [ ! -f  $work_folder"/"$interpro_annotation"/"$donor".index" ]
      then
        echo "Sojourn Audio Drama! Listen to it!!"
        echo "Also the index for "$donor" annotation is not there. Please run INTERPRO_PER_GENE_INDEX"
        exit
      fi

      all_genes=$(awk -F "\t" '{print $1}' $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"* | grep -w -v -F "GeneID" | sort -u )
      echo "Heeeyyyy... What if the gene is a duplicate found by Liftoff?" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/duplicate.tmp"
      grep -v -w "GeneID"  $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"*  | awk -F "\t" '{if ($11>=1) print}' | grep -v -w Original | awk -F "\t" '{print $1}' | awk -F ":" '{print $2}' | grep . | sort -u >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/duplicate.tmp"

      for genid in $all_genes
      do
        total_assemblies=$(grep -w -F $genid $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"* | awk -F ":" '{print $1}' | sort -u | grep -c .)
        check_duplicate=$(grep -c -w -F $genid $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/duplicate.tmp")
        if [ $check_duplicate -gt 0 ]
        then
          search_genid=$(echo $genid | sed "s/^gene-//")
          remove_copy_num=$(echo $search_genid | sed "s/./\n&/g" | grep . | grep -n "_" | tail -n 1 | awk -F ":" '{print $1-1}')
          search_genid=$(echo $search_genid | sed "s/./\n&/g" | grep . | head -n $remove_copy_num | tr -d "\n" | sed "s/$/\n/")
        else
          search_genid=$(echo $genid | sed "s/^gene-//")
        fi

        N_Annotation_Files=$(awk -F "\t" -v gen=$search_genid '{if ($1==gen) print $2}'  $work_folder"/"$interpro_annotation"/"$donor".index" | sort -u | grep -c .)
        echo $genid" -- "$search_genid" -- "$total_assemblies" -- "$N_Annotation_Files

        if [ $N_Annotation_Files -gt 0 ]
        then
          annotation_files=$(awk -F "\t" -v gen=$search_genid '{if ($1==gen) print $2}'  $work_folder"/"$interpro_annotation"/"$donor".index" | sort -u )
          for annot_part in $annotation_files
          do
            awk -F "\t" -v gen=$search_genid -v print_gen=$genid -v N_Assemblies=$total_assemblies '{if ($1==gen) print print_gen"\t"N_Assemblies"\t"$5"\t"$6"\t"$7"\t"$13"\t"$14}' $work_folder"/"$interpro_annotation"/"$annot_part  >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/annot.tmp"
          done
          sort -u $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/annot.tmp" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/Gene_annotation_"$donor".txt"
          rm -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/annot.tmp"
        else
          echo $genid" "$total_assemblies" NA NA NA NA NA" | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/Gene_annotation_"$donor".txt"
        fi
      done
    done
  done
fi

################################################################################

if [ $REGION_GENE_ANOTATION_SIGNALP_TARGETP_TMHMM == TRUE ]
then
  confirm_input_file
  confirm_region_file
  confirm_final_region_file

  storage=Definitive_Regions
  signalP_annotation=Ref_Genes_Annotation/SignalP6
  targetP_annotation=Ref_Genes_Annotation/TargetP
  tmhmm_annotation=Ref_Genes_Annotation/TMHMM

  all_regions=$(sed "1d" $final_region_recovery | awk -F "\t" '{print $1}' | sort -u)
  for region_ID in $all_regions
  do
    all_samples=$(awk -v region=$region_ID '{if ($1==region) print $2}' $final_region_recovery | sort -u)
    all_annotation_donors=$(ls $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"* | sed "s/.*Source_//" | sed "s/_transfered_to_.*//" | sort -u)

    if [ ! -d $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation" ]
    then
      mkdir $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation"
    fi

    if [ ! -d $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals" ]
    then
      mkdir $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals"
    fi

    for donor in $all_annotation_donors
    do
      echo "Heeeyyyy... What if the gene is a duplicate found by Liftoff?" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/duplicate.tmp"
      echo "GeneID PortID Prediction OTHER SP mTP CS_Position" | tr " " "\t"  > $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/"$donor"_TargetP.txt"
      echo "GeneID PortID Prediction OTHER SP(Sec/SPI) CS_Position" | tr " " "\t" > $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/"$donor"_SignalP6.txt"
      echo "GeneID PortID N_TMhelix" | tr " " "\t" > $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/"$donor"_tmhmm.txt"

      grep -v -w "GeneID"  $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"*  | awk -F "\t" '{if ($11>=1) print}' | grep -v -w Original | awk -F "\t" '{print $1}' | awk -F ":" '{print $2}' | grep . | sort -u > $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/duplicate.tmp"
      grep -v -w "GeneID"  $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"*  | awk -F "\t" '{if ($11==0) print $1}' | awk -F "\t" '{print $1}' | awk -F ":" '{print $2}' | grep . | sort -u | sed "s/gene-//" > $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/genid.tmp"

      all_genes=$(cat $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/duplicate.tmp")
      for genid in $all_genes
      do
        echo "Temp file make: "$genid
        total_assemblies=$(grep -w -F $genid $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"* | awk -F ":" '{print $1}' | sort -u | grep -c .)
        check_duplicate=$(grep -c -w -F $genid $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/duplicate.tmp")
        if [ $check_duplicate -gt 0 ]
        then
          search_genid=$(echo $genid | sed "s/^gene-//")
          remove_copy_num=$(echo $search_genid | sed "s/./\n&/g" | grep . | grep -n "_" | tail -n 1 | awk -F ":" '{print $1-1}')
          search_genid=$(echo $search_genid | sed "s/./\n&/g" | grep . | head -n $remove_copy_num | tr -d "\n" | sed "s/$/\n/")
        else
          search_genid=$(echo $genid | sed "s/^gene-//")
        fi
        echo $search_genid >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/genid.tmp"
      done

      work_search_geneid=$(sort -u $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/genid.tmp")
      for search_genid in $work_search_geneid
      do
        # this one will require
        check_target_p_data=$(awk -F "\t" -v gen=$search_genid '{if ($1==gen) print}' $work_folder"/"$targetP_annotation"/"$donor"_TargetP.txt" | grep -c . )
        if [ $check_target_p_data -gt 0 ]
        then
          awk -F "\t" -v gen=$search_genid '{if ($1==gen) print}' $work_folder"/"$targetP_annotation"/"$donor"_TargetP.txt" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/"$donor"_TargetP.txt"
          prot_IDs=$(awk -F "\t" -v gen=$search_genid '{if ($1==gen) print $2}' $work_folder"/"$targetP_annotation"/"$donor"_TargetP.txt" | sort -u)
          for prot in $prot_IDs
          do
            echo "Checking: "$donor" -- "$search_genid" ("$prot")"
            awk -F "\t" -v gen=$search_genid -v prot=$prot '{if ($1==prot) print gen"\t"$0}' $work_folder"/"$signalP_annotation"/"$donor"/prediction_results.txt" | sed "s/\t$/\tNA/" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/"$donor"_SignalP6.txt"
            grep "# "$prot" Number of predicted TMHs:" $work_folder"/"$tmhmm_annotation"/"$donor"_tmhmm.txt" | awk -F "  " -v gen=$search_genid -v prot=$prot '{print gen"\t"prot"\t"$2}' >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/"$donor"_tmhmm.txt"
          done
        fi
      done
      rm -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/duplicate.tmp"
      rm -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/genid.tmp"
    done
  done

  echo "Known issues:"
  echo "Because this is an add on and I need to actually start writting papers and justify my existence: this block of code relies on TargetP results being provided and reformated, so that that file has the protein IDs for every gene."
  echo "Sorry for that."
fi

################################################################################

if [ $REGION_GENE_ANOTATION_TRACKS == TRUE ]
then
  confirm_input_file
  confirm_region_file
  confirm_final_region_file

  storage=Definitive_Regions
  scaffolded_assemblies=UNMASKED_SCAFFOLDING
  liftoff_annotation=Annotation_Transfer_Liftoff
  syn_block_inspection=Raw_Syn_Block_Inspection_Mk2

  all_regions=$(sed "1d" $final_region_recovery | awk -F "\t" '{print $1}' | sort -u)
  for region_ID in $all_regions
  do
    if [ -d $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks" ]
    then
      rm -rf $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks"
    fi
    mkdir $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks"

    annotation_donors=$(ls $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/Gene_annotation_"*".txt" | sed "s/.*\///" | sed "s/Gene_annotation_//" | sed "s/.txt//")

    all_assembly=$(sed "1d" $final_region_recovery | awk -F "\t" '{print $2}' | sort -u)

    for assembly_ID in $all_assembly
    do
      mkdir $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$assembly_ID
      mkdir $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$assembly_ID"/Annotation_tracks"

      for donor in $annotation_donors
      do
        echo "Working on: "$region_ID" -- "$assembly_ID" -- "$donor

        take_item=$(grep -w -F $donor $annotation_keystone | awk -F "\t" '{print $4}')
        awk -F "\t" '{print $1}' $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"$assembly_ID"_info.txt" | sed 1d | sed "s/^gene-//" | awk -v take=$take_item '{print take$1}' >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/assembly_genes.tmp"

        duplicated_genes=$(sed 1d $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"$assembly_ID"_info.txt" | sed "s/Original/0/g"| awk -F "\t" '{if ($11 > 0) print $1}')
        run_duplicated=FALSE
        for dup in $duplicated_genes
        do
          run_duplicated=TRUE
          search_genid=$(echo $dup | sed "s/^gene-//")
          remove_copy_num=$(echo $search_genid | sed "s/./\n&/g" | grep . | grep -n "_" | tail -n 1 | awk -F ":" '{print $1-1}')
          search_genid=$(echo $search_genid | sed "s/./\n&/g" | grep . | head -n $remove_copy_num | tr -d "\n" | sed "s/$/\n/")

          echo $dup" "$search_genid | tr " " "\t" >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/duplicated.tmp"
        done

        annotation_file=$(awk -F "\t" -v assembly=$assembly_ID -v donor=$donor '{if ($1==assembly && $2==donor) print $5}' $work_folder"/"$storage"/"$region_ID"/Region_work_summary.txt")

        current_annotation_file=$annotation_file
        wanted_genes=$work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/assembly_genes.tmp"
        out_track=$work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$assembly_ID"/All_genes_from_"$donor".gff"

        print_current=All_genes
        detailed_gene_track

        current_annotation_file=$out_track

        count_domain=2
        walker_domain=$(grep -c . $interpro_domain_combo)

        sed 1d $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Location_Signals/"$donor"_tmhmm.txt" | awk -F "\t" '{if ($3>0) print $1}' | sort -u >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/Transmembrane_genes.tmp"
        wanted_genes=$work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/Transmembrane_genes.tmp"
        out_track=$work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$assembly_ID"/Annotation_tracks/Transmembrane_genes_from_"$donor".gff"

        print_current=TRANSMEMBRANE
        detailed_gene_track

        awk -F "\t" '{if ($3=="Coils") print $1}' $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/Gene_annotation_"$donor".txt" | sed "s/^gene-//" | awk -v take=$take_item '{print take$1}' >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/Coils_genes.tmp"

        wanted_genes=$work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/Coils_genes.tmp"
        out_track=$work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$assembly_ID"/Annotation_tracks/Coils_genes_from_"$donor".gff"
        print_current=Coils
        detailed_gene_track

        awk -F "\t" '{if ($3=="MobiDBLite") print $1}' $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/Gene_annotation_"$donor".txt" | sed "s/^gene-//" | awk -v take=$take_item '{print take$1}' >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/Disorder_genes.tmp"

        wanted_genes=$work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/Disorder_genes.tmp"
        out_track=$work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$assembly_ID"/Annotation_tracks/Disorder_genes_from_"$donor".gff"
        print_current=Disorder
        detailed_gene_track


        while [ $count_domain -le $walker_domain ]
        do
          identifier_ID=$(sed -n $count_domain"p" $interpro_domain_combo | awk -F "\t" '{print $1}')
          sed -n $count_domain"p" $interpro_domain_combo | awk -F "\t" '{print $2}' | tr ";" "\n" | awk -F "\t" '{print "_ooo_"$1"_ooo_"}' >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$identifier_ID".tmp"

          sed "s/\t/_ooo_/g" $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation/Interpro/Gene_annotation_"$donor".txt" | grep -F -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$identifier_ID".tmp" | awk -F "_ooo_" '{print $1}' >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$identifier_ID"_matching_genes.tmp"
          awk -F "\t" '{print $1}' $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Source_"$donor"_transfered_to_"$assembly_ID"_info.txt" | grep -w -F -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$identifier_ID"_matching_genes.tmp" | sed "s/^gene-//" | awk -v take=$take_item '{print take$1}' >> $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$identifier_ID"_assembly_genes.tmp"

          print_current=$identifier_ID
          wanted_genes=$work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$identifier_ID"_assembly_genes.tmp"
          out_track=$work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$assembly_ID"/Annotation_tracks/"$identifier_ID"_from_"$donor".gff"
          detailed_gene_track

          rm -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$identifier_ID"_matching_genes.tmp"
          rm -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/"$identifier_ID".tmp"

          count_domain=$(($count_domain + 1))
        done

        if [ -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/duplicated.tmp" ]
        then
          rm -f $work_folder"/"$storage"/"$region_ID"/Genes_on_Regions/Annotation_tracks/duplicated.tmp"
        fi
      done
    done
  done

  echo "Code Block complete!!"
  echo "Known issue: Liftoff doesn't updates all of the information on the gff and to extract the annotated features I'm using either the assigned gene from the original GFF, or the locus tag."
  echo "I haven't managed to make gffread to work as intended to extract based on gene IDs, so here is the trade off: "
  echo "Genes identified as extra coppies from the original annotation come as a package. if there are 40 copies the gff tracks will have all 40 copies, regardless of where they are. Since the annotations of all copies are assumed to be the same, it doesn't change anything beyond hurting my pride."
fi
