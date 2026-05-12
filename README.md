## Identification and re-annotation of regions associated with resistance against Schistosoma mansoni in Biomphalaria glabrata
This repository contains the in-house scripts utilized for the publication [XXXXXX]
The goal was to identify, locate and re-annotate across multiple genomes of the species *Biomphalaria glabrata* 2 new regions associated with resistance against *Schistosoma mansoni*. While designed with re-usability in mind, adapting the scrips to work in a different set up might prove challenging as many design desissions respond to local system architecture. For example many analysis are not run directly but instead the script generates input files to be uploaded into a cluster to run in paralel. 

We will provide assistance upon request.

## BASH scripts
### Main_Script.sh
In short, first the region of interest must be defined in one reference assembly defined by the user. Then all assemblies are scaffolded relative to each of the reference genomes using RagTag (https://github.com/malonge/RagTag) and and annotation transfered to them by using Liftoff (https://github.com/agshumate/Liftoff). In paralel large scale alterations affecting the regions of interest are evaluated by synteny block conservation through Syri (https://github.com/schneebergerlab/syri). The researcher then needs to define the definitive location of their region of interest: 1) The transfered on locations that are expected to be part of the region; and 2) Syri's synteny blocks containing the region. Then Prepare a different input table that contains the exact coordinates to include for each assembly.

The script uses the same reference for scaffolding as the genome assembly where the region is defined, but at the end of this stage this one can be switched for a different one for one or more assemblies if deemed beneficial. Since different references can be assembled in different orientations, this can be set up individually for each assembly. Furthermore, a region can be set up to include various fragments but caution is advised. If a region was recovered fragmented in an assembly it is so for a reasos, that be biological or as a technical artifacto.

Next (or in paralel) the reference proteins need to be re-annotated using Interproscan (https://www.ebi.ac.uk/interpro/about/interproscan/), Signalp (https://services.healthtech.dtu.dk/services/SignalP-6.0/), tmhmm (https://services.healthtech.dtu.dk/services/TMHMM-2.0) and TargetP (https://services.healthtech.dtu.dk/services/TargetP-2.0/). The latter could not be installed in our systems for reasons beyond our control, so instead the script utilizes Seqkit (https://bioinf.shenwei.me/seqkit/) to prepare intput files for their web service. Then a separate code block sorts the output to be used for the rest of the pipeline by merging all input tables and adding Gene IDs to each result.

Lastly, the script retrieves and sorts out all available information for the genes located in the region, according to each annotation.

### Gene_Model_Re-annotation.sh
This script takes the results of the main script, sorts genes into gene models or groups following a table provided by the user, and then runs interproscan to produce json files for visual inspection. 

### Dotplot_Analysis_1.sh and Dotplot_Analysis_2.sh
The first script takes the table with the refined region coordinates and orientation and extracts the sequences for latter use by other tools. In our analysis we used D-genies for an initial exploration of the region and quickly evaluate what reference was better suited for the final reconstruction of the region. The second script takes all sequnces in a folder and runs Blastn2dotplots (https://github.com/mokuno3430/blastn2dotplots) through all combinations. Highlights can be added to a single tab files that is parsed bu the script.

### Structural_Variant_Analysis.sh
This script various genomes and defines variants with Minimap2 and paftools.js call, then it preditcs their physiological impact according to each of the annotations transfered to the reference assembly using SnpEFF (https://pcingola.github.io/SnpEff/). Then the script sorts them out and count them using 2 tables provided by the user: 1) Gene Models IDs (same as Gene_Model_Re-annotation.sh), and 2) a table specifying the gene groups that should be analyzed together on what analysis and what is the treshold for counting variants asigned to each group.

Because of the high complexity of potential overlapping gene models caused by artifacts from the annotation transfers plus real overlapping genes. The researcher is required to review the Raw counts and set on these cases if the variant is going to be counted by one of the gene models, both or neither. While daunting on the surface given the files size, in our expirience working on these regions this is only necesary on very concrete places.  

### Other code quirks
1) Liftoff proved hard to paralelize, and with references to use for both scaffolding and annotation transfers, this proved prohibitively time consuming (1 query assembly translates into 16 annotation transfers). Our incomplete work arround was to have 5 copies of each genome sequence file and annotation and stagger their use so not two files are used at the same time. On its current implementation the problem still persist in a low frequency, but the solution is to just try again.

2) For the purposes of variants counts belonging to one group or another (script Structural_Variant_Analysis.sh), a lot of internal back and forth took place if it was better to use a flat number or a percentage of the group size. To accomodate the script takes on both. If the value provided is bellow 1 it will be treated as a percentage, above 1 it will be a set threshold. Please do not provide a float number bigger than 1.



