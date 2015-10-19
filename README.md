# Sensitivity_Calculator

Overall workflow is to calculate sensitivity using the recovered contigs and sequences used for probe design.  
First step is to blast the probe design sequences against the contigs of each sample. The coordinates
recovered are merged if they overlap for particular contigs, then the total bases per contig receiving a blast hit
are calculated.  The corresponding contig bases are counted in the original sequence design file, and the
fraction of bases receiving hits are calculated (=sensitivity).  This is dependent on the outputs of 
the https://github.com/MVZSEQ/denovoTargetCapturePhylogenomics workflow.


STEP 1 - 1_contig_blaster_v1.py
Usage: python 1_contig_blaster_v1.py [directory with JMSR008_index49_targetedRegionAndFlanking.fasta file(s)] 
[Path to probes file {ex. 'H_riggenbachi_probes.FASTA'}] [output directory name]

Create an output directory to store blast outputs before starting (will be used again in second script).

The 'Sample_indexXX_targetedRegionAndFlanking.fasta' files are located in the 
'intarget_assemblies' folder resulting from script 6 (6-TransExonCapPhylo 'contig' option) 
of the https://github.com/MVZSEQ/denovoTargetCapturePhylogenomics workflow. 

The file used for probe design needs to be in fasta format, and will be used over again
for every sample fasta file included.  This file should contain the SEQUENCES used for the probe 
design, NOT the probe sequences themselves.  If multiple species were used for probe design 
(ie more than one sequence per 'gene'), reduce the fasta file to only include one species.

Script will make a blastn database from each Sample_indexXX_targetedRegionAndFlanking.fasta file,
and query the probes fasta file to the sample database.  Cleans up files each step and moves
the output to the directory that you make and indicate.

'Sample_indexXX_targetedRegionAndFlanking.fasta' format used:
>JMPD002_index1_Contig1
GATGACTCTGTTATTATGAAAGGCAATGCGGGAGTACCATGTCAGCTTTTTCTACAAACTTTGTTGTT
>JMPD002_index1_Contig10
AGGGCATCCCGGCCCCTGATTGGTCCTCCAGCAGAGGGAATCTCTGTAGCGTGCTGTATGGGTGTCAC

Probes fasta file format:
>Hyp_riggenbachi_Contig2871_ENSXETP00000004106_F7
GGACTTTGGAGTGGCAGCATTCTCCTGAGAGCAGAGCTGAAATAGGAAAACCTCCCCATCATGGGGCT
>Hyp_riggenbachi_Contig2983_ENSXETP00000022229_CCT6A
CTAGATTGTTCTGGAAGAGCCGGACGCCGCGTGGAGGAAAGGCCGGGATCGTTGCTTTCTGTACACGT

'blastout_Sample_indexXX.txt' format (blast output style 6):
Hyp_riggenbachi_Contig2871_ENSXETP00000004106_F7	JMPD003_index9_Contig1	95.32	171	8	0	127	297	772	942	3e-73	  272
Hyp_riggenbachi_Contig2871_ENSXETP00000004106_F7	JMPD003_index9_Contig1	97.10	138	4	0	437	574	2235	2372	1e-61	  233


STEP 2 - 2_sensitivity_calculator_v1.py

Usage: python 2_sensitivity_calculator_v1.py [blastouts directory] [path to '7-ExonCaptureEvaluationPhylo' script] [Path to probes file {ex. H_riggenbachi_probes.FASTA}]

Use the output directory from 1_contig_blaster_v1.py to sift through the 'blastout_Sample_indexXX.txt' files.  
This script requires the 7-ExonCaptureEvaluationPhylo script from the https://github.com/MVZSEQ/denovoTargetCapturePhylogenomics workflow.
Use the same probes fasta file as in the last step.  Creates a 'Trimmed_blastouts' directory to
move updated blast files, merge overlapping blast coordinates, then calculate sensitivity for each
sample (producing 'Sample_indexXX_lengths.out' output file) and a master output file 'Out_Master_Recovery.txt'.

starting files from previous script 'blastout_Sample_indexXX.txt' format (blast output style 6):
Hyp_riggenbachi_Contig2871_ENSXETP00000004106_F7	JMPD003_index9_Contig1	95.32	171	8	0	127	297	772	942	3e-73	  272
Hyp_riggenbachi_Contig2871_ENSXETP00000004106_F7	JMPD003_index9_Contig1	97.10	138	4	0	437	574	2235	2372	1e-61	  233

intermediate file 'trimmed_blastout_Sample_indexXX.txt_final' (from 7-ExonCaptureEvaluationPhylo
processing) shows merged blast coordinates:
Hyp_riggenbachi_Contig1000_ENSXETP00000062181_SLC7A7	37	850
Hyp_riggenbachi_Contig1006_ENSXETP00000061034_NOP56	1	567
Hyp_riggenbachi_Contig1006_ENSXETP00000061034_NOP56	569	850

Final individual sensitivity calculation per sample ('Sample_indexXX_lengths.out'):
ENSXETP00000036704	833	850	98.0
ENSXETP00000038652	513	526	97.53
ENSXETP00000059911	847	850	99.65
ENSXETP00000019319	476	687	69.29
ENSXETP00000048091	587	850	69.06
ENSXETP00000014695	348	507	68.64

And the final output file you'll want to inspect ('Out_Master_Recovery.txt'):
Sample	total_bases	recovered_bases	perc_recovery
JMPD002_index10	995246	867840	87.2
JMPD002_index11	995246	886496	89.07
JMPD002_index12	995246	874511	87.87
JMPD002_index13	995246	913305	91.77
JMPD002_index14	995246	891390	89.56
JMPD002_index15	995246	890485	89.47
JMPD002_index16	995246	868389	87.25

The 'perc_recovery' column above indicates the sensitivity (%). 

##############
DEPENDENCIES:
numpy - Numerical Python
##############
------------------------
written for Python 2.7.3
Dan Portik
daniel.portik@berkeley.edu
October 2015
------------------------
