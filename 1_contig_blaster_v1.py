import sys
import os
import subprocess as sp
import shutil
import numpy as np
'''
Usage: python 1_contig_blaster_v1.py [directory with JMSR008_index49_targetedRegionAndFlanking.fasta file(s)] 
[Path to probes file {ex. 'H_riggenbachi_probes.FASTA'}] [output directory name]

Create an output directory to store blast outputs before starting (will be used again in second script).

The 'Sample_indexXX_targetedRegionAndFlanking.fasta' files are located in the 
'intarget_assemblies' folder resulting from script 6 (6-TransExonCapPhylo 'contig' option) 
of the https://github.com/MVZSEQ/denovoTargetCapturePhylogenomics workflow. 

The file used for probe design needs to be in fasta format, and will be used over again
for every sample fasta file included.  This file should contain the SEQUENCES used for the probe 
design, NOT the probe sequences themselves. If multiple species were used for probe design 
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



##############
DEPENDENCIES:
numpy - Numerical Python
blastn - needs to be in path to call from command line
##############
------------------------
written for Python 2.7.3
Dan Portik
daniel.portik@berkeley.edu
October 2015
------------------------
'''

fasta_directory = sys.argv[1]
os.chdir(fasta_directory)

probes_file = sys.argv[2]
fh_probes = open(probes_file, 'r')

output_directory = sys.argv[3]

for fl in os.listdir('.'):
    if fl.endswith('targetedRegionAndFlanking.fasta'):
        print "Beginning blast database for {}...".format(fl)

        db_string = "makeblastdb -in {} -dbtype nucl".format(fl)
        proc_makedb = sp.call(db_string, shell=True)
        
        print "Database finished, beginning blastn search."
        
        names = fl.split('_')
        prefix = names[0]+'_'+names[1]
        out_string = 'blastout_'+prefix+'.txt'

        call_string = "blastn -db {0} -query {1} -evalue 1e-10 -outfmt 6 -out {2}".format(fl,probes_file,out_string)
        proc_blastn = sp.call(call_string, shell=True)

        print "Cleaning up files...", '\n'
        for newfile in os.listdir('.'):
            if newfile.endswith('.fasta.nhr'):
                os.remove(newfile)
            if newfile.endswith('.fasta.nin'):
                os.remove(newfile)
            if newfile.endswith('.fasta.nsq'):
                os.remove(newfile)
            if newfile.startswith('blastout_'):
                shutil.move(newfile, output_directory)
