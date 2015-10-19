import sys
import os
import subprocess as sp
import shutil
import numpy as np

'''
python 2_sensitivity_calculator_v1.py [blastouts directory] [path to '7-ExonCaptureEvaluationPhylo' script] [Path to probes file {ex. H_riggenbachi_probes.FASTA}]

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
'''

blast_directory = sys.argv[1]
os.chdir(blast_directory)

perl_script = sys.argv[2]

fasta = sys.argv[3]
fh_fasta = open(fasta, 'r')

contig_dict = {}

for line in fh_fasta:
    line=line.strip()
    if line.startswith('>'):
        gene_name=line[1:]
        #*********this script may fail for you here, because of the way your fasta lines are written
        #Notice the underscores in the above examples, this is key for continuing
        gene_names = gene_name.split('_')
        name = gene_names[3]
        contig_dict[name]=''
    else:
        bp_count = int(0)
        for base in line:
            bp_count += 1
        contig_dict[name] += str(bp_count)

print contig_dict     

def perc_calc(x,y):
    percentage = float( ( float(x) / float(y)) * float(100))
    percentage = str(np.around(percentage, decimals = 2))
    return percentage

def quickstats(x):
    x_array = np.asarray(x,dtype=np.float64)
    x_avg = np.average(x_array)
    x_avg = str(np.around(x_avg, decimals = 1))
    return x_avg

out_directory = 'Trimmed_blastouts'

if not os.path.exists(out_directory):
    os.mkdir(out_directory)

for fl in os.listdir('.'):
    if fl.startswith('blastout'):
        print "Trimming down {} to correct merging format...".format(fl)

        out_name = "trimmed_"+fl
        fh_out = open(out_name, 'a')

        fh_temp = open(fl, 'r')
        lines = fh_temp.readlines()

        for line in lines:
            line = line.strip()
            line = line.split('\t')
            qname = line[0]
            start = line[6]
            stop = line[7]
            fh_out.write(qname+'\t'+start+'\t'+stop+'\n')

        fh_out.close()
        fh_temp.close()

        for fl in os.listdir('.'):
            if fl.startswith('trimmed_'):
                shutil.move(fl, out_directory)
            
print '\n', "Beginning merge portion of script", '\n'

os.chdir(out_directory)

for fl in os.listdir('.'):
    if fl.startswith('trimmed_'):
        print "starting merge script for {}.".format(fl)
        perl_string = "perl {0} MakeBed -b {1} -r 0".format(perl_script, fl)
        proc_perl = sp.call(perl_string, shell=True)

contig_set = set()

for fl in os.listdir('.'):
    if fl.endswith('.txt_final'):
        temp_fh = open(fl, 'r')
        for line in temp_fh:
            line = line.strip()
            line = line.split('\t')
            names = line[0].split('_')
            ens = names[3]
            contig_set.add(ens)
        temp_fh.close()

contig_list = list(contig_set)

fh_contigout = open('out_contig_recovery.txt','a')
fh_contigout.write("Contig"+'\t'+"Avg_recovery"+'\n')

counter = int(0)
for contig in contig_list:
    counter += 1
    print "{0}. Starting analysis of {1}.".format(counter, contig)
    contig_avgs = []
    for fl in os.listdir('.'):
        if fl.endswith('.txt_final'):
            prefix = fl.split('.')
            prefixes = prefix[0].split('_')
            out_name = prefixes[2]+'_'+prefixes[3]+'_lengths.out'

            out_fh = open(out_name, 'a')
            temp_fh = open(fl, 'r')

            add_list = []
            for line in temp_fh:
                line = line.strip()
                line = line.split('\t')
                start = line[1]
                end = line[2]
                names = line[0].split('_')
                ens = names[3]
                if contig == ens:
                    bp = int(end) - int(start)
                    add_list.append(bp)
                    
            fragslen = int(0)
            for l in add_list:
                fragslen += int(l)
            frag_perc = perc_calc(fragslen, contig_dict[contig])
            contig_avgs.append(frag_perc)
            out_fh.write(contig+'\t'+"{}".format(fragslen)+'\t'+"{}".format(contig_dict[contig])+'\t'+frag_perc+'\n')
            out_fh.close()
            
    contig_avg = quickstats(contig_avgs)
    fh_contigout.write(contig+'\t'+contig_avg+'\n')

fh_contigout.close()

fh_master = open('Out_Master_Recovery.txt', 'a')
fh_master.write("Sample"+'\t'+"total_bases"+'\t'+"recovered_bases"+'\t'+"perc_recovery"+'\n')

count = int(0)
for fl in os.listdir('.'):
    if fl.endswith('_lengths.out'):
        count += 1
        print "{0}. Final calculation of recovered fragments for {1}.".format(count,fl)
        names = fl.split('_lengths.')
        sample = names[0]
        t_bases = int(0)
        r_bases = int(0)
        fh_temp = open(fl, 'r')
        for line in fh_temp:
            line = line.strip()
            line = line.split('\t')
            t_bases += int(line[2])
            r_bases += int(line[1])
        rec_perc = perc_calc(r_bases,t_bases)
        fh_master.write(sample+'\t'+"{}".format(t_bases)+'\t'+"{}".format(r_bases)+'\t'+"{}".format(rec_perc)+'\n')

fh_master.close()          

