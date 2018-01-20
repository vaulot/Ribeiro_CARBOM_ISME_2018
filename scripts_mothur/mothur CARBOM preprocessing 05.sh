#==========================================================================================================================
#      CARBOM - On the server  
#==========================================================================================================================


#==========================================================================================================================
#      ================  PROCESSING THE RAW SEQUENCES ===================================================
#==========================================================================================================================

#===================================================
#      # 1a -  Save the data to the carbom folder
#===================================================

cp -r /projet/umr7144/dipo/cgerikasribeiro/data /projet/sbr/carbom/data

#===================================================
#      # 1b -  Unzip and clean the sequences
#===================================================
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/fastq   

#  Rename all the files

rename "18S-" ""  *.gz
rename "_L001" ""  *.gz
rename "_001" ""  *.gz

#  Unzip all the files

gzip -d *.gz

# Remove all sequences shorter than 200 bp
# Quality limit is 20 and at least 75% of bases must be above this quality

R1="R1"
R2="R2"
for NAME in $(cat carbom_file_names.txt); do
	cat $NAME$R1.fastq|fastq_quality_trimmer -Q 33 -t 25 -l 200 -v -o  $NAME$R1.clean1.fastq >$NAME$R1.clean1.rpt
	cat $NAME$R1.clean1.fastq|fastq_quality_filter -Q 33 -q 20 -p 75 -v -o $NAME$R1.clean2.fastq>$NAME$R1.clean2.rpt
	cat $NAME$R2.fastq|fastq_quality_trimmer -Q 33 -t 25 -l 200 -v -o  $NAME$R2.clean1.fastq >$NAME$R2.clean1.rpt
	cat $NAME$R2.clean1.fastq|fastq_quality_filter -Q 33 -q 20 -p 75 -v -o $NAME$R2.clean2.fastq>$NAME$R2.clean2.rpt
	/usr/local/genome2/scripts/get_pairs.py $NAME$R1.clean2.fastq $NAME$R2.clean2.fastq>$NAME.get_pairs.rpt
done;


	
#===================================================
# cat 10n_S47_L001_R1_001.fastq|fastq_quality_trimmer -Q 33 -t 25 -l 200 -v|fastq_quality_filter -Q 33 -q 20 -p 75 -v -o 10n_S47_L001_R1_001.clean.fastq
#[-t N] = Quality threshold - nucleotides with lower #quality will be trimmed (from the end of the sequence).
#[-l N] = Minimum length - sequences shorter than this (after trimming) will be discarded. Default = 0 = no minimum length.
#===================================================
FASTA/Q Trimmer

	$ fastx_trimmer -h
	usage: fastx_trimmer [-h] [-f N] [-l N] [-z] [-v] [-i INFILE] [-o OUTFILE]

	version 0.0.6
	   [-h]         = This helpful help screen.
	   [-f N]       = First base to keep. Default is 1 (=first base).
	   [-l N]       = Last base to keep. Default is entire read.
	   [-z]         = Compress output with GZIP.
	   [-i INFILE]  = FASTA/Q input file. default is STDIN.
	   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.

FASTQ Quality Filter

$ fastq_quality_filter -h
	usage: fastq_quality_filter [-h] [-v] [-q N] [-p N] [-z] [-i INFILE] [-o OUTFILE]

	version 0.0.6
	   [-h]         = This helpful help screen.
	   [-q N]       = Minimum quality score to keep.
	   [-p N]       = Minimum percent of bases that must have [-q] quality.
	   [-z]         = Compress output with GZIP.
	   [-i INFILE]  = FASTA/Q input file. default is STDIN.
	   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.
	   [-v]         = Verbose - report number of sequences.
			  If [-o] is specified,  report will be printed to STDOUT.
			  If [-o] is not specified (and output goes to STDOUT),
			  report will be printed to STDERR.

# to remove sequences if no pair in the corresponding sample, for each sample (R1 = forward, R2 = reverse):
/usr/local/genome2/scripts/get_pairs.py 10n_S47_L001_R1_001.clean.fastq 10n_S47_L001_R2_001.clean.fastq


################################################################################
##3. Make contigs and remove sequences with any ambiguity - to redo from here with all the sequences
################################################################################

# The oligo file can only be used with barcodes

cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/fastq_clean   

/usr/local/genome2/mothur-1.35.1/mothur "#make.contigs(file=carbom_file_names_clean.txt, processors=32)"

/usr/local/genome2/mothur-1.35.1/mothur "#summary.seqs(fasta=carbom_file_names_clean.trim.contigs.fasta)"

/usr/local/genome2/mothur-1.35.1/mothur "#screen.seqs(fasta=carbom_file_names_clean.trim.contigs.fasta,group=carbom_file_names_clean.contigs.groups,maxambig=0,processors=32)"

################################################################################
##3a. Find primers with PCR seqs
################################################################################
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/contigs  

# rename and move file each time !!!!! to the direcotries 18S / 16S /nifH

/usr/local/genome2/mothur-1.35.1/mothur "#pcr.seqs(fasta=carbom_file_names_clean.trim.contigs.good.fasta, group=carbom_file_names_clean.contigs.good.groups, oligos=oligos18szig.oligos, pdiffs=2, processors=8)"
/usr/local/genome2/mothur-1.35.1/mothur "#pcr.seqs(fasta=carbom_file_names_clean.trim.contigs.good.fasta, group=carbom_file_names_clean.contigs.good.groups, oligos=oligos16s.oligos, pdiffs=2, processors=8)"
mothur "#pcr.seqs(fasta=carbom_file_names_clean.trim.contigs.good.fasta, group=carbom_file_names_clean.contigs.good.groups, oligos=oligosnifH.oligos, pdiffs=0, processors=64)"

# Removed 4708978 sequences from your group file.
# Output File Names:
# carbom_file_names_clean.trim.contigs.good.pcr.fasta
# carbom_file_names_clean.trim.contigs.good.bad.accnos
# carbom_file_names_clean.trim.contigs.good.scrap.pcr.fasta
# carbom_file_names_clean.contigs.good.pcr.groups
# It took 2588 secs to screen 6199867 sequences.

mv carbom_file_names_clean.trim.contigs.good.pcr.fasta carbom_file_names_clean.trim.contigs.good.pcr.pd0.nifH.fasta
mv carbom_file_names_clean.contigs.good.pcr.groups carbom_file_names_clean.contigs.good.pcr.pd0.nifH.groups

mothur "#count.groups(group=carbom_file_names_clean.contigs.good.pcr.pd0.nifH.groups)"

################################################################################
##3b. Find primers with trim.seqs
#    -- Much faster but need to remove the sequences that are not selected from the group file -- also gives a bit more sequences....
################################################################################

/usr/local/genome2/mothur-1.35.1/mothur "#trim.seqs(fasta=carbom_file_names_clean.trim.contigs.good.fasta, oligos=oligosnifH.oligos, pdiffs=2, processors=1)"
# rename output
mv carbom_file_names_clean.trim.contigs.good.trim.fasta carbom_file_names_clean.trim.contigs.good.trim.nifH.fasta
# remove unwanted sequences from the group file
# Step 1. produce list of sequences accession -> carbom_file_names_clean.trim.contigs.good.trim.nifH.accnos
# Step 1. use this list to create a new group file with only the nifH sequences
mothur "#list.seqs(fasta=carbom_file_names_clean.trim.contigs.good.trim.nifH.fasta)"
mothur "#get.seqs(accnos=carbom_file_names_clean.trim.contigs.good.trim.nifH.accnos, group=carbom_file_names_clean.contigs.good.groups)"
mv carbom_file_names_clean.contigs.good.pick.groups carbom_file_names_clean.trim.contigs.good.trim.nifH.groups

mothur "#count.groups(group=carbom_file_names_clean.trim.contigs.good.trim.nifH.groups)"











