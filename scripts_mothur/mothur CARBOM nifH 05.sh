#==========================================================================================================================
#                ================   nifH ================ 
#==========================================================================================================================

#==========================================================================================================================
#      ================ BUILDING THE REFERENCE DATABASE AND ALIGNEMENTS ================ 
#==========================================================================================================================

# Start ARB
/usr/local/genome2/arb-6.0.2/arb/bin/arb

#===================================================
#       # nifH -  Use alignement from 2014 Gaby/Buckley et al. paper
#===================================================
cd /projet/umr7144/dipo/vaulot/mothur/database/nifH

# NOT DONE Extract region corresponding to V4 Euks with 0 differences for primers

# Intitial 32954 Final 9654
pcr.seqs(fasta=nifH_gaby_align.fasta, oligos=oligosnifH.oligos, pdiffs=2, keepdots=F, processors=8)

# Rename output files to nifH_gaby_align.pcr.nifh.pd2.align.fasta

# The next line is needed because all sequences do not have same length final : 8042
screen.seqs(fasta=nifH_gaby_align.pcr.nifh.pd2.align.fasta, start=1, end=679, maxlength=500, maxambig=0)

# Remove column having only gaps
filter.seqs(fasta=nifH_gaby_align.pcr.nifh.pd2.align.good.fasta)

# Reduce size of data by having only unique seqs : 4740
unique.seqs(fasta=nifH_gaby_align.pcr.nifh.pd2.align.good.filter.fasta)


#===================================================
#       # nifH -  Use the Zehr database - NOT USED
#===================================================

gzip -d nifH_zehr.gz

# Remove column having only gaps
mothur "#filter.seqs(fasta=nifH_zehr.fasta)"
# [ERROR]: Sequences are not all the same length, please correct.

mothur "#summary.seqs(fasta=nifH_zehr.fasta)"
                # Start   End     NBases  Ambigs  Polymer NumSeqs
# Minimum:        237     1998    156     0       3       1
# 2.5%-tile:      958     2331    276     0       3       275
# 25%-tile:       1105    2415    324     0       4       2748
# Median:         1132    2418    327     0       4       5495
# 75%-tile:       1132    2433    360     0       4       8242
# 97.5%-tile:     1165    2886    594     1       6       10714
# Maximum:        1216    52308   2241    66      66      10988
# Mean:   1110.56 2722.47 357.798 0.277212        4.35821
# of Seqs:      10988



#==========================================================================================================================
#      ================ PROCESS nifH ================ !! USED with sequences obtained with trim.seqs - Redo with pcr.seqs
#==========================================================================================================================

# rename the files
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/nifH

mv carbom_file_names_clean.trim.contigs.good.pcr.nifH.fasta carbom.nifH.fasta 
mv carbom_file_names_clean.contigs.good.pcr.nifH.groups carbom.nifH.groups

################################################################################
## 3. Remove sequences shorter than 300 bp
################################################################################

mothur "#screen.seqs(fasta=carbom.nifH.fasta, group=carbom.nifH.groups, minlength=300, maxambig=0, processors=32)"
# Output File Names:
# carbom.nifH.good.fasta
# carbom.nifH.bad.accnos
# carbom.nifH.good.groups
# It took 531 secs to screen 1490889 sequences.

mv carbom.nifH.good.fasta carbom.nifH.fasta 
mv carbom.nifH.good.groups carbom.nifH.groups

################################################################################
##4a. unique seqs
################################################################################
mothur "#unique.seqs(fasta=carbom.nifH.fasta)"

# Output File Names:
# carbom.nifH.names
# carbom.nifH.unique.fasta

mothur "#summary.seqs(fasta=carbom.nifH.unique.fasta, processors=32)"
################################################################################
##4b. count sequences
################################################################################

# counting sequences
mothur "#count.seqs(name=carbom.nifH.names, group=carbom.nifH.groups, processors=32)"
# Output File Names:
# carbom.nifH.count_table

# counting sequences per sample
mothur "#count.groups(count=carbom.nifH.count_table)"
# Output File Names:
# carbom.nifH.count.summary

################################################################################
##4c. remove sequences singletons 
################################################################################

cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/nifH

# remove singletons
mothur "#split.abund(count=carbom.nifH.count_table, fasta=carbom.nifH.unique.fasta, cutoff=1, accnos=true)"
# Output File Names:
# carbom.nifH.rare.count_table
# carbom.nifH.abund.count_table
# rare.accnos
# abund.accnos
# carbom.nifH.unique.rare.fasta
# carbom.nifH.unique.abund.fasta

# compute the abundance of the non-singletons
mothur "#count.groups(count=carbom.nifH.abund.count_table)"
wc -l carbom.nifH.unique.abund.fasta

################################################################################
## 5. Align to nifH  and filter
################################################################################
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/nifH

mothur "#align.seqs(fasta=carbom.nifH.unique.abund.fasta, reference=/projet/umr7144/dipo/vaulot/mothur/database/nifH/nifH_buckley.fasta, flip=T, processors=32)"
# [WARNING]: Some of your sequences generated alignments that eliminated too many bases, a list is provided in carbom.nifH.unique.abund.flip.accnos. If the reverse compliment proved to be better it was reported.
# It took 11 secs to align 34046 sequences.
# Output File Names:
# carbom.nifH.unique.abund.align
# carbom.nifH.unique.abund.align.report
# carbom.nifH.unique.abund.flip.accnos


# display the first 250 sequences to check the aliognment
head -500 carbom.nifH.unique.abund.align > carbom.nifH.unique.abund.align.extract

# Filter bad characters
mothur "#filter.seqs(fasta=carbom.nifH.unique.abund.align, processors=32)"
# Length of filtered alignment: 453
# Number of columns removed: 1319
# Length of the original alignment: 1772
# Number of sequences used to construct filter: 34046
# Output File Names:
# carbom.filter
# carbom.nifH.unique.abund.filter.fasta

head -500 carbom.nifH.unique.abund.filter.fasta > carbom.nifH.unique.abund.filter.fasta.extract
# Note still has the . character at the end of the file

################################################################################
## 6. Precluster with 2 differences the whole sample - Do by group much faster
################################################################################

mothur "#pre.cluster(fasta=carbom.nifH.unique.abund.filter.fasta, count=carbom.nifH.abund.count_table, diffs=2, processors=32)"

# Output File Names:
# carbom.nifH.unique.abund.filter.precluster.fasta
# carbom.nifH.unique.abund.filter.precluster.count_table
# carbom.nifH.unique.abund.filter.precluster.10n.map
# carbom.nifH.unique.abund.filter.precluster.10p.map...

wc -l carbom.nifH.unique.abund.filter.precluster.fasta

################################################################################
## 7. Remove chimeras - Two methods uchime and vsearch (new)
#    vsearch did not work
################################################################################

# uchime
mothur "#chimera.uchime(fasta=carbom.nifH.unique.abund.filter.precluster.fasta, count=carbom.nifH.unique.abund.filter.precluster.count_table, processors=32)"

# Output File Names:
# carbom.nifH.unique.abund.filter.precluster.denovo.uchime.chimeras
# carbom.nifH.unique.abund.filter.precluster.denovo.uchime.accnos

# vsearch - NOT USED
mothur "#chimera.vsearch(fasta=carbom.nifH.unique.abund.filter.precluster.fasta, count=carbom.nifH.unique.abund.filter.precluster.count_table, processors=32)"

# remove chimeras
mothur "#remove.seqs(fasta=carbom.nifH.unique.abund.filter.precluster.fasta, accnos=carbom.nifH.unique.abund.filter.precluster.denovo.uchime.accnos, count=carbom.nifH.unique.abund.filter.precluster.count_table)"

# Removed 155 sequences from your fasta file.
# Removed 1693 sequences from your count file.
# Output File Names:
# carbom.nifH.unique.abund.filter.precluster.pick.fasta
# carbom.nifH.unique.abund.filter.precluster.pick.count_table


mothur "#summary.seqs(fasta=carbom.nifH.unique.abund.filter.precluster.pick.fasta)"
# Output File Names:
# carbom.nifH.unique.abund.filter.precluster.pick.summary

# recount sequence in each sample
mothur "#count.groups(count=carbom.nifH.unique.abund.filter.precluster.pick.count_table)"
# Output File Names:
# carbom.nifH.unique.abund.filter.precluster.pick.count.summary

################################################################################
## 8 Split sequences between controls and samples
################################################################################
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/nifH

# samples
mothur "#get.groups(fasta=carbom.nifH.unique.abund.filter.precluster.pick.fasta,count=carbom.nifH.unique.abund.filter.precluster.pick.count_table, accnos=carbom_groups_samples.txt)"
# Selected 2490 sequences from your fasta file.
# Selected 1198040 sequences from your count file.
# Output File names:
# carbom.nifH.unique.abund.filter.precluster.pick.pick.fasta
# carbom.nifH.unique.abund.filter.precluster.pick.pick.count_table

mv carbom.nifH.unique.abund.filter.precluster.pick.pick.fasta carbom.nifH.unique.abund.filter.precluster.pick.samples.fasta
mv carbom.nifH.unique.abund.filter.precluster.pick.pick.count_table carbom.nifH.unique.abund.filter.precluster.pick.samples.count_table
mothur "#count.groups(count=carbom.nifH.unique.abund.filter.precluster.pick.samples.count_table)"

# controls
mothur "#get.groups(fasta=carbom.nifH.unique.abund.filter.precluster.pick.fasta,count=carbom.nifH.unique.abund.filter.precluster.pick.count_table, accnos=carbom_groups_controls.txt)"
# Output File names:
# carbom.nifH.unique.abund.filter.precluster.pick.pick.fasta
# carbom.nifH.unique.abund.filter.precluster.pick.pick.count_table

mv carbom.nifH.unique.abund.filter.precluster.pick.pick.fasta carbom.nifH.unique.abund.filter.precluster.pick.controls.fasta
mv carbom.nifH.unique.abund.filter.precluster.pick.pick.count_table carbom.nifH.unique.abund.filter.precluster.pick.controls.count_table

################################################################################
##                             SAMPLES  
################################################################################


################################################################################
## 9  Remove low abundances sequences from samples - in our case we remove sequences with less than < 10  
################################################################################
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/nifH

mothur "#split.abund(count=carbom.nifH.unique.abund.filter.precluster.pick.samples.count_table, fasta=carbom.nifH.unique.abund.filter.precluster.pick.samples.fasta, cutoff=10, accnos=true)"
# Output File Names:
# carbom.nifH.unique.abund.filter.precluster.pick.samples.rare.count_table
# carbom.nifH.unique.abund.filter.precluster.pick.samples.abund.count_table
# rare.accnos
# abund.accnos
# carbom.nifH.unique.abund.filter.precluster.pick.samples.rare.fasta
# carbom.nifH.unique.abund.filter.precluster.pick.samples.abund.fasta

cp carbom.nifH.unique.abund.filter.precluster.pick.samples.abund.fasta        carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.fasta
cp carbom.nifH.unique.abund.filter.precluster.pick.samples.abund.count_table carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.count_table

mothur "#summary.seqs(fasta=carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.fasta)"

# Output File Names:
# carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.summary
# It took 0 secs to summarize 296 sequences.

# recount number of sequences in each sample after removing low abundance sequences
mothur "#count.groups(count=carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.count_table)"
# Output File Names:
# carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.count.summary


################################################################################
##11 Cluster and do OTUs - Only on abundant sequences (nseq > 10)
################################################################################

mothur "#dist.seqs(fasta=carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.fasta, processors=32)"
# Output File Names:
# carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.dist
# It took 1 seconds to calculate the distances for 296 sequences.


mothur "#cluster(column=carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.dist, count=carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.count_table)"
# It took 20 seconds to cluster
# Output File Names:
# carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.list


mothur "#get.oturep(fasta=carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.fasta, column=carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.dist, count=carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.count_table, list=carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.list, cutoff=0.02)"
# Output File Names:
# carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.0.01.rep.count_table
# carbom.nifH.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.0.01.rep.fasta




################################################################################
##                             CONTROLS  
################################################################################


################################################################################
## 9  Remove low abundances sequences from controls - in our case we remove sequences with less than < 10  
################################################################################
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/nifH

mothur "#split.abund(count=carbom.nifH.unique.abund.filter.precluster.pick.controls.count_table, fasta=carbom.nifH.unique.abund.filter.precluster.pick.controls.fasta, cutoff=10, accnos=true)"
# Output File Names:
# carbom.nifH.unique.abund.filter.precluster.pick.controls.rare.count_table
# carbom.nifH.unique.abund.filter.precluster.pick.controls.abund.count_table
# rare.accnos
# abund.accnos
# carbom.nifH.unique.abund.filter.precluster.pick.controls.rare.fasta
# carbom.nifH.unique.abund.filter.precluster.pick.controls.abund.fasta


cp carbom.nifH.unique.abund.filter.precluster.pick.controls.abund.fasta        carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.fasta
cp carbom.nifH.unique.abund.filter.precluster.pick.controls.abund.count_table  carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.count_table

mothur "#summary.seqs(fasta=carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.fasta)"

# Output File Names:
# carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.summary
# It took 0 secs to summarize 45 sequences.

# recount number of sequences in each sample after removing low abundance sequences
mothur "#count.groups(count=carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.count_table)"
# Total seqs: 2802466.
# Output File Names:
# carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.count.summary

################################################################################
##11 Cluster and do OTUs - only on abundant sequences (nseq > 10)
################################################################################

mothur "#dist.seqs(fasta=carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.fasta, processors=32)"
# Output File Names:
# carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.dist
# It took 0 seconds to calculate the distances for 45 sequences.

mothur "#cluster(column=carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.dist, count=carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.count_table)"
# It took 38 seconds to cluster
# Output File Names:
# carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.an.unique_list.list

mothur "#get.oturep(fasta=carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.fasta, column=carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.dist, count=carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.count_table, list=carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.an.unique_list.list, cutoff=0.01, sorted=bin)"
# Output File Names:
# carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.an.unique_list.0.01.rep.count_table
# carbom.nifH.uniq.preclust.no_chim.more_than_10.controls.an.unique_list.0.01.rep.fasta































################################################################################
##4. Etapes similaires au 18S
################################################################################
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/nifH

unique.seqs(fasta=carbom_file_names_clean.trim.contigs.good.trim.nifH.fasta)

align.seqs(fasta=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.fasta, reference=/projet/umr7144/dipo/cgerikasribeiro/mothur/database/nifH/nifH_gaby_align.pcr.nifh.pd2.final.align.fasta, flip=T, processors=8)

filter.seqs(fasta=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.align)

 pre.cluster(fasta=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.fasta, name=carbom_file_names_clean.trim.contigs.good.trim.nifH.names, group=carbom_file_names_clean.trim.contigs.good.trim.nifH.groups, diffs=2)
 
 chimera.uchime(fasta=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.fasta, name=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.names, group=carbom_file_names_clean.trim.contigs.good.trim.nifH.groups)
 
 remove.seqs(fasta=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.fasta, accnos=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.uchime.accnos, name=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.names, group=carbom_file_names_clean.trim.contigs.good.trim.nifH.groups)

################################################################################
##9a Count sequences and use count file
################################################################################
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/nifH

/usr/local/genome2/mothur-1.35.1/mothur "#count.seqs(name=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.names, group=carbom_file_names_clean.trim.contigs.good.trim.nifH.pick.groups, processors=8)"


################################################################################
##8.1. Split sequences between controls and samples
################################################################################
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/nifH

/usr/local/genome2/mothur-1.35.1/mothur "#get.groups(fasta=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.fasta, count=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.count_table, accnos=carbom_groups_samples.txt)"

# Must rename file to .sample.

/usr/local/genome2/mothur-1.35.1/mothur "#get.groups(fasta=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.fasta, count=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.count_table, accnos=carbom_groups_controls.txt)"

# Must rename file to .controls.

################################################################################
##9a Remove low abundances sequences - in our case we remove sequences with less than < 10  
################################################################################

/usr/local/genome2/mothur-1.35.1/mothur "#split.abund(count=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.samples.count_table, fasta=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.samples.fasta, cutoff=10, accnos=true)"

################################################################################
##10 Cluster and do OTUs - ONLY on SAMPLES and only on abundant sequences (nseq > 10)
################################################################################

/usr/local/genome2/mothur-1.35.1/mothur "#dist.seqs(fasta=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.samples.abund.fasta, processors=8)"

/usr/local/genome2/mothur-1.35.1/mothur "#cluster(column=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.samples.abund.dist, count=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.samples.abund.count_table)"

/usr/local/genome2/mothur-1.35.1/mothur "#get.oturep(fasta=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.samples.abund.fasta, column=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.samples.abund.dist, count=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.samples.abund.count_table, list=carbom_file_names_clean.trim.contigs.good.trim.nifH.unique.filter.precluster.pick.samples.abund.an.unique_list.list, cutoff=0.01, sorted=size)"

