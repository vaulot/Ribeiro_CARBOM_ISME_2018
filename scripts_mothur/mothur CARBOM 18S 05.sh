#==========================================================================================================================
#                ================   18S ================ 
#==========================================================================================================================

#==========================================================================================================================
#     *** BUILDING THE REFERENCE DATABASE AND ALIGNEMENTS ***
#==========================================================================================================================




#======================================================================================================
#       18S  - PR2 - for taxonomy assignement _ Use gb203 version 4.4 - 178 096 sequences
#             Final files are pr2_v4_zingone.taxo and pr2_v4_zingone.fasta
#======================================================================================================

cd /projet/umr7144/dipo/vaulot/mothur/database/pr2_gb203_v_4.4/


# Need to have at least pd = 2 in order to catch the Haptophytes... 
# pdiffs=2 / 79227 sequences  - miss the haptophyta completely
# pdiffs=2 (forward) + rdiffs=2 (reverse) /  129026
mothur "#pcr.seqs(fasta=pr2_gb203_version_4.4.fasta, taxonomy=pr2_gb203_version_4.4.taxo, oligos=oligos18szig.oligos, pdiffs=2,  rdiffs=2, keepdots=F, processors=32)"

rename "pcr" "pcr.pd2.rd2.V4.zingone" * 

# remove sequences shorter than 300 and longer than 600 and containing any ambig - remain 114006
mothur "#screen.seqs(fasta=pr2_gb203_version_4.4.pcr.pd2.rd2.V4.zingone.fasta, minlength=300, maxlength=600, maxambig=0)"
mothur "#summary.seqs(fasta=pr2_gb203_version_4.4.pcr.pd2.rd2.V4.zingone.good.fasta)"

# only keep unique seqs - remain 61200
mothur "#unique.seqs(fasta=pr2_gb203_version_4.4.pcr.pd2.rd2.V4.zingone.good.fasta)"
mothur "#summary.seqs(fasta=pr2_gb203_version_4.4.pcr.pd2.rd2.V4.zingone.good.unique.fasta)"

# get the accession number of the unique sequences
# Output : pr2_gb203_version_4.4.pcr.pd2.rd2.V4.zingone.good.unique.accnos
mothur "#list.seqs(fasta=pr2_gb203_version_4.4.pcr.pd2.rd2.V4.zingone.good.unique.fasta)"

# get the taxonomy of the unique sequences
# Output File Names:  pr2_gb203_version_4.4.pcr.pd2.rd2.V4.zingone.pick.taxo
mothur "#get.seqs(accnos=pr2_gb203_version_4.4.pcr.pd2.rd2.V4.zingone.good.unique.accnos, taxonomy=pr2_gb203_version_4.4.pcr.pd2.rd2.V4.zingone.taxo)"

# copy the files to simpler name
cp pr2_gb203_version_4.4.pcr.pd2.rd2.V4.zingone.pick.taxo pr2_v4_zingone.taxo
cp pr2_gb203_version_4.4.pcr.pd2.rd2.V4.zingone.good.unique.fasta pr2_v4_zingone.fasta


#================================================================================================
#       # 18S -  Silva seed - for alignement and compute OTUs
#================================================================================================

cd /projet/umr7144/dipo/vaulot/mothur/database/silva_123

# uncompress archive
# files produced are 
# - silva.seed_v123.tax
# - silva.seed_v123.align
gzip -d  silva.seed_v123.tgz
tar xvf silva.seed_v123.tar

# Extract only eukaryotes
mothur "#get.lineage(taxonomy=silva.seed_v123.tax, taxon=Eukaryota, fasta=silva.seed_v123.align)"

# Remove column having only gaps
# Length of filtered alignment: 6886
# Number of columns removed: 43114
# Length of the original alignment: 50000
# Number of sequences used to construct filter: 2537
mothur "#filter.seqs(fasta=silva.seed_v123.pick.align, processors=32)"
# Remove column having only dots
mothur "#filter.seqs(fasta=silva.seed_v123.pick.filter.fasta, trump=., processors=32)"

mv silva.seed_v123.pick.filter.filter.fasta silva.seed_v123.euk.align
mv silva.seed_v123.pick.tax silva.seed_v123.euk.tax

# Extract region corresponding to V4 Euks with differences pdiffs=1,  rdiffs=2
mothur "#pcr.seqs(fasta=silva.seed_v123.euk.align, taxonomy=silva.seed_v123.euk.tax, oligos=oligos18szig.oligos, pdiffs=1,  rdiffs=2, keepdots=F, processors=32)"

# Here it is necessary to go to Geneious to remove bad sequences as well as sequences that are offset (if first base of primer is not present)
# A couple of bad seqeunces are also removed
# Files are renamed as 
#    silva.seed_v123.euk.V4.zingone.align
#    silva.seed_v123.euk.V4.zingone.tax


# Reduce size of data by having only unique seqs
mothur "#unique.seqs(fasta=silva.seed_v123.euk.V4.zingone.align)"

mothur "#summary.seqs(fasta=silva.seed_v123.euk.V4.zingone.unique.align)"


#                Start   End     NBases  Ambigs  Polymer NumSeqs
# Minimum:        1       2051    262     0       5       1
# 2.5%-tile:      1       2052    359     0       5       50
# 25%-tile:       1       2052    376     0       5       494
# Median:         1       2052    377     0       6       988
# 75%-tile:       1       2052    378     0       6       1482
# 97.5%-tile:     1       2052    407     1       6       1926
# Maximum:        2       2052    591     5       9       1975
# # of Seqs:      1975


#==========================================================================================================================
#                ================   18S ================ 
#==========================================================================================================================

# rename the files
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/18S

mv carbom_file_names_clean.trim.contigs.good.pcr.18S.fasta carbom.18S.V4.fasta 
mv carbom_file_names_clean.contigs.good.pcr.18S.groups carbom.18S.V4.groups



################################################################################
##4a. unique seqs
################################################################################
mothur "#unique.seqs(fasta=carbom.18S.V4.fasta)"

# Output File Names:
# carbom.18S.V4.names
# carbom.18S.V4.unique.fasta

mothur "#summary.seqs(fasta=carbom.18S.V4.unique.fasta, processors=32)"

################################################################################
##4b. remove sequences singletons : goes from 3788366 reads to 2937528 reads with 150 252 unique
#  If one use a threshold at 10 (goes from 3788366 reads to 2501737 reads with 20300 unique sequences - too much reads removed)
#  IF one does not remove any sequence then 1 001 090 unique sequences....)
################################################################################

cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/18S

# counting unique sequences
mothur "#count.seqs(name=carbom.18S.V4.names, group=carbom.18S.V4.groups, processors=32)"
# Output File Names:
# carbom.18S.V4.count_table

# counting unique sequences per sample
mothur "#count.groups(count=carbom.18S.V4.count_table)"
# Output File Names:
# carbom.18S.V4.count.summary

# remove singletons
mothur "#split.abund(count=carbom.18S.V4.count_table, fasta=carbom.18S.V4.unique.fasta, cutoff=1, accnos=true)"
# Output File Names:
# carbom.18S.V4.rare.count_table
# carbom.18S.V4.abund.count_table
# rare.accnos
# abund.accnos
# carbom.18S.V4.unique.rare.fasta
# carbom.18S.V4.unique.abund.fasta

# compute the abundance of the non-singletons
mothur "#count.groups(count=carbom.18S.V4.abund.count_table)"
wc -l carbom.18S.V4.unique.abund.fasta

################################################################################
## 5. Align to Silva and filter
################################################################################
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/18S

mothur "#align.seqs(fasta=carbom.18S.V4.unique.abund.fasta, reference=/projet/umr7144/dipo/vaulot/mothur/database/silva_123/silva.seed_v123.euk.V4.zingone.unique.align, flip=T, processors=32)"
# Output File Names:
# carbom.18S.V4.unique.abund.align
# carbom.18S.V4.unique.abund.align.report
# carbom.18S.V4.unique.abund.flip.accnos

# display the first 250 sequences to check the aliognment
head -500 carbom.18S.V4.unique.abund.align > carbom.18S.V4.unique.abund.align.extract

# Filter bad characters
mothur "#filter.seqs(fasta=carbom.18S.V4.unique.abund.align, processors=32)"
# Length of filtered alignment: 635
# Number of columns removed: 1418
# Length of the original alignment: 2053
# Number of sequences used to construct filter: 104254

# Output File Names:
# carbom.filter
# carbom.18S.V4.unique.abund.filter.fasta


head -500 carbom.18S.V4.unique.abund.filter.fasta > carbom.18S.V4.unique.abund.filter.fasta.extract
# Note still has the . character at the end of the file

################################################################################
## 6. Precluster with 2 differences the whole sample - Do by group much faster
################################################################################

mothur "#pre.cluster(fasta=carbom.18S.V4.unique.abund.filter.fasta, count=carbom.18S.V4.abund.count_table, diffs=2, processors=32)"

# Output File Names:
# carbom.18S.V4.unique.abund.filter.precluster.fasta
# carbom.18S.V4.unique.abund.filter.precluster.count_table
# carbom.18S.V4.unique.abund.filter.precluster.10n.map
# carbom.18S.V4.unique.abund.filter.precluster.10p.map...

wc -l carbom.18S.V4.unique.abund.filter.precluster.fasta

################################################################################
## 7. Remove chimeras - Two methods uchime and vsearch (new)
#    vsearch did not work
################################################################################

# uchime
mothur "#chimera.uchime(fasta=carbom.18S.V4.unique.abund.filter.precluster.fasta, count=carbom.18S.V4.unique.abund.filter.precluster.count_table, processors=32)"

# Output File Names:
# carbom.18S.V4.unique.abund.filter.precluster.denovo.uchime.chimeras
# carbom.18S.V4.unique.abund.filter.precluster.denovo.uchime.accnos

# vsearch - NOT USED
mothur "#chimera.vsearch(fasta=carbom.18S.V4.unique.abund.filter.precluster.fasta, count=carbom.18S.V4.unique.abund.filter.precluster.count_table, processors=32)"

# remove chimeras
mothur "#remove.seqs(fasta=carbom.18S.V4.unique.abund.filter.precluster.fasta, accnos=carbom.18S.V4.unique.abund.filter.precluster.denovo.uchime.accnos, count=carbom.18S.V4.unique.abund.filter.precluster.count_table)"
# Removed 3771 sequences from your fasta file.
# Removed 95796 sequences from your count file.
Output File Names:
# carbom.18S.V4.unique.abund.filter.precluster.pick.fasta
# carbom.18S.V4.unique.abund.filter.precluster.pick.count_table

mothur "#summary.seqs(fasta=carbom.18S.V4.unique.abund.filter.precluster.pick.fasta)"
# Output File Names:
# carbom.18S.V4.unique.abund.filter.precluster.pick.summary
# It took 2 secs to summarize 30709 sequences.

# recount sequence in each sample
mothur "#count.groups(count=carbom.18S.V4.unique.abund.filter.precluster.pick.count_table)"
# Total seqs: 2889838.
# Output File Names:
# carbom.18S.V4.unique.abund.filter.precluster.pick.count.summary

################################################################################
## 8 Split sequences between controls and samples
################################################################################
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/18S

# samples
mothur "#get.groups(fasta=carbom.18S.V4.unique.abund.filter.precluster.pick.fasta,count=carbom.18S.V4.unique.abund.filter.precluster.pick.count_table, accnos=carbom_groups_samples.txt)"
# Selected 5395 sequences from your fasta file.
# Selected 2956531 sequences from your count file.
# Output File names:
# carbom.18S.V4.unique.abund.filter.precluster.pick.pick.fasta
# carbom.18S.V4.unique.abund.filter.precluster.pick.pick.count_table

mv carbom.18S.V4.unique.abund.filter.precluster.pick.pick.fasta carbom.18S.V4.unique.abund.filter.precluster.pick.samples.fasta
mv carbom.18S.V4.unique.abund.filter.precluster.pick.pick.count_table carbom.18S.V4.unique.abund.filter.precluster.pick.samples.count_table
mothur "#count.groups(count=carbom.18S.V4.unique.abund.filter.precluster.pick.samples.count_table)"

# controls
mothur "#get.groups(fasta=carbom.18S.V4.unique.abund.filter.precluster.pick.fasta,count=carbom.18S.V4.unique.abund.filter.precluster.pick.count_table, accnos=cce_groups_controls.txt)"
# Output File names:
# carbom.18S.V4.unique.abund.filter.precluster.pick.pick.fasta
# carbom.18S.V4.unique.abund.filter.precluster.pick.pick.count_table

mv carbom.18S.V4.unique.abund.filter.precluster.pick.pick.fasta carbom.18S.V4.precluster.no_chimeras.controls.fasta
mv carbom.18S.V4.unique.abund.filter.precluster.pick.pick.count_table carbom.18S.V4.precluster.no_chimeras.controls.count_table

################################################################################
## 9  Remove low abundances sequences from samples - in our case we remove sequences with less than < 10  
################################################################################
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/18S

mothur "#split.abund(count=carbom.18S.V4.unique.abund.filter.precluster.pick.samples.count_table, fasta=carbom.18S.V4.unique.abund.filter.precluster.pick.samples.fasta, cutoff=10, accnos=true)"
# Output File Names:
# carbom.18S.V4.unique.abund.filter.precluster.pick.samples.rare.count_table
# carbom.18S.V4.unique.abund.filter.precluster.pick.samples.abund.count_table
# rare.accnos
# abund.accnos
# carbom.18S.V4.unique.abund.filter.precluster.pick.samples.rare.fasta
# carbom.18S.V4.unique.abund.filter.precluster.pick.samples.abund.fasta

################################################################################
## 10 - Remove short sequences  < 200 bp
################################################################################

mothur "#screen.seqs(fasta=carbom.18S.V4.unique.abund.filter.precluster.pick.samples.abund.fasta, count=carbom.18S.V4.unique.abund.filter.precluster.pick.samples.abund.count_table, minlength=200, processors=32)"
# Output File Names:
# carbom.18S.V4.unique.abund.filter.precluster.pick.samples.abund.good.fasta
# carbom.18S.V4.unique.abund.filter.precluster.pick.samples.abund.bad.accnos
# carbom.18S.V4.unique.abund.filter.precluster.pick.samples.abund.good.count_table

cp carbom.18S.V4.unique.abund.filter.precluster.pick.samples.abund.good.fasta        carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.fasta
cp carbom.18S.V4.unique.abund.filter.precluster.pick.samples.abund.good.count_table carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.count_table

mothur "#summary.seqs(fasta=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.fasta)"

# Output File Names:
# carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.summary
# It took 0 secs to summarize 871 sequences.

# recount number of sequences in each sample after removing low abundance sequences
mothur "#count.groups(count=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.count_table)"
# Total seqs: 2802466.
# Output File Names:
# carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.count.summary

################################################################################
## 11. Classify sequences
################################################################################
cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/18S

# Classify Controls - NOT DONE
mothur "#classify.seqs(fasta=carbom.18S.V4.precluster.no_chimeras.controls.fasta,count=carbom.18S.V4.precluster.no_chimeras.controls.count_table, reference=/projet/umr7144/dipo/vaulot/mothur/database/pr2_gb203_v_4.4/pr2_gb203_version_4.4.fasta, taxonomy=/projet/umr7144/dipo/vaulot/mothur/database/pr2_gb203_v_4.4/pr2_gb203_version_4.4.taxo, processors=32, probs=T)"

# Output File Names:
# carbom.18S.V4.precluster.no_chimeras.controls.4.wang.taxonomy
# carbom.18S.V4.precluster.no_chimeras.controls.4.wang.tax.summary

# Classify Samples
mothur "#classify.seqs(fasta=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.fasta, count=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.count_table, reference=/projet/umr7144/dipo/vaulot/mothur/database/pr2_gb203_v_4.4/pr2_gb203_version_4.4.fasta, taxonomy=/projet/umr7144/dipo/vaulot/mothur/database/pr2_gb203_v_4.4/pr2_gb203_version_4.4.taxo, processors=32, probs=T)"
# It took 132 secs to classify 871 sequences.
# Output File Names:
# carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.4.wang.taxonomy
# carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.4.wang.tax.summary


################################################################################
##11 Cluster and do OTUs - ONLY on SAMPLES and only on abundant sequences (nseq > 10)
################################################################################

mothur "#dist.seqs(fasta=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.fasta, processors=32)"
# Output File Names:
# carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.dist
# It took 7 seconds to calculate the distances for 871 sequences.

mothur "#cluster(column=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.dist, count=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.count_table)"
# It took 38 seconds to cluster
# Output File Names:
# carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.list

cd /projet/umr7144/dipo/cgerikasribeiro/mothur/carbom/18S
mothur "#classify.otu(taxonomy=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.4.wang.taxonomy, count=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.count_table, list=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.list, label=0.02, probs=F, basis=sequence)"
# 0.02  287  
# Output File Names:
# carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.0.02.cons.taxonomy
# carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.0.02.cons.tax.summary

mothur "#get.oturep(fasta=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.fasta, column=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.dist, count=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.count_table, list=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.list, cutoff=0.01)"
# Output File Names:
# carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.0.02.rep.count_table
# carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.0.02.rep.fasta


mothur "#create.database(list=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.list, count=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.0.02.rep.count_table, label=0.02, repfasta=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.0.02.rep.fasta , constaxonomy=carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.0.02.cons.taxonomy)"
# Output File Names:
# carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.database
mv carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.database carbom.18S.V4.uniq.preclust.no_chim.more_than_10.samples.an.unique_list.0.02.database


