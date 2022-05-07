###Download the raw sequencing data

for i in SRR*
do
echo $i
/home/xiaomin/16S-meta/sratoolkit.2.11.0-ubuntu64/bin/fastq-dump --gzip --split-3 $i
done

####qiime2-2020.6
###generate manifest.txt
cd /home/xiaomin/16S-meta/PRJNA767939/raw
echo "sample-id,forward-absolute-filepath" >> manifest1.csv
for i in `ls | awk '{print$0}'`;
do
ls $i | grep 1.fastq |awk -v a=$i -v b=$PWD '{print a","b"/"a""}'>> manifest1.csv;
done

echo "reverse-absolute-filepath" >> manifest2.csv
for i in `ls | awk '{print$0}'`;
do
ls $i | grep 2.fastq |awk -v a=$i -v b=$PWD '{print ""b"/"a""}'>> manifest2.csv;
done

paste -d, manifest1.csv manifest2.csv >> manifest.csv
rm manifest1.csv manifest2.csv

###importing data
conda activate qiime2-2020.6
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2 --input-path manifest.txt --output-path demux.qza 
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv

###Based on the plots in demux.qzv, choose the values for --p-trunc-len and --p-trim-left in this case

qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza --p-trunc-len-f 275 --p-trunc-len-r 227 --p-trim-left-f 23 --p-trim-left-r 27 --p-max-ee-f 2 --p-max-ee-r 4 --o-representative-sequences rep-seqs.qza --o-table table.qza --o-denoising-stats stats.qza
qiime metadata tabulate --m-input-file stats.qza --o-visualization stats.qzv
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file mapping.txt
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

###Taxonomic analysis
qiime feature-classifier classify-sklearn --i-classifier ~/16S-meta/classifier_gg_13_8_99.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

###Filter out Choloroplast and Mitochondira 
qiime taxa filter-table --i-table table.qza --i-taxonomy taxonomy.qza --p-exclude mitochondria,chloroplast,Archaea,Unassigned --o-filtered-table table-no-mitochondria-no-chloroplast.qza
qiime taxa filter-seqs --i-sequences rep-seqs.qza --i-taxonomy taxonomy.qza --p-exclude mitochondria,chloroplast,Archaea,Unassigned --o-filtered-sequences rep-seqs-no-mitochondria-no-chloroplast.qza
mv table-no-mitochondria-no-chloroplast.qza table.qza
mv rep-seqs-no-mitochondria-no-chloroplast.qza rep-seqs.qza

###
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file mapping.txt
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file mapping.txt --o-visualization taxa-bar-plots.qzv

###Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

###Alpha and beta diversity analysis
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth ${minnum} --m-metadata-file mapping.txt --output-dir core-metrics-results
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file mapping.txt --o-visualization core-metrics-results/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/evenness_vector.qza --m-metadata-file mapping.txt --o-visualization core-metrics-results/evenness-group-significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/shannon_vector.qza --m-metadata-file mapping.txt --o-visualization core-metrics-results/shannon-group-significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/observed_features_vector.qza --m-metadata-file mapping.txt --o-visualization core-metrics-results/observed_features-group-significance.qzv

for category_1 in Group; do echo $category_1; qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file mapping.txt --p-method permanova --m-metadata-column $category_1 --o-visualization core-metrics-results/unweighted_unifrac-permanova-${category_1}-significance.qzv--p-pairwise; qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza --m-metadata-file mapping.txt --p-method permanova --m-metadata-column $category_1 --o-visualization core-metrics-results/weighted_unifrac-permanova-${category_1}-significance.qzv--p-pairwise; qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza --m-metadata-file mapping.txt --p-method permanova --m-metadata-column $category_1 --o-visualization core-metrics-results/bray_curtis-permanova-${category_1}-significance.qzv--p-pairwise; done

####Exporting a feature table
for f in rep-seqs.qza table.qza taxonomy.qza ; do echo $f; qiime tools export --input-path $f --output-path exported; done
biom add-metadata -i exported/feature-table.biom -o exported/feature-table.taxonomy.biom --observation-metadata-fp exported/taxonomy.tsv --observation-header OTUID,taxonomy,confidence
biom convert -i exported/feature-table.taxonomy.biom -o exported/feature-table.taxonomy.txt --to-tsv --header-key taxonomy

###relative abundance
perl /home/xiaomin/16S-meta/github/Bayegy/stat_otu_tab.pl -unif min exported/feature-table.taxonomy.txt -prefix exported/Relative/otu_table --even exported/Relative/otu_table.even.txt -spestat exported/Relative/classified_stat_relative.xls
mv exported/Relative/otu_table.p.relative.mat exported/Relative/otu_table.Phylum.relative.txt
mv exported/Relative/otu_table.c.relative.mat exported/Relative/otu_table.Class.relative.txt
mv exported/Relative/otu_table.o.relative.mat exported/Relative/otu_table.Order.relative.txt
mv exported/Relative/otu_table.f.relative.mat exported/Relative/otu_table.Family.relative.txt
mv exported/Relative/otu_table.g.relative.mat exported/Relative/otu_table.Genus.relative.txt
mv exported/Relative/otu_table.s.relative.mat exported/Relative/otu_table.Species.relative.txt

###picrust2 
conda activate picrust2
cd exported
picrust2_pipeline.py -s dna-sequences.fasta -i feature-table.biom -o picrust2_out_pipeline -p 1


######################################meta-analysis of alpha diversity
setwd("E:/Programming/R/R_CHXM/16S_meta_manuscript/alpha_diversity/HC_CP")
###install.packages("meta")
###install.packages("metafor")
library(meta)
library(metafor)

###Pielou’s evenness
test_data<-read.table("alpha_HC_CP_evenness.txt", header = TRUE,sep="\t")
head(test_data)
# meta-analysis with continuout outcome
# comb.fixed/comb.random: indicator whether a fix/random effect mata-analysis to be conducted.
# sm: Three different types of summary measures to choose,standardized mean difference (SMD),mean difference (MD), ratio of means (ROM)
res.flesiss =  metacont(Ne, Me, Se, Nc, Mc, Sc,
                        fixed = T, random = T, studlab=paste(author, year), 
                        data = test_data, sm = "SMD") 
res.flesiss

forest (res.flesiss, leftcols = c('studlab', 'Ne', 'Nc'), family="sans" , fontsize=9.5,
        lwd=2,col.diamond.fixed="#20854EFF",
        col.diamond.lines.fixed="black" ,
        col.diamond.random="#BC3C29FF" , col.diamond.lines.random="black" ,col.square="skyblue" , col.study="lightslategray" ,
        lty.fixed=4,plotwidth="6cm" , colgap.forest.left="0.5cm",
        colgap.forest.right="0.5cm" , just.forest="right", colgap.left="0.5cm",colgap.right="0.5cm")

###funnel plot, Begg’s correlation test and Egger’s regression test 
funnel(res.flesiss)
metabias(res.flesiss,method.bias = "rank")  # begg
metabias(res.flesiss,method.bias = "linreg") # egger

###Sensitivity analysis that removal of either of the included studies
m2=metainf(res.flesiss,pooled = "random")
m2
forest(m2, leftcols = c('studlab','I2'),col.square = "skyblue",col.diamond = "#BC3C29FF" )

###Faith’s phylogenetic diversity
rm(list=ls()) # clear all
test_data<-read.table("alpha_HC_CP_FD.txt", header = TRUE,sep="\t")
head(test_data)
# meta-analysis with continuout outcome
# comb.fixed/comb.random: indicator whether a fix/random effect mata-analysis to be conducted.
# sm: Three different types of summary measures to choose,standardized mean difference (SMD),mean difference (MD), ratio of means (ROM)
res.flesiss =  metacont(Ne, Me, Se, Nc, Mc, Sc,
                        fixed = T, random = T, studlab=paste(author, year),
                        data = test_data, sm = "SMD") 
res.flesiss

forest (res.flesiss, leftcols = c('studlab', 'Ne', 'Nc'), family="sans" , fontsize=9.5,
        lwd=2,col.diamond.fixed="#20854EFF",
        col.diamond.lines.fixed="black" ,
        col.diamond.random="#BC3C29FF" , col.diamond.lines.random="black" ,col.square="skyblue" , col.study="lightslategray" ,
        lty.fixed=4,plotwidth="6cm" , colgap.forest.left="0.5cm",
        colgap.forest.right="0.5cm" , just.forest="right", colgap.left="0.5cm",colgap.right="0.5cm")

### funnel plot, Begg’s correlation test and Egger’s regression test
funnel(res.flesiss)
metabias(res.flesiss,method.bias = "rank")  # begg
metabias(res.flesiss,method.bias = "linreg") # egger

### Sensitivity analysis that removal of either of the included studies
m2=metainf(res.flesiss,pooled = "random")
m2
forest(m2, leftcols = c('studlab','I2'),col.square = "skyblue",col.diamond = "#BC3C29FF" )

### observed species
rm(list=ls()) # clear all
test_data<-read.table("alpha_HC_CP_observed.txt", header = TRUE,sep="\t")
head(test_data)
# meta-analysis with continuout outcome
# comb.fixed/comb.random: indicator whether a fix/random effect mata-analysis to be conducted.
# sm: Three different types of summary measures to choose,standardized mean difference (SMD),mean difference (MD), ratio of means (ROM)
res.flesiss =  metacont(Ne, Me, Se, Nc, Mc, Sc,
                        fixed = T, random = T, studlab=paste(author, year),
                        data = test_data, sm = "SMD") 
res.flesiss

forest (res.flesiss, leftcols = c('studlab', 'Ne', 'Nc'), family="sans" , fontsize=9.5,
        lwd=2,col.diamond.fixed="#20854EFF",
        col.diamond.lines.fixed="black" ,
        col.diamond.random="#BC3C29FF" , col.diamond.lines.random="black" ,col.square="skyblue" , col.study="lightslategray" ,
        lty.fixed=4,plotwidth="6cm" , colgap.forest.left="0.5cm",
        colgap.forest.right="0.5cm" , just.forest="right", colgap.left="0.5cm",colgap.right="0.5cm")

###funnel plot, Begg’s correlation test and Egger’s regression test
funnel(res.flesiss)
metabias(res.flesiss,method.bias = "rank")  # begg
metabias(res.flesiss,method.bias = "linreg") # egger

###Sensitivity analysis that removal of either of the included studies

m2=metainf(res.flesiss,pooled = "random")
m2
forest(m2, leftcols = c('studlab','I2'),col.square = "skyblue",col.diamond = "#BC3C29FF" )

### Shannon’s diversity index
rm(list=ls()) # clear all
test_data<-read.table("alpha_HC_CP_shannon.txt", header = TRUE,sep="\t")
head(test_data)
# meta-analysis with continuout outcome
# comb.fixed/comb.random: indicator whether a fix/random effect mata-analysis to be conducted.
# sm: Three different types of summary measures to choose,standardized mean difference (SMD),mean difference (MD), ratio of means (ROM)
res.flesiss =  metacont(Ne, Me, Se, Nc, Mc, Sc,
                        fixed = T, random = T, studlab=paste(author, year),
                        data = test_data, sm = "SMD") 
res.flesiss

forest (res.flesiss, leftcols = c('studlab', 'Ne', 'Nc'), family="sans" , fontsize=9.5,
        lwd=2,col.diamond.fixed="#20854EFF",
        col.diamond.lines.fixed="black" ,
        col.diamond.random="#BC3C29FF" , col.diamond.lines.random="black" ,col.square="skyblue" , col.study="lightslategray" ,
        lty.fixed=4,plotwidth="6cm" , colgap.forest.left="0.5cm",
        colgap.forest.right="0.5cm" , just.forest="right", colgap.left="0.5cm",colgap.right="0.5cm")

### funnel plot, Begg’s correlation test and Egger’s regression test
funnel(res.flesiss)
metabias(res.flesiss,method.bias = "rank")  # begg
metabias(res.flesiss,method.bias = "linreg") # egger

### Sensitivity analysis that removal of either of the included studies
m2=metainf(res.flesiss,pooled = "random")
m2
forest(m2, leftcols = c('studlab','I2'),col.square = "skyblue",col.diamond = "#BC3C29FF")

######################################meta-analysis of microbial composition
# BiocManager::install("metamicrobiomeR")
setwd("E:/Programming/R/R_CHXM/16S_meta_manuscript/taxa/genus/CP_HC/data_HC_CP/all")
library(metamicrobiomeR)
library(data.table)
library(stringr)
library(plyr)

rm(list=ls()) # clear all
###genus
####Directory of data
filedir <- "E:/Programming/R/R_CHXM/16S_meta_manuscript/taxa/genus/CP_HC/data_HC_CP/all"

myfile <- list.files(filedir)
myfile <- str_replace_all(myfile,".txt","")
taxacompare_all <- c()
for (j in myfile) {
  mydata <- read.csv(paste(filedir,"/",j,".txt",sep = ""),header = TRUE,row.names = 1,sep = "\t")
  personid <- as.factor(row.names(mydata))
  x.sampleid <- as.factor(row.names(mydata))
  Group <- mydata$Group
  mydata <- mydata[,-1]
  newdata <- as.data.frame(matrix(data = 0,nrow = nrow(mydata),ncol = ncol(mydata)))
  names(newdata) <- names(mydata)
  for (i in 1:ncol(mydata)) {
    newdata[,i] <- as.numeric(mydata[,i])
  }
  newdata$personid <- personid
  newdata$Group <- Group
  newdata$x.sampleid <- x.sampleid
  names(newdata) <- str_replace_all(names(newdata),"[|]",".")
  names(newdata) <- str_replace_all(names(newdata),"[.]$","")
  taxacompare <- taxa.compare(taxtab=newdata,propmed.rel="gamlss",comvar="Group",adjustvar="Group",longitudinal="no",p.adjust.method="fdr", percent.filter = 0.05, relabund.filter = 5e-05)
  taxacompare$study <- j
  taxacompare$pop <- j
  taxacompare_all <- rbind.fill(taxacompare_all,taxacompare)
}


metataxaaa <- meta.taxa(taxcomdat=taxacompare_all, summary.measure="RR", pool.var="id", studylab="study", backtransform=FALSE, percent.meta=0.5, p.adjust.method="fdr")

###table
metataxa <- metatab.show(metatab=metataxaaa$random,com.pooled.tab=taxacompare_all, highest.lev = "g",
                         tax.lev="l6",showvar=names(metataxaaa$random),p.cutoff.type="p", p.cutoff=0.05,display="table")
write.table (metataxa, file ="metataxa.csv", sep =",", row.names =FALSE)


###heatmap
metataxa <- metatab.show(metatab=metataxaaa$random,com.pooled.tab=taxacompare_all,
                         tax.lev="l6",showvar=names(metataxaaa$random),p.cutoff.type="p", p.cutoff=0.05,display="data",plot="heatmap")


meta.niceplot(metadat=metataxa,sumtype="taxa",level="sub",p="p",
              p.adjust="p.adjust",phyla.col="black",p.sig.heat="no",leg.key.size=0.4,leg.text.size=5,heat.text.x.size=8,heat.text.x.angle=45,
              forest.axis.text.y=7,forest.axis.text.x=10,heat.forest.width.ratio=c(1,1.2),point.ratio=c(3,1.5),line.ratio=c(2,1))


dev.copy(device = x11)