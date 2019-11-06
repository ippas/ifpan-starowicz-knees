require(edgeR)
require(preprocessCore)
require(magrittr)
require(gplots)
require(RColorBrewer)
require(tidyverse)

#create sample info
samples <- data.frame(read.table('/home/ifpan/projects/ifpan-starowicz-knees/cuffnorm/samples.table',header = TRUE, sep = "\t"))[,c(1,2)]
samples$file <- substr(samples$file,13,17)
samples <- samples %>% mutate(file = str_replace(file, "/", ""))

groups <- data.frame(read.table('/home/ifpan/projects/ifpan-starowicz-knees/analysis/groups.csv',header = TRUE, sep = ";"))
groups <- rbind.data.frame(groups, groups)
groups$side <- c(rep("Cai", 36),rep("Cac", 36))
groups$sample <- paste(groups$side,groups$Numer, sep="")

samples$side <- groups$side[match(samples$file, groups$sample)]
samples$drug <- groups$treatment[match(samples$file, groups$sample)]
samples$damage <- groups$type[match(samples$file, groups$sample)]
samples$animal <- groups$Numer[match(samples$file, groups$sample)]

#changed mistake in naming of groups
samples$damage[samples$animal == 33 | samples$animal == 34 | samples$animal == 35 | samples$animal == 36] <- "NaCl"

#change ordering and fix drug info
samples <- samples[order(samples$animal),]
samples$drug <- factor(c(rep("intact", 12), rep("jwh", 11), rep("veh", 9), rep("intact", 11), rep("jwh", 12), rep("veh", 12)))

samples$group <- paste(samples$side, samples$damage, samples$drug, sep=".")


#remove incomplete samples:
samples.paired <- samples[samples$animal != 11 & samples$animal != 15 &  samples$animal != 17 & samples$animal != 18 & samples$animal != 24,]

rm(groups)


#get data and rat genes annotation
fpkm <- data.frame(read.table('/home/ifpan/projects/ifpan-starowicz-knees/cuffnorm/genes.fpkm_table',header = TRUE, sep = "\t"))
download.file('http://149.156.177.112/projects/ifpan-annaradli-ldopa/rproject/mart_export_rnor.txt','rnor.genes.tsv')
rnor.gene.list <- data.frame(read.delim("rnor.genes.tsv"))
fpkm$gene.name <- rnor.gene.list$Gene.name[match(fpkm$tracking_id, rnor.gene.list$Gene.stable.ID)]
fpkm <- data.frame(fpkm$gene.name, fpkm[,c(1:68)])
rownames(fpkm) <- fpkm$tracking_id
colnames(fpkm) <- c("gene.name","transcrip.ID",samples$file[match(colnames(fpkm[,c(3:69)]), samples$sample_id)])
fpkm <- data.frame(fpkm[,c(1,2)], fpkm[,c(3:69)][,samples$file])

rm(rnor.gene.list)

#normalize the distribution 
fpkms.normalised <- data.matrix(fpkm[,c(-1,-2)])
normalize.quantiles(fpkms.normalised,copy=FALSE)
fpkms.log <- log2(fpkms.normalised + 1)

rm(fpkms.normalised)

#remove fpkms.log that have rowMeans < 1
fpkms.log[rowMeans(fpkms.log) < 1,] <- NA

#exploratory analysis:

stat.one.way <- function(counts, groups) {
  ifelse(
    is.na(counts[1]),
    NA,
    unlist(summary(aov(counts ~ as.factor(groups))))[9]) }


results <- data.frame(fpkm$gene.name, fpkm$transcrip.ID, row.names = rownames(fpkm))

#perform one-way stat on all groups:

apply(fpkms.log[,match(samples$file, colnames(fpkms.log))],
      1,
      stat.one.way,
      groups=samples$group) %>% 
  p.adjust(method="fdr") -> results$p.one.way

#three-way on all paired samples:

stat <- function(counts,
                 side,
                 damage,
                 drug,
                 animal) {
  if (is.na(counts[1])) {
    c(rep(NA, 7))
  } else {
    unlist(summary(aov(counts ~ 
                         as.factor(side)
                       * as.factor(damage)
                       * as.factor(drug) + 
                         Error(as.factor(animal)/as.factor(side)))))[c(17,18,19,41,42,43,44)]
                         }
                          }

apply(fpkms.log[,match(samples.paired$file, colnames(fpkms.log))],
      1,
      stat,
      side = samples.paired$side,
      damage = samples.paired$damage,
      drug = samples.paired$drug,
      animal = samples.paired$animal
      ) %>% t %>% apply(.,
        2,
        p.adjust,
        method="fdr") %>% 
  data.frame() %>%
  bind_cols(results,.) -> results


colnames(results) <- c("gene.name", "transcipt.ID", "p.one.way", "p.damage", "p.drug", "p.drug.damage", "p.side", "p.side.damage", "p.side.drug", "p.side.drug.damage")

results <- bind_cols(results, data.frame(fpkms.log))

write.csv(results, "all-genes-aov-adj.csv", row.names = FALSE)
write.csv(samples, "samples-all.csv", row.names = FALSE)


#plot heatmap of top damage induced-genes:

mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)

to.plot <- fpkms.log[match(
  results$transcipt.ID[which(results$p.damage < 0.00003)],
  rownames(fpkms.log)),
  samples.paired$file[order(samples.paired$group)]]

group.names <- unique(samples.paired$group[order(samples.paired$group)])

col.labels <- c(rep("", 3), group.names[1], rep(" ", 4), 
                group.names[2], rep(" ", 5),
                group.names[3], rep(" ", 5),
                group.names[4], rep(" ", 4),
                group.names[5], rep(" ", 4),
                group.names[6], rep(" ", 3),
                group.names[7], rep(" ", 4),
                group.names[8], rep(" ", 4),
                group.names[9], rep(" ", 5),
                group.names[10],rep(" ", 5),
                group.names[11],rep(" ", 4),
                group.names[12])



cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}

to.plot %>% 
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 3) %>%
  t %>%
  `colnames<-`(colnames(to.plot)) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    Colv = FALSE,
    main="",
    scale="row",
    colsep = c(6,11,14,19,25,31,37,42,45,50,56),
    sepwidth = c(0.3,0.3),
    labRow=fpkm$gene.name[match(rownames(.), rownames(fpkm))],
    labCol=col.labels,         
    srtCol = 45,
    cexRow = 0.3,
    offsetCol = 0.1
  )

#plot a heatmap of top drug-induced genes:

to.plot <- fpkms.log[match(
  results$transcipt.ID[which(results$p.drug < 0.192)],
  rownames(fpkms.log)),
  samples.paired$file[order(samples.paired$group)]]


to.plot %>% 
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 3) %>%
  t %>%
  `colnames<-`(colnames(to.plot)) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    Colv = FALSE,
    main="",
    scale="row",
    colsep = c(6,11,14,19,25,31,37,42,45,50,56),
    sepwidth = c(0.3,0.3),
    labRow=fpkm$gene.name[match(rownames(.), rownames(fpkm))],
    labCol=col.labels,         
    srtCol = 45,
    cexRow = 0.3,
    offsetCol = 0.1
  )

#calculate fold-changes:

fold.change <- function(x, group, ctrl) {
  abs(mean(as.numeric((x[match(samples$file[samples$group == ctrl],
                                      colnames(results))])))
    -
      mean(as.numeric((x[match(samples$file[samples$group == group],
                                      colnames(results))]))))
}
  
results %>% mutate(
  fold.change.jwh.vs.veh.ipsi=
  apply(results, 1, fold.change, group = "Cai.MIA.jwh", ctrl = "Cai.MIA.veh"),
  fold.change.jwh.vs.veh.contra=
    apply(results, 1, fold.change, group = "Cac.MIA.jwh", ctrl = "Cac.MIA.veh")
) -> results
  
# find genes that are significantly affected by dmaage and have a <1.3 fold change
# and plot them

results %>%
  filter(p.damage < 0.05) %>%
  filter(fold.change.jwh.vs.veh.ipsi > 1.3 |
         fold.change.jwh.vs.veh.contra > 1.3
         ) %>% data.frame() -> results.fold.filtered

rownames(results.fold.filtered) <- results.fold.filtered$gene.name


write.csv(results.fold.filtered, "filtered-results-fold.csv", row.names = FALSE)
write.csv(results, "results-all-fold.csv", row.names = FALSE)

#results.fold.filtered$gene.name %>% 
#write.table(quote = FALSE, row.names = FALSE, sep = "\t")


results.fold.filtered[,samples.paired$file[order(samples.paired$group)]] %>%
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 3) %>%
  t %>%
  `colnames<-`(colnames(to.plot)) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    Colv = FALSE,
    main="",
    scale="row",
    colsep = c(6,11,14,19,25,31,37,42,45,50,56),
    sepwidth = c(0.3,0.3),
    labCol=col.labels,         
    srtCol = 45
  )

#export summary table:

results %>% mutate(
  mean.exp = rowMeans(results[,samples$file])
) -> results

results.summary <- results[-c(3,11:77)]

colnames(results.summary) <- c(
  "gene.name", "transcript.ID", 
  "p.type", "p.treatment", 
  "p.type.treatment", "p.side", 
  "p.side.type", "p.side.treatment",
  "p.side.type.treatment", "fold.mia.jwh.vs.veh.ipsi",
  "fold.mia.jwh.vs.veh.contra", "mean.exp"
)

write.csv(results.summary, "knees-results-summary.csv", row.names = FALSE)
