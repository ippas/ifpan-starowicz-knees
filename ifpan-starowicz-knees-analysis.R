require(edgeR)
require(preprocessCore)
require(magrittr)
require(gplots)
require(RColorBrewer)
require(tidyverse)
require(stringi)

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

#write.csv(results, "all-genes-aov-adj.csv", row.names = FALSE)
#write.csv(samples, "samples-all.csv", row.names = FALSE)

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
    apply(results, 1, fold.change, group = "Cac.MIA.jwh", ctrl = "Cac.MIA.veh"),
  fold.change.mia.vs.nacl.veh.ipsi=
    apply(results, 1, fold.change, group = "Cai.NaCl.veh", ctrl = "Cai.MIA.veh"),
  fold.change.mia.vs.nacl.intact.ipsi=
    apply(results, 1, fold.change, group = "Cai.NaCl.intact", ctrl = "Cai.MIA.intact"),
  fold.change.intact.vs.jwh.mia.ipsi=
    apply(results, 1, fold.change, group = "Cai.MIA.jwh", ctrl = "Cai.MIA.intact"),
  fold.change.veh.vs.jwh.mia.ipsi=
    apply(results, 1, fold.change, group = "Cai.MIA.jwh", ctrl = "Cai.MIA.veh")
) -> results


write.csv(results, "results-all-fold.csv", row.names = FALSE)


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

results.summary %>% mutate(
  fold.change.mia.vs.nacl.veh.ipsi = results$fold.change.mia.vs.nacl.veh.ipsi,
  fold.change.mia.vs.nacl.intact.ipsi= results$fold.change.mia.vs.nacl.intact.ipsi
) -> results.summary

write.csv(results.summary, "knees-results-summary.csv", row.names = FALSE)

### adding mean exp from requested groups and creating a new table
### Requested groups included: 

selected_means <- results[,c(1,2)]

selected_means %>% mutate(
  mean.ipsi.mia.veh = rowMeans(results[,samples$file[samples$group == "Cai.MIA.veh"]]),
  mean.ipsi.mia.jwh = rowMeans(results[,samples$file[samples$group == "Cai.MIA.jwh"]]),
  mean.ipsi.nacl.veh = rowMeans(results[,samples$file[samples$group == "Cai.NaCl.veh"]]),
  mean.ipsi.nacl.jwh = rowMeans(results[,samples$file[samples$group == "Cai.NaCl.jwh"]])
) -> selected_means

selected_means <- na.omit(selected_means)
write.csv(selected_means, "selected_means.csv", row.names = FALSE)

### calculate a t-test for Cai.MIA.veh vs Cai.MIA.jwh

a <- as.factor(samples$group[which(samples$group == "Cai.MIA.veh" | 
                                     samples$group == "Cai.MIA.jwh")])


stat.paired.t <- function(x) {
  if (is.na(x[1])) { NA
  } else {
    tryCatch((t.test(x ~ a))$p.value, error=function(err) NA)
  }}



fpkms.log[,samples$file[which(samples$group == "Cai.MIA.veh" | 
                                samples$group == "Cai.MIA.jwh")]] %>% 
  apply(1, stat.paired.t) -> results$t.test.mia.cai.veh.vs.jwh

write.csv(results, "all-genes-with-jwh-t-test-unadjusted.csv", row.names = FALSE)


#gen a list of genes provided by Kuba and give means for individual samples

up.from.kuba <- scan('up-from-kuba.txt', character())
down.from.kuba <- scan('down-from-kuba.txt', character())
changed.from.kuba <- scan('changed-from-kuba.txt', character())

results %>%
  filter(gene.name %in% changed.from.kuba) %>%
  select(gene.name, transcipt.ID, samples$file[which(samples$group == "Cai.MIA.veh" | 
                                                       samples$group == "Cai.MIA.jwh" |
                                                       samples$group == "Cai.NaCl.veh" |
                                                       samples$group == "Cai.NaCl.jwh")]) -> selected.fpkm.logs

write.csv(selected.fpkm.logs, "selected-fpkms-log.csv", row.names = FALSE)

# identify responders based on preprepared list of genes (from this list https://link.springer.com/article/10.1007/s10142-017-0576-6/figures/3?shared-article-renderer)
pre_genes = c("Pck1", "Adipoq", "Fabp4", "Lipe", "Itga8", "Itgb8", "Lamb2", "Prdk1",
              "Prkca", "Fzd8", "Fgf2", "Agt", "Frzb", "Ace", "Bmp6", "Bmp4", "Mmp14",
              "Mmp2", "Inhba", "Tlr7", "Oplah", "Ahcy", "Tlr2", "Tlr1")

results %>% filter(results$gene.name %in% pre_genes) %>% filter(!is.na(p.one.way)) -> pre_results

#Samples 

to.remove = c("Cai4", "Cai6", "Cai8", "Cai10")
samples.no.outliers <- samples[!(samples$file %in% to.remove),]

results.no.outliers <- results[c(1,2)]

#perform one-way stat on all groups without outliers:

apply(fpkms.log[,match(samples.no.outliers$file, colnames(fpkms.log))],
      1,
      stat.one.way,
      groups=samples.no.outliers$group) %>% 
  p.adjust(method="fdr") -> results.no.outliers$p.one.way

### calculate a t-test for Cai.MIA.veh vs Cai.MIA.jwh

a <- as.factor(samples.no.outliers$group[which(samples.no.outliers$group == "Cai.MIA.veh" | 
                                     samples.no.outliers$group == "Cai.MIA.jwh")])


fpkms.log[,samples.no.outliers$file[which(samples.no.outliers$group == "Cai.MIA.veh" | 
                                samples.no.outliers$group == "Cai.MIA.jwh")]] %>% 
  apply(1, stat.paired.t) -> results.no.outliers$t.test.mia.cai.veh.vs.jwh

## add counts

results.no.outliers <- bind_cols(results.no.outliers, data.frame(fpkms.log))

## add fold change 

fold.change.2 <- function(x, group, ctrl) {
  abs(mean(as.numeric((x[match(samples.no.outliers$file[samples.no.outliers$group == ctrl],
                               colnames(results.no.outliers))])))
      -
        mean(as.numeric((x[match(samples.no.outliers$file[samples.no.outliers$group == group],
                                 colnames(results.no.outliers))]))))
}

results.no.outliers %>% mutate(fold.change.jwh.vs.veh.ipsi=
  apply(results.no.outliers, 1, fold.change.2, group = "Cai.MIA.jwh", ctrl = "Cai.MIA.veh")) -> results.no.outliers


### HEATMAP PLOTTING:


mypalette <- brewer.pal(11,"RdBu")
morecols <- colorRampPalette(mypalette)

group.names <- unique(samples.paired$group[order(samples.paired$group)])

col.labels <- c(rep("", 3), group.names[1], rep(" ", 4), 
                group.names[2], rep(" ", 3),
                group.names[3], rep(" ", 4),
                group.names[4], rep(" ", 4),
                group.names[5], rep(" ", 5),
                group.names[6], rep(" ", 5),
                group.names[7], rep(" ", 5),
                group.names[8], rep(" ", 5),
                group.names[9], rep(" ", 5),
                group.names[10],rep(" ", 5),
                group.names[11],rep(" ", 5),
                group.names[12])


cut.threshold <- function(x, threshold = 2.5) {
  x[x > threshold] <- threshold
  x[x < -threshold] <- -threshold
  x
}


#plot heatmap of top damage induced-genes
to.plot <- fpkms.log[match(
  results$transcipt.ID[which(results$p.damage < 0.00003)],
  rownames(fpkms.log)),
  samples$file[order(samples$group)]]

#plot a heatmap of top drug-induced genes:
to.plot <- fpkms.log[match(
  results$transcipt.ID[which(results$p.drug < 0.192)],
  rownames(fpkms.log)),
  samples$file[order(samples$group)]]

#plot a heatmap of top drug-damage interaction genes:
to.plot <- fpkms.log[match(
  results$transcipt.ID[which(results$p.drug.damage < 0.3)],
  rownames(fpkms.log)),
  samples$file[order(samples$group)]]


#plot a heatmap of top drug-damage interaction genes:
to.plot <- fpkms.log[match(
  results$transcipt.ID[which(results$p.drug.damage < 0.3)],
  rownames(fpkms.log)),
  samples$file[order(samples$group)]]

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
    colsep = c(6,11,14,19,25,31,37,43,49,55,61),
    sepwidth = c(0.3,0.3),
    labRow=fpkm$gene.name[match(rownames(.), rownames(fpkm))],
    labCol=col.labels,         
    srtCol = 45,
    cexRow = 0.3,
    offsetCol = 0.1
  )

# custom filter go see jwh-reversed genes:
results %>%
  filter(p.one.way < 0.05) %>%
  filter(fold.change.jwh.vs.veh.ipsi > 0.6) %>% 
  filter(t.test.mia.cai.veh.vs.jwh < 0.1) %>%
  filter_at(samples$file, all_vars(. > 2)) -> to.plot

write.csv(to.plot, "selected-oneway-0-05-fold-0-6-ttest-0-1-all-counts-2.csv", row.names = FALSE)

# custom filter go see jwh-reversed genes:
results.no.outliers %>%
  filter(p.one.way < 0.1) %>%
  filter(fold.change.jwh.vs.veh.ipsi > 0.7) %>% 
  filter(t.test.mia.cai.veh.vs.jwh < 0.2) %>%
  filter_at(samples$file, all_vars(. > 2)) -> to.plot


to.plot %>%
  select(samples$file[order(samples$group)]) %>%
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 3) %>%
  t %>%
  `colnames<-`(samples$file[order(samples$group)]) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    Colv = FALSE,
    cexRow = 1,
    main="",
    scale="row",
    colsep = c(6,11,14,19,25,31,37,43,49,55,61),
    sepwidth = c(0.3,0.3),
    labCol=col.labels,         
    srtCol = 45,
    labRow = to.plot$gene.name,
    offsetCol = 0
  )


#plot heatmaps for genes from Kuba:

results %>%
  filter(gene.name %in% changed.from.kuba) -> to.plot.all

to.plot <- to.plot.all[complete.cases(to.plot),]

#plot heatmap of the PCR genes
pcr.genes <- c("Ccl2", "Col2a1", "Il6", "Timp1", "Timp4", "Comp","Mmp2", "Mmp3", "Mmp9", "Mmp13")

results %>%
  filter(gene.name %in% pcr.genes) -> to.plot.all
 
to.plot <- to.plot.all[complete.cases(to.plot),]


#Heatmap of preselected genes:

read.table("mia_selected_genes.csv") %>% as.vector() %>% t() -> preselected_genes
preselected_genes_groups <- read.csv("gene_groups.csv", sep = "\t")

results %>%
  filter(gene.name %in% preselected_genes) %>% select (-c(fold.change.jwh.vs.veh.ipsi, fold.change.jwh.vs.veh.contra)) -> results.preselected



#run the 3-way ANOVA again:

apply(results.preselected[,match(samples.paired$file, colnames(results.preselected))],
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
  data.frame() -> temp.results


colnames(temp.results) <- c("p.damage", "p.drug", "p.drug.damage", "p.side", "p.side.damage", "p.side.drug", "p.side.drug.damage")

results.preselected[,4:10] <- temp.results


write.csv(results.preselected, "results_preselected_genes.csv")


# add info about gene groups that each selected gene belongs to

check_gene <- function(x) { 
  if (x %in% preselected_genes_groups$MMPs) {return("MMPs")}
  if (x %in% preselected_genes_groups$Ils) {return("Ils")}
  if (x %in% preselected_genes_groups$CCls) {return("CCls")}
  if (x %in% preselected_genes_groups$Cols) {return("Cols")}
  if (x %in% preselected_genes_groups$ECS) {return("ECS")}
  }

results.preselected$gene_group <-lapply(results.preselected$gene.name, check_gene)


results.preselected %>% filter(gene_group == "ECS") %>% filter(p.damage < 0.05) %>% nrow() 

results.preselected %>% filter(gene_group == "ECS") %>% filter(t.test.mia.cai.veh.vs.jwh < 0.3) %>% nrow()

results.preselected %>%filter(gene_group == "ECS") %>%  filter(p.damage < 0.05) %>% filter(t.test.mia.cai.veh.vs.jwh < 0.3) %>% nrow()

results.preselected %>%
  filter(p.drug < 0.05) -> to.plot


 
 to.plot %>%
  select(samples$file[order(samples$group)]) %>%
  apply(1, scale) %>%
  t %>%
  apply(1, cut.threshold, threshold = 2.5) %>%
  t %>% `colnames<-`(samples$file[order(samples$group)]) %>%
  heatmap.2(
    distfun = function(x) as.dist(1-cor(t(x))),
    col=rev(morecols(50)),trace="none",
    Colv = FALSE,
    cexRow = 1,
    main="",
    scale="row",
    colsep = c(6,11,14,19,25,31,37,43,49,55,61),
    sepwidth = c(0.3,0.3),
    labCol=col.labels,         
    srtCol = 45,
    labRow = to.plot.all$gene.name,
    offsetCol = 0
  )

 
 
 
 
 
 
 
