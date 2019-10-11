# ifpan-starowicz-knees
### RNAseq analysis

#### STEP 1: Quality Control with fastqc v0.11.8

generate json input files based on the file_1.fq.gz file_2.fq.gz naming convention:

```
ls *1.fq.gz | xargs -i bash -c 'BASENAME2=$(echo {} | cut -d "." -f 1 | cut -d "_" -f 1); echo $BASENAME2' | xargs -i bash -c 'echo "{\"quality_check_fastqc_workflow.quality_check_fastqc.fastq_1\":\"{}_1.fq.gz\",\"quality_check_fastqc_workflow.quality_check_fastqc.fastq_2\":\"{}_2.fq.gz\"}">{}-input.json'
```

run qc on an example file:

```
java -jar /opt/tools/cromwell-44.jar run https://raw.githubusercontent.com/gosborcz/workflows/master/intelliseq-quality-check-fastqc.wdl -i Cac10-input.json > Cac10.txt
```

To run the workflow on all files in the folder (under screen):

```
ls *.json | xargs -i bash -c 'java -jar /opt/tools/cromwell-44.jar run https://gitlab.com/intelliseq/workflows/raw/master/src/main/wdl/tasks/quality-check-fastqc/v0.1/quality-check-fastqc.wdl -i {} > log-{}.txt'
```
To generate the MultiQC 1.7 report:

```
docker run --rm -v $PWD:/data ewels/multiqc:latest multiqc /data -o /data
```
[MultiQC report can be found under this link](http://149.156.177.112/projects/ifpan-starowicz-knees/fq/multiqc_report.html#fastqc)

#### STEP 2: Alignment to the rat genome with Hisat2

Generate input jsons:

```
ls *1.fq.gz | xargs -i bash -c 'BASENAME=$(echo {} | cut -d "." -f 1 | cut -d "_" -f 1); echo $BASENAME' | xargs -i bash -c 'echo "{\"align_to_rat_genome.align_with_hisat2.fastq1\":\"{}_1.fq.gz\",\"align_to_rat_genome.align_with_hisat2.sample_name\":\"{}\",\"align_to_rat_genome.align_with_hisat2.fastq2\":\"{}_2.fq.gz\"}">{}-input.json'
```
Run the workflow (Hisat2 2.1.0.)
```
ls *1.fq.gz | xargs -i bash -c 'BASENAME=$(echo {} | cut -d "." -f 1 | cut -d "_" -f 1); echo $BASENAME' | xargs -i bash -c 'java -jar /opt/tools/cromwell-44.jar run https://raw.githubusercontent.com/gosborcz/workflows/master/align-with-hisat2-to-rat-genome -i {}-input.json > log-{}.txt'
```
#### STEP 3: Transcript abundance estimation (cufflinks package v.2.2.1)

1. Create folders to store cuffquant output: 
```
ls ./bam | xargs -i basename {} .bam | xargs -i mkdir ./cuffquant/{}
```
2. Run cuffquant
```
ls ./bam | xargs -i basename {} .bam | xargs \
-i bash -c 'docker run -d --rm \
-v $PWD:/data octavianus90/cufflinks_final:latest \
cuffquant -o /data/cuffquant/{} \
/data/rn6/Rattus_norvegicus.Rnor_6.0.90.gtf /data/bam/{}.bam'
```
3. Run cuffnorm
```
docker run -rm -it -v $PWD:/data octavianus90/cufflinks_final:latest \
/bin/bash
  $ LIST=`ls /cuffquant/*/*`
  $ cuffnorm -o /data/cuffnorm /data/rn6/Rattus_norvegicus.Rnor_6.0.90.gtf $LIST
```
