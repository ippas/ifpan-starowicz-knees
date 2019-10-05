# ifpan-starowicz-knees
RNAseq analysis

STEP 1: Quality Control with fastqc VERSION!

generate json input files based on the file_1.fq.gz file_2.fq.gz naming convention:

```
ls *1.fq.gz | xargs -i bash -c 'BASENAME2=$(echo {} | cut -d "." -f 1 | cut -d "_" -f 1); echo $BASENAME2' | xargs -i bash -c 'echo "{\"quality_check_fastqc_workflow.quality_check_fastqc.fastq_1\":\"{}_1.fq.gz\",\"quality_check_fastqc_workflow.quality_check_fastqc.fastq_2\":\"{}_2.fq.gz\"}">{}-input.json'
```

run qc on an example file:

```
java -jar /opt/tools/cromwell-44.jar run https://gitlab.com/intelliseq/workflows/raw/master/src/main/wdl/tasks/quality-check-fastqc/v0.1/quality-check-fastqc.wdl -i Cac10-input.json > Cac10.txt
```

To run the workflow on all files in the folder (under screen):

```
ls *.json | xargs -i bash -c 'java -jar /opt/tools/cromwell-44.jar run https://gitlab.com/intelliseq/workflows/raw/master/src/main/wdl/tasks/quality-check-fastqc/v0.1/quality-check-fastqc.wdl -i {} > log-{}.txt'
```
To generate the MultiQC 1.7 report:

```
docker run --rm -v $PWD:/data ewels/multiqc:latest multiqc /data -o /data
```


