# ifpan-starowicz-knees
RNAseq analysis

STEP 1: Quality Control with fastqc

generate json input files based on the file_1.fq.gz file_2.fq.gz naming convention:

```ls *1.fq.gz | xargs -i bash -c 'BASENAME2=$(echo {} | cut -d "." -f 1 | cut -d "_" -f 1); echo $BASENAME2' | xargs -i bash -c 'echo "{\"generate_fastqc_report_workflow.generate_fastqc_report.fastq_1\":\"{}_1.fq.gz\",\"generate_fastqc_report_workflow.generate_fastqc_report.fastq_2\":\"{}_2.fq.gz\"}">{}-input.json'
```

run qc on an example file:

```
java -jar /opt/tools/cromwell-44.jar run https://gitlab.com/intelliseq/workflows/raw/master/src/main/wdl/tasks/quality-check-fastqc/v0.1/quality-check-fastqc.wdl -i Cac10-input.json > Cac10.txt
```

To run the workflow on all files in the folder: ls *.json | xargs -i bash -c 'java -jar /opt/tools/cromwell-44.jar run https://gitlab.com/intelliseq/workflows/raw/master/src/main/wdl/tasks/quality-check-fastqc/v0.1/quality-check-fastqc.wdl -i {}-input.json > log-{}.txt' & *the '&' sign at the end of line tells bash to run whatever command in the background

