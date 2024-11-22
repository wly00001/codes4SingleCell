#!/bin/bash


cellranger count --id=PBS \
   --fastqs=GSE160450/rawfastq \
   --sample=SRR12937940 \
   --transcriptome=reference/CellRanger/refdata-gex-mm10-2020-A



cellranger count --id=LPS \
   --fastqs=GSE160450/rawfastq \
   --sample=SRR12937941 \
   --transcriptome=reference/CellRanger/refdata-gex-mm10-2020-A


