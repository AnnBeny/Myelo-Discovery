# nextflow kuberun xsvato01/TP53_nf -r main -head-image 'cerit.io/nextflow/nextflow:22.11.1' -resume -with-report -c nextflow.config -params-file samplesheet.json 

nextflow run /home/ciri/workflow/myelo/src/project/xsvato01/archer_nf -c nextflow.config -profile standard -resume -with-report