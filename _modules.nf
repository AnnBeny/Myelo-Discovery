// ### UMI a read preprocessing ###
// # adaptors to trim: AACCGCCAGGAGT
// # UMI 1-8 bases in read position
process UMI_extract {
	tag "UMI_extract on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/temp/", mode:'copy'
	//label "s_cpu"
	//label "xxs_mem"
	container "quay.io/biocontainers/mulled-v2-452184f0fcbe6b0806405adb8f0b3873a7bd70a8:9c5efc89686d718e262836409bbf8541c6c5a159-0"
	
	input:
	tuple val(name), val(sample), path(fwd), path(rev)

	output:
	tuple val(name), val(sample), path("${name}.umi*.fastq.gz")

	script:
	""" 
	echo UMI_extract $name

	umi_tools extract -I ${fwd} --bc-pattern=NNNNNNNN --read2-in=${rev} --stdout=${name}.umi.R1.fastq.gz --read2-out=${name}.umi.R2.fastq.gz || \
    umi_tools extract --ignore-read-pair-suffixes -I ${fwd} --bc-pattern=NNNNNNNN --read2-in=${rev} --stdout=${name}.umi.R1.fastq.gz --read2-out=${name}.umi.R2.fastq.gz
	"""
	// --ignore-read-pair-suffixes // use this for MGI
}

//        source /home/kristina/miniconda3/bin/activate umi_tools

process TRIMMING {
	tag "TRIMMING on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/temp/", mode:'copy'

	//label "s_cpu"
	//label "xxs_mem"
	
	input:
	tuple val(name), val(sample), path(reads)

	output:
	tuple val(name), val(sample), path("*.fastq.gz")

	script:
	"""
	conda activate 'MYELO-ann'
	echo "=== Debug info ==="
    echo "Current PATH: \$PATH"
    echo "Cutadapt location: \$(which cutadapt || echo 'not found')"
	echo "Cutadapt version: \$(cutadapt --version || echo 'not found')"
    echo "Python version: \$(python --version)"
	echo "Conda env: \$CONDA_DEFAULT_ENV"
    
	echo TRIMMING $name
	
	cutadapt -g AACCGCCAGGAGT -m 50 \
		-o ${name}.trimmed1.R1.fastq.gz \
		-p ${name}.trimmed1.R2.fastq.gz $reads
	"""
}

//	source /home/kristina/miniconda3/bin/activate cutadapt

process FIRST_ALIGN_BAM {
	tag "FIRST_ALIGN_BAM on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/temp/", mode:'copy'
	conda 'MYELO-ann'
	//label "m_cpu"
	//label "l_mem"

	input:
	tuple val(name), val(sample), path(reads)

	output:
	tuple val(name), val(sample), path("${name}.bam")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	echo FIRST_ALIGN_BAM $name

	bwa mem -R ${rg} -t $task.cpus ${params.refindex} $reads > ${name}.sam
	samtools view -Sb ${name}.sam -o ${name}.bam
	"""
}

//	source /home/kristina/miniconda3/bin/activate bwa

process SORT_INDEX {
	tag "Sort index on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/temp/", mode:'copy'
	conda 'MYELO-ann'
	//label "m_mem"
	//label "xs_cpu"

	input:
	tuple val(name), val(sample), path(bam)

	output:
	tuple val(name), val(sample), path("${name}.sorted.bam"), path("${name}.sorted.bai")


	script:
	"""
	echo SORT_INDEX $name
	samtools sort $bam -o ${name}.sorted.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bai
	"""
}


process DEDUP {
	tag "DEDUP on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
	conda 'MYELO-ann'
	//label "s_cpu"
	//label "m_mem"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.deduped.bam"), path("${name}.deduped.bai")

	script:
	"""
	echo DEDUP $name

	umi_tools dedup -I $bam --paired -S ${name}.deduped.bam
	samtools index ${name}.deduped.bam ${name}.deduped.bai
	"""
}

//	source /home/kristina/miniconda3/bin/activate umi_tools

process MUTECT2 {
	tag "MUTECT2 on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/vcfs/", mode:'copy'
	conda 'gatk4'

	//label "m_mem"
	//label "s_cpu"

	input:
	tuple val(name), val(sample), path(bam), path(bai)
	
	output:
	tuple val(name), val(sample), path ("${sample.name}.mutect.vcf"), path ("${sample.name}.mutect.vcf.stats")

	script:
	"""
    # Set JAVA_HOME explicitly to use Java 8
	unset JAVA_HOME

    export JAVA_HOME=/anna/miniconda3/envs/gatk4/jre
    export PATH=\$JAVA_HOME/bin:\$PATH
    
    echo "=== Debug info ==="
    echo "Java version: \$(java -version 2>&1)"
    echo "JAVA_HOME: \$JAVA_HOME"

	echo MUTECT2 $name
	/home/anna/SOFT/gatk4/gatk-4.1.0.0/gatk Mutect2 \
		--reference ${params.ref}.fa \
		--input ${bam} \
		--tumor-sample $name \
		--annotation StrandArtifact \
		--min-base-quality-score 20 \
		--intervals $params.ivl \
		--output ${sample.name}.mutect.vcf

	STATS_FILE=\$(find ${task.workDir} -name "${sample.name}.mutect.vcf.stats")
	"""
}

/* 
	echo MUTECT2 $name
	echo "Reference: ${params.ref}.fa"
    echo "BAM: ${bam}"
    echo "BAI: ${bai}"
    echo "Intervals: ${params.bed}"


	/home/anna/SOFT/gatk/gatk Mutect2 \
		--R ${params.ref}.fa \
		-I ${bam} \
		-tumor-sample ${sample.name} \
		-O ${sample.name}.mutect.vcf \
		-L $params.bed \
		--min-base-quality-score 20 \
		--genotype-pon-sites false \
		--genotype-germline-sites false \
		--af-of-alleles-not-in-resource -1.0 \
		--tumor-lod-to-emit 3.0 \
		--initial-tumor-lod 2.0 \
		--pcr-indel-model CONSERVATIVE \
		--max-reads-per-alignment-start 50 \
		--native-pair-hmm-threads $task.cpus

	gatk FilterMutectCalls \
		-V ${sample.name}.mutect.vcf \
		-R ${params.ref}.fa \
		-O ${sample.name}.filtered.vcf

	echo "Checking for stats file..."
    ls -lh ${task.workDir}
*/

/*
	/home/anna/SOFT/gatk/gatk Mutect2 \
		--reference ${params.ref}.fa \
		--input ${bam} \
		--tumor-sample $name \
		--min-base-quality-score 20 \
		--intervals $params.bed \
		--output ${sample.name}.mutect.vcf 
*/
//	/home/anna/SOFT/gatk/gatk Mutect2 --reference ${params.ref}.fa --input ${bam} --tumor-sample $name --min-base-quality-score 20 --intervals $params.ivl --output ${sample.name}.mutect.vcf
//	source /home/kristina/miniconda3/bin/activate gatk4

process FILTER_MUTECT {
	tag "FILTER_MUTECT on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/temp/", mode:'copy'
	conda 'MYELO-ann'
	//label "s_mem"
	//label "s_cpu"

	input:
	tuple val(name), val(sample), path(vcf_input), path(stats)
	
	output:
	tuple val(name), val(sample), path ("${sample.name}.mutect.filt.vcf")

	script:
	"""
	echo FILTER_MUTECT $name

	/home/anna/SOFT/gatk4/gatk-4.1.0.0/gatk FilterMutectCalls \
		-V $vcf_input  \
		-O ${sample.name}.mutect.filt.vcf \
		--reference ${params.ref}.fa \
		--stats ${stats}
	"""
}

//	source /home/kristina/miniconda3/bin/activate gatk4

process NORMALIZE_MUTECT {
	tag "NORMALIZE_MUTECT on $name using $task.cpus CPUs $task.memory"
	conda 'MYELO-ann'
	//label "xxs_mem"
	//label "s_cpu"
	//container "staphb/bcftools:1.10.2"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(name), val(sample), path ("${sample.name}.mutect.filt.norm.vcf")

	script:
	"""
	echo NORMALIZE_MUTECT $name
	source /home/kristina/miniconda3/bin/activate bcftools
	bcftools norm -m-both $vcf_input > ${sample.name}.mutect.filt.norm.vcf
	"""
}

//        source /home/kristina/miniconda3/bin/activate bcftools

process ANNOTATE_MUTECT {
	tag "ANNOTATE_MUTECT on $name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/vcfs/", mode:'copy'
	conda 'MYELO-ann'

	//label "s_mem"
	//label "s_cpu"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(name), val(sample), path("${sample.name}.mutect2.filt.norm.vep.vcf")

	script:
	"""
	echo ANNOTATE_MUTECT $name
	source /home/kristina/miniconda3/bin/activate  vep95
	vep -i $vcf_input \
		--cache --cache_version 95 --dir_cache $params.vep \
		--fasta ${params.vep_fasta} --offline \
		--vcf --everything -o ${sample.name}.mutect2.filt.norm.vep.vcf
	"""	
}

//vep -i $vcf_input --cache --cache_version 113 --dir_cache $params.vep \
//	source /home/kristina/miniconda3/bin/activate  vep

process FILTER_VCF {
	tag "FILTER_VCF on $name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/vcfs/", mode:'copy'
	conda 'MYELO-ann'
	//container "staphb/bcftools:1.10.2"
	//label "xxs_mem"
	//label "s_cpu"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(name), val(sample), path("${sample.name}.mutect2.filt.norm.vep.filt.vcf")

	script:
	"""
	echo FILTER_VCF $name

	bcftools view -f 'PASS,clustered_events,multiallelic' $vcf_input > ${sample.name}.mutect2.filt.norm.vep.filt.vcf
	"""	
}

//        source /home/kristina/miniconda3/bin/activate bcftools

process VCF2CSV {
	tag "VCF2CSV on $name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/variants/", mode:'copy'
	conda 'MYELO-ann'
	//label "xs_mem"
	//label "s_cpu"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(sample.run), path("${name}.csv")

	script:
	"""
	echo VCF2CSV $name

	python $params.vcf2csv simple --build GRCh37 -i $vcf_input -o ${name}.csv
	"""	
}

//	source /home/kristina/miniconda3/bin/activate vcf2csv

process MERGE_TABLES {
	tag "MERGE_TABLES on $run using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${run}/variants/", mode:'copy'
	conda 'MYELO-ann'
	//label "s_mem"
	//label "s_cpu"
	debug true

	input:
	tuple val(run), path(all_annotated_normed)
	
	output:
	path "${run}.merged_variant.table_NEW.tsv"

	script:
	"""
	echo MERGE_TABLES $run

	Rscript --vanilla $params.mergescript $run
	"""	
}

//	source /home/kristina/miniconda3/bin/activate erko

process FLT3 {
	tag "FLT3 on $name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/FLT3/", mode:'copy'
	conda 'MYELO-ann'
	//label "m_mem"
	//label "s_cpu"
	errorStrategy { task.exitStatus in [143,137,104,134,139,247,null,'', '-'] ? 'retry' : 'ignore' }

	input:
	tuple val(name), val(sample), path(bam), path(bai)
	
	output:
	tuple val(name), val(sample), path("${name}.deduped_FLT3_ITD_summary.txt")

	script:
	"""

	echo FLT3 $name
    tar -C /home/anna/SOFT -xf $params.flt3tar
	perl /home/anna/SOFT/FLT3/FLT3_ITD_ext/FLT3_ITD_ext.pl --bam $bam --output ./ --ngstype amplicon --genome hg19 --fgbiojar $params.fgbio --picardjar $params.picard --refindex /mnt/hdd2/anna/Myelo/src/project/xsvato01/archer_nf/bin/FLT3/FLT3_bwaindex/FLT3_dna_e1415
	touch ${name}.deduped_FLT3_ITD_summary.txt
	ls -alh
	"""	
}

// source /home/anna/miniconda3/bin/activate perl
// perl /tmp/FLT3/FLT3_ITD_ext/FLT3_ITD_ext.pl
//	source /home/kristina/miniconda3/bin/activate perl


process BAMQC {
	tag "BAMQC on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/temp/", mode:'copy'
	conda 'MYELO-ann'
	//label "s_mem"
	//label "s_cpu"
	//container 'quay.io/biocontainers/mulled-v2-b0664646864bfdb46c5343b1b2b93fc05adb4b77:39a005770a3e30fb6aa3bf424b57ddf52bae7ece-0'

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("*")

	script:
	"""
	echo BAMQC $name

	samtools flagstat $bam > ${name}.flagstat
	samtools stats $bam > ${name}.samstats
	picard CollectHsMetrics I=$bam BAIT_INTERVALS=$params.ivl TARGET_INTERVALS=$params.ivl R=${params.ref}.fa O=${name}.hs_metrics
	picard CollectAlignmentSummaryMetrics I=$bam R=${params.ref}.fa O=${name}.aln_metrics
	"""
}

//        source /home/kristina/miniconda3/bin/activate picard

process COVERAGE {
	tag "COVERAGE on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/temp/", mode:'copy'
	conda 'MYELO-ann'
	//label "xl_mem"
	//label "s_cpu"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.PBcov.cons.txt")

	script:
	"""
	echo COVERAGE $name

	bedtools coverage -abam $params.bed -b $bam -d > ${name}.PBcov.cons.txt
	"""
}

//	source /home/kristina/miniconda3/bin/activate bedtools

process COVERAGE_POSTPROCESS {
	tag "COVERAGE_POSTPROCESS on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/Cov/", mode:'copy'
	conda 'MYELO-ann'
	//label "s_mem"
	//label "xxs_mem"

	input:
	tuple val(name), val(sample), path(txt)

	output:
	tuple val(name), val(sample), path("${name}.perexon_stat.txt")

	script:
	"""
	echo COVERAGE_POSTPROCESS $name

	Rscript --vanilla $params.covscript $txt
	ls -al
	"""
}

//	source /home/kristina/miniconda3/bin/activate erko

process MULTIQC {
	tag "MultiQC on all samples from $run using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${run}/QC/", mode:'copy'
	conda 'MYELO-ann'
	//container 'ewels/multiqc:v1.18'
	//label "xs_mem"
	//label "s_cpu"

	input:
	tuple val(run), path(all_annotated_normed)

	output:
	path "*"

	script:
	"""
	echo MULTIQC $run 
	PYTHONPATH=/home/anna/SOFT/multiqc/lib/python3.9/site-packages /home/anna/SOFT/multiqc/bin/multiqc . -n MultiQC.html
	"""

}

