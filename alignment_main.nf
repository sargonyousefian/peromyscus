#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Create a channel for input reads - handling single-end FASTQ files
reads_ch = Channel
    .fromPath(params.input)
    .map { file ->
        // Extract sample ID from filename (e.g., Peromyscus99.fq -> Peromyscus99)
        def sample_id = file.baseName
        return tuple(sample_id, file)
    }

// Debug: View the contents of input channels
reads_ch.view{ "Sample: $it" }

// Trim adapters using fastp
process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*.fq.gz"
    publishDir "${params.outdir}/fastp_reports", mode: 'copy', pattern: "*.{html,json}"

    module 'StdEnv/2023:fastp/0.24.0'
    
    cpus 2
    memory '4 GB'
    time '4h'

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.html", emit: html
    path "${sample_id}_fastp.json", emit: json

    script:
    """
    fastp \
        -i ${read} \
        -o ${sample_id}_trimmed.fq.gz \
        --thread ${task.cpus} \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --adapter_sequence AGATCGGAAGAGC \
        --trim_front1 6 \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json
    
    echo "Trimming completed for ${sample_id}"
    """
}

// Run Quality Control on FastQ Files (original)
process FASTQC_RAW {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    
    module 'StdEnv/2023:fastqc/0.12.1'
    
    cpus 1
    memory '4 GB'
    time '2h'

    input:
    tuple val(sample_id), path(read)

    output:
    path "fastqc_${sample_id}_raw"

    script:
    """
    mkdir fastqc_${sample_id}_raw
    fastqc -o fastqc_${sample_id}_raw -f fastq -q ${read}
    """
}

// Run Quality Control on trimmed reads
process FASTQC_TRIMMED {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy'
    
    module 'StdEnv/2023:fastqc/0.12.1'
    
    cpus 1
    memory '4 GB'
    time '2h'

    input:
    tuple val(sample_id), path(read)

    output:
    path "fastqc_${sample_id}_trimmed"

    script:
    """
    mkdir fastqc_${sample_id}_trimmed
    fastqc -o fastqc_${sample_id}_trimmed -f fastq -q ${read}
    """
}

// Align TRIMMED reads using BWA and sort with Samtools
process ALIGN_AND_SORT {
    tag "$sample_id"
    publishDir "${params.outdir}/aligned", mode: 'copy', pattern: "*.{bam,bai}"
    publishDir "${params.outdir}/stats", mode: 'copy', pattern: "*.{flagstat,stats}"

    module 'StdEnv/2023:bwa/0.7.17:samtools/1.22.1'
    
    cpus 4
    memory '16 GB'
    time '20h'

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam
    path "${sample_id}.flagstat", emit: flagstat
    path "${sample_id}.stats", emit: stats

    script:
    def read_group = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA"
    
    """
    echo "Processing TRIMMED single-end sample: ${sample_id}"
    echo "Using BWA index: ${params.bwa_index_prefix}"
    
    # Decompress if needed and align
    gunzip -c ${read} | \
    bwa mem -t ${task.cpus} -R "${read_group}" ${params.bwa_index_prefix} - | \
    samtools sort -@ ${task.cpus} -o ${sample_id}.bam -
    
    samtools index ${sample_id}.bam
    
    # Quality checks
    samtools flagstat ${sample_id}.bam > ${sample_id}.flagstat
    samtools stats ${sample_id}.bam > ${sample_id}.stats
    
    echo "Alignment completed for ${sample_id}"
    
    # Print alignment summary
    echo "Summary for ${sample_id}:"
    samtools flagstat ${sample_id}.bam | grep "mapped ("
    """
}

// Generate depth per position
process SAMTOOLS_DEPTH {
    tag "$sample_id"
    publishDir "${params.outdir}/depth", mode: 'copy'

    module 'StdEnv/2023:samtools/1.22.1'
    
    cpus 2
    memory '4 GB'
    time '2h'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}_depth.txt"

    script:
    """
    samtools depth ${bam} > ${sample_id}_depth.txt

    if [ ! -s ${sample_id}_depth.txt ]; then
        echo "Error: ${sample_id}_depth.txt is empty"
        exit 1
    fi
    
    echo "Depth file created for ${sample_id}"
    wc -l ${sample_id}_depth.txt
    """
}

// Define the main workflow
workflow {
    log.info ""
    log.info "========================================"
    log.info "  Peromyscus Alignment Pipeline"
    log.info "  WITH FASTP ADAPTER TRIMMING"
    log.info "========================================"
    log.info "Input reads       : ${params.input}"
    log.info "BWA index prefix  : ${params.bwa_index_prefix}"
    log.info "Output directory  : ${params.outdir}"
    log.info "========================================"
    log.info ""

    // Run FastQC on raw reads
    FASTQC_RAW(reads_ch)
    
    // Trim adapters with fastp
    FASTP(reads_ch)
    
    // Run FastQC on trimmed reads
    FASTQC_TRIMMED(FASTP.out.trimmed_reads)
    
    // Align trimmed reads
    aligned = ALIGN_AND_SORT(FASTP.out.trimmed_reads)

    // Generate depth
    SAMTOOLS_DEPTH(aligned.bam)
}
