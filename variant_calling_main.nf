#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Variant Calling Pipeline using BCFtools
 * Per-Sample Variant Calling - Each BAM gets its own VCF
 * Input: BAM files from multiple library directories
 * Output: Individual filtered VCF files per sample
 * 
 * Sample naming with library suffix to avoid collisions:
 * - library1: Peromyscus123.bam -> Peromyscus123_1.vcf.gz
 * - library2: Peromyscus123.bam -> Peromyscus123_2.vcf.gz
 * - library3: Peromyscus123.bam -> Peromyscus123_3.vcf.gz
 * 
 * Duplicates across libraries are kept separate for QC comparison
 */

// Parameters
params.bam_dirs = [
    '/home/sydt/scratch/alignment/library1/aligned',
    '/home/sydt/scratch/alignment/library2/aligned',
    '/home/sydt/scratch/alignment/library3/aligned'
]
params.outdir = '/home/sydt/scratch/variant_calling/results'
params.reference = '/home/sydt/scratch/alignment/genome_index/GCF_049852395.1_HU_Pman_BW_mat_3.1_genomic.fna'

// BCFtools parameters (matching your example)
params.min_qual = 19
params.min_gq = 9
params.platform = 'ILLUMINA'

log.info """
=========================================
Variant Calling Pipeline - BCFtools
Per-Sample Mode (Individual VCFs)
=========================================
BAM directories    : ${params.bam_dirs.join(', ')}
Reference genome   : ${params.reference}
Output directory   : ${params.outdir}
Min QUAL score     : ${params.min_qual}
Min GQ score       : ${params.min_gq}
Platform           : ${params.platform}
Output             : One VCF per sample (~1585 VCFs)
Library suffix     : _1, _2, _3 added to filenames
=========================================
"""

/*
 * Process: Index reference genome for BCFtools
 */
process INDEX_REFERENCE {
    tag "reference"
    
    input:
    path reference
    
    output:
    tuple path(reference), path("${reference}.fai")
    
    script:
    """
    samtools faidx ${reference}
    """
}

/*
 * Process: BCFtools mpileup - Per-sample variant calling
 * Runs individually on each BAM file
 */
process BCFTOOLS_MPILEUP {
    tag "${sample_id}"
    
    cpus 4
    memory '16 GB'
    time '12h'
    
    input:
    tuple val(sample_id), path(bam), path(bai), path(reference), path(reference_fai)
    
    output:
    tuple val(sample_id), path("${sample_id}.bcf")
    
    script:
    """
    # Run bcftools mpileup on single BAM file
    bcftools mpileup \\
        -a DP,AD \\
        --skip-indels \\
        -P ${params.platform} \\
        -f ${reference} \\
        ${bam} \\
        -o ${sample_id}.bcf
    
    echo "BCF file generated for ${sample_id}"
    """
}

/*
 * Process: BCFtools call and filter - Per-sample VCF generation
 * Generates individual filtered VCF for each sample
 */
process BCFTOOLS_CALL_AND_FILTER {
    tag "${sample_id}"
    publishDir "${params.outdir}/vcf_files", mode: 'copy', pattern: "${sample_id}.vcf.gz"
    publishDir "${params.outdir}/vcf_files", mode: 'copy', pattern: "${sample_id}.vcf.gz.csi"
    publishDir "${params.outdir}/stats", mode: 'copy', pattern: "${sample_id}_stats.txt"
    
    cpus 2
    memory '8 GB'
    time '4h'
    
    input:
    tuple val(sample_id), path(bcf)
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.csi")
    path "${sample_id}_stats.txt"
    
    script:
    """
    # Call variants, filter, and output VCF for single sample
    bcftools call \\
        -m \\
        --variants-only \\
        --format-fields GQ \\
        --skip-variants indels \\
        ${bcf} \\
    | bcftools filter \\
        --set-GTs . \\
        --include 'QUAL > ${params.min_qual} && FMT/GQ > ${params.min_gq}' \\
    | bcftools view \\
        --min-alleles 2 \\
        --max-alleles 2 \\
        --types snps \\
        --apply-filter "PASS" \\
        --output-type v \\
        --output-file ${sample_id}.vcf
    
    # Compress and index
    bgzip -c ${sample_id}.vcf > ${sample_id}.vcf.gz
    bcftools index ${sample_id}.vcf.gz
    
    # Generate variant statistics for this sample
    echo "Variant Calling Statistics - ${sample_id}" > ${sample_id}_stats.txt
    echo "========================================" >> ${sample_id}_stats.txt
    echo "" >> ${sample_id}_stats.txt
    echo "Total variants:" >> ${sample_id}_stats.txt
    bcftools view -H ${sample_id}.vcf | wc -l >> ${sample_id}_stats.txt
    echo "" >> ${sample_id}_stats.txt
    echo "Variants per chromosome:" >> ${sample_id}_stats.txt
    bcftools view -H ${sample_id}.vcf | cut -f1 | sort | uniq -c >> ${sample_id}_stats.txt
    
    # Clean up uncompressed VCF
    rm ${sample_id}.vcf
    
    echo "Completed variant calling for ${sample_id}"
    """
}

/*
 * Main workflow
 */
workflow {
    // Index reference genome once
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
    indexed_reference = INDEX_REFERENCE(reference_ch)
    
    // Create channels for each library with library suffix in sample_id
    // Library 1 - add _1 suffix
    library1_ch = Channel
        .fromFilePairs("${params.bam_dirs[0]}/*.{bam,bam.bai}", size: 2, checkIfExists: true)
        .map { sample_id, files ->
            def bam = files.find { it.name.endsWith('.bam') && !it.name.endsWith('.bam.bai') }
            def bai = files.find { it.name.endsWith('.bam.bai') }
            def new_id = "${sample_id}_1"
            tuple(new_id, bam, bai)
        }
    
    // Library 2 - add _2 suffix
    library2_ch = Channel
        .fromFilePairs("${params.bam_dirs[1]}/*.{bam,bam.bai}", size: 2, checkIfExists: true)
        .map { sample_id, files ->
            def bam = files.find { it.name.endsWith('.bam') && !it.name.endsWith('.bam.bai') }
            def bai = files.find { it.name.endsWith('.bam.bai') }
            def new_id = "${sample_id}_2"
            tuple(new_id, bam, bai)
        }
    
    // Library 3 - add _3 suffix
    library3_ch = Channel
        .fromFilePairs("${params.bam_dirs[2]}/*.{bam,bam.bai}", size: 2, checkIfExists: true)
        .map { sample_id, files ->
            def bam = files.find { it.name.endsWith('.bam') && !it.name.endsWith('.bam.bai') }
            def bai = files.find { it.name.endsWith('.bam.bai') }
            def new_id = "${sample_id}_3"
            tuple(new_id, bam, bai)
        }
    
    // Combine all libraries (each sample keeps its library suffix)
    all_samples = library1_ch.mix(library2_ch, library3_ch)
    
    // Add reference to each sample and flatten the tuple
    samples_with_ref = all_samples
        .combine(indexed_reference)
        .map { sample_id, bam, bai, ref, ref_fai ->
            tuple(sample_id, bam, bai, ref, ref_fai)
        }
    
    // Run mpileup on each sample individually
    bcf_files = BCFTOOLS_MPILEUP(samples_with_ref)
    
    // Call variants and filter for each sample
    BCFTOOLS_CALL_AND_FILTER(bcf_files)
}

workflow.onComplete {
    log.info """
    =========================================
    Pipeline completed!
    =========================================
    Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration  : ${workflow.duration}
    Output dir: ${params.outdir}
    =========================================
    That's a wrap! Great work on set today guys!
    """
}
