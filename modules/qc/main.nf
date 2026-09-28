process RAW_NANOPLOT {

    tag "${sample_id}"

    publishDir "${params.outdir}/01_qc/01_raw_nanoplot/${sample_id}",
        mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path(reads),
          path("${sample_id}_raw_nanoplot"),
          emit: qc

    script:
    """
    mkdir -p ${sample_id}_raw_nanoplot

    NanoPlot \
        --fastq ${reads} \
        --outdir ${sample_id}_raw_nanoplot \
        --threads ${task.cpus}

    test -s ${sample_id}_raw_nanoplot/NanoStats.txt
    """
}


process FILTLONG {

    tag "${sample_id}"

    publishDir "${params.outdir}/01_qc/02_filtered_fastq",
        mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(raw_qc)

    output:
    tuple val(sample_id),
          path("${sample_id}.filtered.fastq.gz"),
          emit: filtered

    path "${sample_id}.filtlong.log"

    script:
    """
    filtlong \
        --min_length 1000 \
        --keep_percent 90 \
        ${reads} \
        2> ${sample_id}.filtlong.log \
    | pigz -p ${task.cpus} \
    > ${sample_id}.filtered.fastq.gz

    pigz -t ${sample_id}.filtered.fastq.gz
    test -s ${sample_id}.filtered.fastq.gz
    """
}


process FILTERED_NANOPLOT {

    tag "${sample_id}"

    publishDir "${params.outdir}/01_qc/03_filtered_nanoplot/${sample_id}",
        mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path(reads),
          path("${sample_id}_filtered_nanoplot"),
          emit: qc

    script:
    """
    mkdir -p ${sample_id}_filtered_nanoplot

    NanoPlot \
        --fastq ${reads} \
        --outdir ${sample_id}_filtered_nanoplot \
        --threads ${task.cpus}

    test -s ${sample_id}_filtered_nanoplot/NanoStats.txt
    """
}


process MULTIQC {

    tag "multiqc"

    publishDir "${params.outdir}/01_qc/04_multiqc",
        mode: 'copy'

    input:
    path qc_directories

    output:
    path "multiqc_report.html",
         emit: report

    script:
    """
    multiqc \
        ${qc_directories} \
        --filename multiqc_report.html \
        --force

    test -s multiqc_report.html
    """
}
