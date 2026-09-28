process KRAKEN2_CLASSIFY {

    tag "${sample_id}"

    publishDir "${params.outdir}/04_read_taxonomy/kraken2/${sample_id}",
        mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val kraken_db

    output:
    tuple val(sample_id),
          path("${sample_id}.kraken2.report"),
          emit: report

    tuple val(sample_id),
          path("${sample_id}.kraken2.output"),
          emit: classifications

    script:
    """
    kraken2 \
        --db ${kraken_db} \
        --threads ${task.cpus} \
        --report ${sample_id}.kraken2.report \
        --output ${sample_id}.kraken2.output \
        ${reads}

    test -s ${sample_id}.kraken2.report
    """
}


process BRACKEN_ABUNDANCE {

    tag "${sample_id}"

    publishDir "${params.outdir}/04_read_taxonomy/bracken/${sample_id}",
        mode: 'copy'

    input:
    tuple val(sample_id), path(kraken_report)
    val kraken_db
    val bracken_read_length
    val bracken_level
    val bracken_threshold

    output:
    tuple val(sample_id),
          path("${sample_id}.bracken.tsv"),
          emit: abundance

    tuple val(sample_id),
          path("${sample_id}.bracken.report"),
          emit: report

    script:
    """
    bracken \
        -d ${kraken_db} \
        -i ${kraken_report} \
        -o ${sample_id}.bracken.tsv \
        -w ${sample_id}.bracken.report \
        -r ${bracken_read_length} \
        -l ${bracken_level} \
        -t ${bracken_threshold}

    test -s ${sample_id}.bracken.tsv
    """
}
