process MERGE_READS {

    tag "${sample_id}"

    publishDir "${params.outdir}/00_merged_reads",
        mode: 'copy'

    input:
    tuple val(sample_id), path(read_files)

    output:
    tuple val(sample_id),
          path("${sample_id}.raw.fastq.gz"),
          emit: merged_reads

    script:
    """
    echo "Sample: ${sample_id}"
    echo "Input FASTQ files:"
    printf '%s\n' ${read_files}

    cat ${read_files} | pigz -p ${task.cpus} > ${sample_id}.raw.fastq.gz

    pigz -t ${sample_id}.raw.fastq.gz

    test -s ${sample_id}.raw.fastq.gz
    """
}
