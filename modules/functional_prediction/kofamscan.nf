process KOFAMSCAN {

    tag "${sample_id}"

    publishDir "${params.outdir}/09_functional_prediction/kofamscan/${sample_id}",
        mode: 'copy'

    conda "${projectDir}/envs/kofamscan.yml"

    cpus 16
    memory '64 GB'
    time '48h'

    input:
    tuple val(sample_id), path(proteins)
    val kofam_profiles
    val kofam_ko_list

    output:
    tuple val(sample_id),
          path("${sample_id}.kofam.tsv"),
          emit: annotations

    script:
    """
    exec_annotation \
        --profile ${kofam_profiles} \
        --ko-list ${kofam_ko_list} \
        --cpu ${task.cpus} \
        --format detail-tsv \
        -o ${sample_id}.kofam.tsv \
        ${proteins}

    test -s ${sample_id}.kofam.tsv
    """
}
