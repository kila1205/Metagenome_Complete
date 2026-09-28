process PARSE_TAXONOMY {

    tag "${sample_id}"

    publishDir "${params.outdir}/08_taxonomy/gtdbtk/${sample_id}",
        mode: 'copy'

    conda "${projectDir}/envs/taxonomy.yml"

    cpus 2
    memory '8 GB'
    time '2h'

    input:
    tuple val(sample_id), path(gtdbtk_summary)
    path parser_script

    output:
    tuple val(sample_id),
          path("${sample_id}.taxonomy.tsv"),
          emit: taxonomy_table

    script:
    """
    python ${parser_script} \
        ${gtdbtk_summary} \
        ${sample_id}.taxonomy.tsv

    test -s ${sample_id}.taxonomy.tsv
    """
}
