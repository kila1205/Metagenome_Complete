process EGGNOG_MAPPER {

    tag "${sample_id}"

    publishDir "${params.outdir}/10_functional_annotation/eggnog/${sample_id}",
        mode: 'copy'

    conda "${projectDir}/envs/eggnog.yml"

    cpus 16
    memory '64 GB'
    time '48h'

    input:
    tuple val(sample_id), path(proteins)
    val eggnog_data_dir

    output:
    tuple val(sample_id),
          path("${sample_id}.emapper.annotations"),
          emit: annotations

    tuple val(sample_id),
          path("${sample_id}.emapper.seed_orthologs"),
          emit: seed_orthologs

    // eggNOG-mapper does not always create a .log file.
    // Keep it optional because downstream analysis only requires
    // the annotation and seed-ortholog outputs.
    path "${sample_id}.emapper.log",
         emit: log,
         optional: true

    script:
    """
    emapper.py \
        -i ${proteins} \
        --itype proteins \
        -m diamond \
        --cpu ${task.cpus} \
        --data_dir ${eggnog_data_dir} \
        --output ${sample_id}

    test -s ${sample_id}.emapper.annotations
    """
}
