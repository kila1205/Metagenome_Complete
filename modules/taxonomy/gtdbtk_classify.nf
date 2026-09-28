process GTDBTK_CLASSIFY {

    tag "${sample_id}"

    conda "${projectDir}/envs/taxonomy.yml"

    input:
    tuple val(sample_id), path(bin_dir)
    val gtdbtk_data_path

    output:
    tuple val(sample_id),
          path("gtdbtk_out"),
          emit: gtdbtk_results

    script:
    """
    mkdir -p gtdbtk_out

    export GTDBTK_DATA_PATH="${gtdbtk_data_path}"

    gtdbtk classify_wf \
        --genome_dir ${bin_dir} \
        --out_dir gtdbtk_out \
        --extension fa \
        --cpus ${task.cpus}

    test -d gtdbtk_out
    """
}
