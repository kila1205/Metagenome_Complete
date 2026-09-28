process COVERM_ABUNDANCE {

    tag "${sample_id}"

    publishDir "${params.outdir}/06_abundance/${sample_id}",
        mode: 'copy'

    input:
    tuple val(sample_id),
          path(reads),
          path(bin_dir)

    output:
    tuple val(sample_id),
          path("${sample_id}.coverm.tsv"),
          emit: abundance

    script:
    """
    export TMPDIR=\$PWD/tmp
    mkdir -p \$TMPDIR

    coverm genome \
        --single ${reads} \
        --genome-fasta-directory ${bin_dir} \
        --genome-fasta-extension fa \
        --mapper minimap2-ont \
        --methods count relative_abundance mean \
        --min-read-percent-identity 80 \
        --min-covered-fraction 0 \
        --threads ${task.cpus} \
        --output-file ${sample_id}.coverm.tsv

    test -s ${sample_id}.coverm.tsv
    """
}
