process PRODIGAL {

    tag "${sample_id}"

    publishDir "${params.outdir}/03_gene_prediction/01_predictions/${sample_id}",
        mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id),
          path("${sample_id}.genes.gff"),
          path("${sample_id}.genes.fna"),
          path("${sample_id}.proteins.faa"),
          emit: predictions

    script:
    """
    prodigal \
        -i ${assembly} \
        -p meta \
        -f gff \
        -o ${sample_id}.genes.gff \
        -d ${sample_id}.genes.fna \
        -a ${sample_id}.proteins.faa

    test -s ${sample_id}.genes.gff
    test -s ${sample_id}.genes.fna
    test -s ${sample_id}.proteins.faa
    """
}


process PROTEIN_STATS {

    tag "${sample_id}"

    publishDir "${params.outdir}/03_gene_prediction/02_protein_stats",
        mode: 'copy'

    input:
    tuple val(sample_id),
          path(gff),
          path(genes),
          path(proteins)

    output:
    path "${sample_id}.protein_stats.txt"

    script:
    """
    seqkit stats -T \
        -a \
        ${proteins} \
        > ${sample_id}.protein_stats.txt

    test -s ${sample_id}.protein_stats.txt
    """
}
