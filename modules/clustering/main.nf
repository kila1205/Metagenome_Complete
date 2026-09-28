nextflow.enable.dsl = 2


process COMBINE_PROTEINS {

    tag "${catalog_id}"

    publishDir "${params.outdir}/06_clustering/01_combined_proteins/${catalog_id}",
        mode: 'copy'

    conda "${projectDir}/envs/clustering.yml"

    cpus 2
    memory '8 GB'

    input:
    tuple val(catalog_id), path(protein_files)

    output:
    tuple val(catalog_id),
          path("${catalog_id}.all_proteins.faa"),
          emit: combined

    script:
    """
    cat ${protein_files} > ${catalog_id}.all_proteins.faa

    test -s ${catalog_id}.all_proteins.faa
    """
}


process MMSEQS_CLUSTER {

    tag "${catalog_id}"

    publishDir "${params.outdir}/06_clustering/02_mmseqs/${catalog_id}",
        mode: 'copy'

    conda "${projectDir}/envs/clustering.yml"

    cpus 16
    memory '48 GB'

    input:
    tuple val(catalog_id), path(proteins)
    val min_seq_id
    val coverage

    output:
    tuple val(catalog_id),
          path("${catalog_id}.representative_proteins.faa"),
          emit: representatives

    tuple val(catalog_id),
          path("${catalog_id}.cluster_membership.tsv"),
          emit: membership

    script:
    """
    mkdir -p tmp

    mmseqs createdb \
        ${proteins} \
        proteins_db

    mmseqs cluster \
        proteins_db \
        cluster_db \
        tmp \
        --min-seq-id ${min_seq_id} \
        -c ${coverage} \
        --cov-mode 0 \
        --threads ${task.cpus}

    mmseqs createtsv \
        proteins_db \
        proteins_db \
        cluster_db \
        ${catalog_id}.cluster_membership.tsv

    mmseqs createsubdb \
        cluster_db \
        proteins_db \
        representative_db

    mmseqs convert2fasta \
        representative_db \
        ${catalog_id}.representative_proteins.faa

    test -s ${catalog_id}.representative_proteins.faa
    test -s ${catalog_id}.cluster_membership.tsv
    """
}


process CATALOG_STATS {

    tag "${catalog_id}"

    publishDir "${params.outdir}/06_clustering/03_catalog_stats/${catalog_id}",
        mode: 'copy'

    conda "${projectDir}/envs/clustering.yml"

    cpus 2
    memory '8 GB'

    input:
    tuple val(catalog_id), path(representatives)

    output:
    tuple val(catalog_id),
          path("${catalog_id}.protein_catalog_stats.txt"),
          emit: stats

    script:
    """
    seqkit stats -a -T \
        ${representatives} \
        > ${catalog_id}.protein_catalog_stats.txt

    test -s ${catalog_id}.protein_catalog_stats.txt
    """
}
