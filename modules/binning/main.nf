nextflow.enable.dsl=2


process MAP_READS {

    tag "${sample_id}"

    publishDir "${params.outdir}/04_binning/01_mapping/${sample_id}",
        mode: 'copy'

    conda "${projectDir}/envs/binning.yml"

    cpus 16
    memory '48 GB'

    input:
    tuple val(sample_id), path(assembly), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai"),
          emit: mapped

    script:
    """
    minimap2 \
        --split-prefix ${sample_id}_mmi \
        -ax map-ont \
        -t ${task.cpus} \
        ${assembly} \
        ${reads} \
    | samtools sort \
        -@ ${task.cpus} \
        -o ${sample_id}.sorted.bam

    samtools index \
        -@ ${task.cpus} \
        ${sample_id}.sorted.bam

    test -s ${sample_id}.sorted.bam
    test -s ${sample_id}.sorted.bam.bai
    """
}


process CONTIG_DEPTH {

    tag "${sample_id}"

    publishDir "${params.outdir}/04_binning/02_depth",
        mode: 'copy'

    conda "${projectDir}/envs/binning.yml"

    cpus 4
    memory '16 GB'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id),
          path("${sample_id}.depth.txt"),
          emit: depth

    script:
    """
    jgi_summarize_bam_contig_depths \
        --outputDepth ${sample_id}.depth.txt \
        ${bam}

    test -s ${sample_id}.depth.txt
    """
}


process METABAT2 {

    tag "${sample_id}"

    publishDir "${params.outdir}/04_binning/03_bins/${sample_id}",
        mode: 'copy'

    conda "${projectDir}/envs/binning.yml"

    cpus 16
    memory '48 GB'

    input:
    tuple val(sample_id), path(assembly), path(depth)
    val min_contig

    output:
    tuple val(sample_id),
          path("${sample_id}_bins"),
          emit: bins

    script:
    """
    mkdir -p ${sample_id}_bins

    metabat2 \
        -i ${assembly} \
        -a ${depth} \
        -o ${sample_id}_bins/${sample_id}.bin \
        -m ${min_contig} \
        -t ${task.cpus} \
        -v \
        2> ${sample_id}.metabat.log

    BIN_COUNT=\$(find ${sample_id}_bins \
        -type f -name "*.fa" | wc -l)

    echo "Number of bins: \$BIN_COUNT" \
        >> ${sample_id}.metabat.log
    """
}


process BIN_STATS {

    tag "${sample_id}"

    publishDir "${params.outdir}/04_binning/04_bin_stats",
        mode: 'copy'

    conda "${projectDir}/envs/binning.yml"

    cpus 2
    memory '8 GB'

    input:
    tuple val(sample_id), path(bin_dir)

    output:
    tuple val(sample_id),
          path("${sample_id}.bin_stats.tsv"),
          emit: stats

    script:
    """
    echo -e "bin\\tlength\\tsequences" \
        > ${sample_id}.bin_stats.tsv

    for BIN in ${bin_dir}/*.fa
    do
        [ -e "\$BIN" ] || continue

        NAME=\$(basename "\$BIN")

        seqkit stats \
            -T "\$BIN" \
        | awk -v name="\$NAME" \
            'NR==2 {print name "\\t" \$5 "\\t" \$4}' \
        >> ${sample_id}.bin_stats.tsv
    done
    """
}
