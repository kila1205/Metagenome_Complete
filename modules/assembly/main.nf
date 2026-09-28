process METAFLYE {

    tag "${sample_id}"

    publishDir "${params.outdir}/02_assembly/01_metaflye/${sample_id}",
        mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("assembly.fasta"),
          emit: assembly

    path "assembly_info.txt",
         optional: true

    script:
    """
    flye \
        --nano-hq ${reads} \
        --meta \
        --threads ${task.cpus} \
        --out-dir flye_output

    cp flye_output/assembly.fasta assembly.fasta

    if [ -f flye_output/assembly_info.txt ]; then
        cp flye_output/assembly_info.txt assembly_info.txt
    fi

    test -s assembly.fasta
    """
}


process ASSEMBLY_STATS {

    tag "${sample_id}"

    publishDir "${params.outdir}/02_assembly/02_stats",
        mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path "${sample_id}.assembly_stats.txt"

    script:
    """
    seqkit stats -a -T ${assembly} \
        > ${sample_id}.assembly_stats.txt
    """
}


process QUAST {

    tag "${sample_id}"

    publishDir "${params.outdir}/02_assembly/03_quast/${sample_id}",
        mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path "${sample_id}_quast"

    script:
    """
    mkdir -p ${sample_id}_quast

    quast.py \
        ${assembly} \
        -o ${sample_id}_quast \
        -t ${task.cpus}
    """
}
