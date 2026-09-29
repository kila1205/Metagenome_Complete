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

    path "${sample_id}.emapper.log",
         emit: log,
         optional: true

    script:
    """
    # eggNOG-mapper 2.1.12 expects DIAMOND inside its package bin
    # directory, while Conda installs DIAMOND in \$CONDA_PREFIX/bin.
    # Create the expected link before running eggNOG-mapper.

    EGGNOG_PKG_BIN=\$(python - <<'PY'
import os
import eggnogmapper
print(os.path.join(os.path.dirname(eggnogmapper.__file__), "bin"))
PY
)

    mkdir -p "\$EGGNOG_PKG_BIN"
    ln -sf "\$CONDA_PREFIX/bin/diamond" "\$EGGNOG_PKG_BIN/diamond"

    echo "========================================"
    echo "eggNOG environment check"
    echo "========================================"
    python --version
    echo "emapper: \$(command -v emapper.py)"
    echo "diamond: \$(command -v diamond)"
    diamond --version

    emapper.py \
        -i ${proteins} \
        --itype proteins \
        -m diamond \
        --cpu ${task.cpus} \
        --data_dir ${eggnog_data_dir} \
        --output ${sample_id}

    test -s ${sample_id}.emapper.annotations
    test -s ${sample_id}.emapper.seed_orthologs
    """
}
