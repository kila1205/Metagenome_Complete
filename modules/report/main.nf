process GENERATE_REPORT {

    tag "${sample_id}"

    publishDir "${params.outdir}/report",
        mode: 'copy'

    conda "${projectDir}/envs/report.yml"

    cpus 2
    memory '8 GB'
    time '2h'

    input:
    tuple val(sample_id),
          // stageAs: NanoPlot SELALU menamai keluarannya NanoStats.txt, jadi
          // statistik sebelum dan sesudah penyaringan punya nama sama dan
          // bertabrakan saat dipentaskan ke satu direktori kerja:
          //   "input file name collision -- multiple input files for each of
          //    the following file names: NanoStats.txt"  (25 Sep 2026)
          // Diberi nama berbeda saat dipentaskan. Skrip di bawah tetap memakai
          // ${raw_nanostats} dan ${filtered_nanostats}; Nextflow otomatis
          // menggantinya dengan nama baru ini.
          path(raw_nanostats, stageAs: 'raw_NanoStats.txt'),
          path(filtered_nanostats, stageAs: 'filtered_NanoStats.txt'),
          path(assembly_stats),
          path(quast_report),
          path(protein_stats),
          path(bin_stats),
          path(coverm),
          path(kraken_report),
          path(bracken),
          path(taxonomy),
          path(catalog_stats),
          path(kofam),
          path(eggnog)

    path report_script

    output:
    path "${sample_id}_metagenome_report.pdf",
         emit: pdf

    path "${sample_id}_metagenome_summary.tsv",
         emit: summary

    path "${sample_id}_report_figures",
         emit: figures

    script:
    """
    python ${report_script} \
        --sample-id "${sample_id}" \
        --raw-nanostats "${raw_nanostats}" \
        --filtered-nanostats "${filtered_nanostats}" \
        --assembly-stats "${assembly_stats}" \
        --quast-report "${quast_report}" \
        --protein-stats "${protein_stats}" \
        --bin-stats "${bin_stats}" \
        --coverm "${coverm}" \
        --kraken-report "${kraken_report}" \
        --bracken "${bracken}" \
        --taxonomy "${taxonomy}" \
        --catalog-stats "${catalog_stats}" \
        --kofam "${kofam}" \
        --eggnog "${eggnog}" \
        --output-pdf "${sample_id}_metagenome_report.pdf" \
        --output-summary "${sample_id}_metagenome_summary.tsv" \
        --figures-dir "${sample_id}_report_figures"

    test -s ${sample_id}_metagenome_report.pdf
    test -s ${sample_id}_metagenome_summary.tsv
    """
}
