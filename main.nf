nextflow.enable.dsl = 2


/*
========================================================
 MODULES
========================================================
*/

include {
    RAW_NANOPLOT;
    FILTLONG;
    FILTERED_NANOPLOT;
    MULTIQC
} from './modules/qc/main'

include {
    METAFLYE;
    ASSEMBLY_STATS;
    QUAST
} from './modules/assembly/main'

include {
    PRODIGAL;
    PROTEIN_STATS
} from './modules/gene_prediction/main'

include {
    MAP_READS;
    CONTIG_DEPTH;
    METABAT2;
    BIN_STATS
} from './modules/binning/main'

include {
    COMBINE_PROTEINS;
    MMSEQS_CLUSTER;
    CATALOG_STATS
} from './modules/clustering/main'

include {
    COVERM_ABUNDANCE
} from './modules/abundance/main'

include {
    KRAKEN2_CLASSIFY;
    BRACKEN_ABUNDANCE
} from './modules/read_taxonomy/kraken_bracken'

include {
    GTDBTK_CLASSIFY
} from './modules/taxonomy/gtdbtk_classify'

include {
    PARSE_TAXONOMY
} from './modules/taxonomy/parse_taxonomy'

include {
    KOFAMSCAN
} from './modules/functional_prediction/kofamscan'

include {
    EGGNOG_MAPPER
} from './modules/functional_annotation/eggnog'

include {
    GENERATE_REPORT
} from './modules/report/main'


/*
========================================================
 PARAMETERS
========================================================
*/

params.input     = null
params.sample_id = "Metagenome"
params.outdir    = "${projectDir}/results"

params.min_contig = 1500

params.min_seq_id = 0.90
params.coverage   = 0.80

params.taxonomy_parser =
    "${projectDir}/scripts/parse_gtdbtk_taxonomy.py"


/*
--------------------------------------------------------
Kraken2 + Bracken
--------------------------------------------------------
*/

params.kraken_db =
    params.kraken_db ?: null

params.bracken_read_length = 300
params.bracken_level       = "S"
params.bracken_threshold   = 10


/*
--------------------------------------------------------
KofamScan
--------------------------------------------------------
*/

params.kofam_profiles =
    params.kofam_profiles ?: null

params.kofam_ko_list =
    params.kofam_ko_list ?: null


/*
--------------------------------------------------------
eggNOG-mapper
--------------------------------------------------------
*/

params.eggnog_data_dir =
    params.eggnog_data_dir ?: null


/*
--------------------------------------------------------
GTDB-Tk
--------------------------------------------------------
*/

params.gtdbtk_data_path =
    params.gtdbtk_data_path ?: null


/*
--------------------------------------------------------
Automatic report
--------------------------------------------------------
*/

params.report_script =
    "${projectDir}/scripts/report/generate_metagenome_report.py"


/*
========================================================
 WORKFLOW
========================================================
*/

workflow {

    /*
    ====================================================
    INPUT
    ====================================================
    */

    if (!params.input) {
        error "Missing required parameter: --input"
    }

    input_file = file(
        params.input,
        checkIfExists: true
    )

    /*
     * The input is already one merged biological sample.
     *
     * Channel structure:
     *
     *   sample_id, FASTQ
     */

    raw_reads = Channel.value(
        tuple(
            params.sample_id,
            input_file
        )
    )


    /*
    ====================================================
    READ QC + FILTERING
    ====================================================
    */

    RAW_NANOPLOT(
        raw_reads
    )

    FILTLONG(
        RAW_NANOPLOT.out.qc
    )

    FILTERED_NANOPLOT(
        FILTLONG.out.filtered
    )


    qc_dirs = RAW_NANOPLOT.out.qc
        .map { sample_id, reads, qc_dir ->
            qc_dir
        }
        .mix(
            FILTERED_NANOPLOT.out.qc
                .map { sample_id, reads, qc_dir ->
                    qc_dir
                }
        )
        .collect()


    MULTIQC(
        qc_dirs
    )


    /*
    ====================================================
    STEP 4
    READ-BASED TAXONOMY
    Kraken2 -> Bracken
    ====================================================
    */

    kraken_db_ch = Channel.value(
        file(
            params.kraken_db,
            checkIfExists: true
        )
    )


    KRAKEN2_CLASSIFY(
        FILTLONG.out.filtered,
        kraken_db_ch
    )


    BRACKEN_ABUNDANCE(
        KRAKEN2_CLASSIFY.out.report,
        kraken_db_ch,
        params.bracken_read_length,
        params.bracken_level,
        params.bracken_threshold
    )


    /*
    ====================================================
    STEP 1
    METAGENOME ASSEMBLY
    ====================================================
    */

    METAFLYE(
        FILTLONG.out.filtered
    )


    assembly_for_stats = METAFLYE.out.assembly
        .map { sample_id, assembly ->
            tuple(sample_id, assembly)
        }


    assembly_for_quast = METAFLYE.out.assembly
        .map { sample_id, assembly ->
            tuple(sample_id, assembly)
        }


    assembly_for_genes = METAFLYE.out.assembly
        .map { sample_id, assembly ->
            tuple(sample_id, assembly)
        }


    ASSEMBLY_STATS(
        assembly_for_stats
    )

    QUAST(
        assembly_for_quast
    )


    /*
    ====================================================
    STEP 2
    GENE PREDICTION
    ====================================================
    */

    PRODIGAL(
        assembly_for_genes
    )


    PROTEIN_STATS(
        PRODIGAL.out.predictions
    )


    /*
    ====================================================
    STEP 6
    PROTEIN CLUSTERING
    ====================================================
    */

    proteins_by_sample = PRODIGAL.out.predictions
        .map { sample_id, gff, genes, proteins ->
            tuple(
                sample_id,
                proteins
            )
        }
        .groupTuple()


    COMBINE_PROTEINS(
        proteins_by_sample
    )


    MMSEQS_CLUSTER(
        COMBINE_PROTEINS.out.combined,
        params.min_seq_id,
        params.coverage
    )


    CATALOG_STATS(
        MMSEQS_CLUSTER.out.representatives
    )


    /*
    ====================================================
    STEP 9
    FUNCTIONAL PREDICTION
    KofamScan
    ====================================================
    */

    kofam_profiles_ch = Channel.value(
        file(
            params.kofam_profiles,
            checkIfExists: true
        )
    )


    kofam_ko_list_ch = Channel.value(
        file(
            params.kofam_ko_list,
            checkIfExists: true
        )
    )


    KOFAMSCAN(
        MMSEQS_CLUSTER.out.representatives,
        kofam_profiles_ch,
        kofam_ko_list_ch
    )


    /*
    ====================================================
    STEP 10
    LARGE-SCALE FUNCTIONAL ANNOTATION
    eggNOG-mapper
    ====================================================
    */

    eggnog_data_ch = Channel.value(
        file(
            params.eggnog_data_dir,
            checkIfExists: true
        )
    )


    EGGNOG_MAPPER(
        MMSEQS_CLUSTER.out.representatives,
        eggnog_data_ch
    )


    /*
    ====================================================
    STEP 3
    BINNING
    ====================================================
    */

    assembly_for_binning = METAFLYE.out.assembly
        .map { sample_id, assembly ->
            tuple(sample_id, assembly)
        }


    reads_for_binning = FILTLONG.out.filtered
        .map { sample_id, reads ->
            tuple(sample_id, reads)
        }


    binning_input = assembly_for_binning
        .join(reads_for_binning)
        .map { sample_id, assembly, reads ->
            tuple(
                sample_id,
                assembly,
                reads
            )
        }


    MAP_READS(
        binning_input
    )


    CONTIG_DEPTH(
        MAP_READS.out.mapped
    )


    metabat_input = CONTIG_DEPTH.out.depth
        .combine(
            assembly_for_binning,
            by: 0
        )
        .map { sample_id, depth, assembly ->
            tuple(
                sample_id,
                assembly,
                depth
            )
        }


    METABAT2(
        metabat_input,
        params.min_contig
    )


    BIN_STATS(
        METABAT2.out.bins
    )


    /*
    ====================================================
    STEP 5
    BIN ABUNDANCE
    ====================================================
    */

    abundance_reads = FILTLONG.out.filtered
        .map { sample_id, reads ->
            tuple(
                sample_id,
                reads
            )
        }


    abundance_bins = METABAT2.out.bins
        .map { sample_id, bin_dir ->
            tuple(
                sample_id,
                bin_dir
            )
        }


    abundance_input = abundance_reads
        .join(abundance_bins)
        .map { sample_id, reads, bin_dir ->
            tuple(
                sample_id,
                reads,
                bin_dir
            )
        }


    COVERM_ABUNDANCE(
        abundance_input
    )


    /*
    ====================================================
    STEP 7
    GENOME-RESOLVED TAXONOMY
    GTDB-Tk
    ====================================================
    */

    gtdbtk_data_ch = Channel.value(
        file(
            params.gtdbtk_data_path,
            checkIfExists: true
        )
    )

    GTDBTK_CLASSIFY(
        METABAT2.out.bins,
        gtdbtk_data_ch
    )


    taxonomy_summary =
        GTDBTK_CLASSIFY.out.gtdbtk_results
        .map { sample_id, result_dir ->

            def summary = file(
                "${result_dir}/classify/gtdbtk.bac120.summary.tsv",
                checkIfExists: true
            )

            tuple(
                sample_id,
                summary
            )
        }


    taxonomy_parser_ch = Channel.value(
        file(
            params.taxonomy_parser,
            checkIfExists: true
        )
    )


    PARSE_TAXONOMY(
        taxonomy_summary,
        taxonomy_parser_ch
    )


    /*
    ====================================================
    STEP 8
    STATISTICAL / ECOLOGICAL ANALYSIS
    ====================================================

    The current dataset contains one biological sample.

    Therefore comparative analyses requiring multiple
    biological samples are intentionally not executed:

      - differential abundance testing
      - between-sample beta diversity
      - PERMANOVA
      - ecological association networks

    The ecology workflow can be applied later when
    multiple biological samples are available.

    ====================================================
    */


    /*
    ====================================================
    AUTOMATIC FINAL REPORT
    ====================================================
    */


    raw_nanostats_ch = RAW_NANOPLOT.out.qc
        .map { sample_id, reads, qc_dir ->
            tuple(
                sample_id,
                file(
                    "${qc_dir}/NanoStats.txt",
                    checkIfExists: true
                )
            )
        }


    filtered_nanostats_ch = FILTERED_NANOPLOT.out.qc
        .map { sample_id, reads, qc_dir ->
            tuple(
                sample_id,
                file(
                    "${qc_dir}/NanoStats.txt",
                    checkIfExists: true
                )
            )
        }


    assembly_stats_ch = ASSEMBLY_STATS.out
        .map { stats ->

            def name =
                stats.baseName
                    .replaceFirst(/\.assembly_stats$/, "")

            tuple(
                name,
                stats
            )
        }


    quast_ch = QUAST.out
        .map { quast_dir ->

            def name =
                quast_dir.baseName
                    .replaceFirst(/_quast$/, "")

            tuple(
                name,
                quast_dir
            )
        }


    protein_stats_ch = PROTEIN_STATS.out
        .map { stats ->

            def name =
                stats.baseName
                    .replaceFirst(/\.protein_stats$/, "")

            tuple(
                name,
                stats
            )
        }


    /*
     * Combine all result streams using sample ID.
     */

    report_inputs = raw_nanostats_ch
        .join(filtered_nanostats_ch)
        .join(assembly_stats_ch)
        .join(quast_ch)
        .join(protein_stats_ch)
        .join(BIN_STATS.out.stats)
        .join(COVERM_ABUNDANCE.out.abundance)
        .join(KRAKEN2_CLASSIFY.out.report)
        .join(BRACKEN_ABUNDANCE.out.abundance)
        .join(PARSE_TAXONOMY.out.taxonomy_table)
        .join(CATALOG_STATS.out.stats)
        .join(KOFAMSCAN.out.annotations)
        .join(EGGNOG_MAPPER.out.annotations)


    report_script_ch = Channel.value(
        file(
            params.report_script,
            checkIfExists: true
        )
    )


    GENERATE_REPORT(
        report_inputs,
        report_script_ch
    )
}
