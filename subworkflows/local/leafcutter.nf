include { REGTOOLS_JUNCTIONSEXTRACT } from '../modules/nf-core/regtools/junctionsextract/main'
include { LEAFCUTTER_CLUSTER } from '../modules/local/leafcutter_cluster'

workflow LEAFCUTTER{

    take:
    ch_genome_bam
    ch_genome_bam_index
    ch_samplesheet

    main:

    ch_versions = Channel.empty()

    REGTOOLS_JUNCTIONSEXTRACT(ch_samplesheet, ch_genome_bam, ch_genome_bam_index)
    ch_versions = ch_versions.mix(REGTOOLS_JUNCTIONSEXTRACT.out.versions)
    ch_juncs = REGTOOLS_JUNCTIONSEXTRACT.out.junc.collect()

    LEAFCUTTER_CLUSTER(ch_samplesheet, ch_juncs)
    ch_versions = ch_versions.mix(LEAFCUTTER_CLUSTER.out.versions)

    emit:
    juncs                   = REGTOOLS_JUNCTIONSEXTRACT.out.junc.collect()
    counts                  = LEAFCUTTER_CLUSTER.out.counts
    perind_couts            = LEAFCUTTER_CLUSTER.out.perind_counts
    cluster_significance    = LEAFCUTTER_CLUSTER.out.cluster_significance
    versions                = ch_versions

}
