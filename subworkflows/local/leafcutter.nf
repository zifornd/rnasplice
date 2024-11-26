include { REGTOOLS_JUNCTIONSEXTRACT } from '../../modules/nf-core/regtools/junctionsextract/main'
include { LEAFCUTTER_CLUSTER } from '../../modules/local/leafcutter_cluster'

workflow LEAFCUTTER{

    take:
    ch_genome_bam
    ch_genome_bam_index
    ch_gtf

    main:

    ch_versions = Channel.empty()

    REGTOOLS_JUNCTIONSEXTRACT(ch_genome_bam.join(ch_genome_bam_index))
    ch_versions = ch_versions.mix(REGTOOLS_JUNCTIONSEXTRACT.out.versions)
    ch_juncs = REGTOOLS_JUNCTIONSEXTRACT.out.junc
        .map { it[1] }
        .collect()

    LEAFCUTTER_CLUSTER(ch_juncs, ch_gtf)
    ch_versions = ch_versions.mix(LEAFCUTTER_CLUSTER.out.versions)

    emit:
    juncs                   = REGTOOLS_JUNCTIONSEXTRACT.out.junc.collect()
    counts                  = LEAFCUTTER_CLUSTER.out.counts
    perind_couts            = LEAFCUTTER_CLUSTER.out.perind_counts
    versions                = ch_versions

}
