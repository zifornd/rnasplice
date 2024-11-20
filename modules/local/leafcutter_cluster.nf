process LEAFCUTTER_CLUSTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::leafcutter=0.2.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/leafcutter:0.2.9--py37r41h399db7b_2' :
        'quay.io/biocontainers/leafcutter:0.2.9--py37r41h399db7b_2' }"

    input:
    tuple val(meta), path(junc_files)
    path(gtf)

    output:
    tuple val(meta), path("*_perind_numers.counts.gz"), emit: counts
    tuple val(meta), path("*_perind.counts.gz")       , emit: perind_counts
    tuple val(meta), path("*_cluster_significance.txt"), emit: cluster_significance
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    //TODO bash line to make txt file listing junctions files, feed this to python script
    """
    python $projectDir/bin/leafcutter_cluster_regtools.py \\
        $args \\
        -j ${junc_files.join(',')} \\
        -o ${prefix} \\
        -r ${gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leafcutter: \$(python -c "import leafcutter; print(leafcutter.__version__)")
    END_VERSIONS
    """
}
