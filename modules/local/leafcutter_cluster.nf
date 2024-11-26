process LEAFCUTTER_CLUSTER {
    label 'process_medium'

    conda "conda-forge::python=3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.7' :
        'quay.io/biocontainers/python:3.7' }"

    input:
    path junc_files
    path gtf

    output:
    path("*perind_numers.counts.gz") , emit: counts
    path("*perind.counts.gz")        , emit: perind_counts
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ls *junc > test_juncfiles.txt
    leafcutter_cluster_regtools.py \\
        $args \\
        -j test_juncfiles.txt \\
        -o lc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        leafcutter: 0.2.9)
    END_VERSIONS
    """
}
