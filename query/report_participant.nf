N = 0
Input_ch = Channel.fromPath(params.patients_csv).splitCsv(header: true).map { it -> [N++, it] }
process report {
    tag "${row.family_id}"
    maxForks 25
    memory '2G'
    errorStrategy "retry"
    maxRetries 2
    clusterOptions "${params.clusterOptionsBasic}"
    input:
        set val(ind), val(row) from Input_ch
    output:
        val(row.family_id) into Report_ch
    script:
        ssh_socket = "${ind}.socket"
        local_port = params.mysql_tunnel_port + ind%10000
        outdir = "${params.outdir}/${row.family_id}"
    """
    set   +eu
    source ~/.bashrc
    set -eu
    source ${params.utilsh_path}
    echo ${HOSTNAME}
    # establish ssh tunnel, and get port
    coproc CAT { cat; }
    ssh_tunnel ${local_port} ${params.ssh_key} ${params.hostname} ${params.hostport} ${ssh_socket} >&\${CAT[1]}
    read newport <&\${CAT[0]}
    kill \$CAT_PID


    error=0
    {
        python ${baseDir}/report_participant.py \
            --family_id ${row.family_id} \
            --outdir ${outdir} \
            --config ${baseDir}/query_config.yml \
            --mysqlyaml ${baseDir}/mysql.yml \
            --port \${newport} 
    } || {
        echo "\$? script error"
        error=1
    }

    ssh -S ${ssh_socket} -O exit ${params.hostname}

    exit \$error

    """
}
/*
process igv_plot {
    tag "${family_id}"
    clusterOptions "${params.clusterOptionsBasic}"
    module "bio/IGV/2.5.0-Java-11"
    memory "8G"
    input:
        val(family_id) from Report_ch
    script:
        wd = "${params.outdir}/${family_id}"
    """
    set   +eu
    source ~/.bashrc
    set -eu
    if [ -s ${wd}/igv_batch.txt ]; then
        xvfb-run --auto-servernum igv.sh \
            -g /public_data_resources/IGV/hg38/hg38_local.genome \
            -b ${wd}/igv_batch.txt || exit 0
    fi
    """
}
*/
