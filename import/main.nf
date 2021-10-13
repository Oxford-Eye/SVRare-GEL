// This goes from creating the tables to prepare files for database import
// Import can't be done in parallel.
process 'create_tables' {
    clusterOptions "${params.clusterOptionsBasic}"
    output:
        val(create_done) into create_tables_ch
    script:
        create_done = true
        ssh_socket = 'sshs.socket'
    """
    set +eu
    source ~/.bashrc
    set -eu
    source ${baseDir}/bin/utils.sh

    # establish ssh tunnel, and get port
    coproc CAT { cat; }
    ssh_tunnel ${params.mysql_tunnel_port} ${params.ssh_key} ${params.hostname} ${params.hostport} ${ssh_socket} >&\${CAT[1]}
    read newport <&\${CAT[0]}
    kill \$CAT_PID

    error=0
    {
        python ${baseDir}/create_tables.py --port \${newport} --mysqlyaml ${baseDir}/mysql.yml
    } || {
        echo "\$? script error"
        error=1
    }
    ssh -S ${ssh_socket} -O exit ${params.hostname}

    exit \$error
    
    """
}

process 'import_gene' {
    clusterOptions "${params.clusterOptionsBasic}"
    input:
        val(create_table_done) from create_tables_ch
    output:
        val(import_gene_done) into import_gene_ch
        file("${outfile_base}.gz") into gtf_ch
    script:
        ssh_socket = "ssh.socket"
        import_gene_done = true
        outfile_base = "sorted.gtf"
    
    """
    set  +eu
    source ~/.bashrc
    set -eu
    source ${baseDir}/bin/utils.sh

    # establish ssh tunnel, and get port
    coproc CAT { cat; }
    ssh_tunnel ${params.mysql_tunnel_port} ${params.ssh_key} ${params.hostname} ${params.hostport} ${ssh_socket} >&\${CAT[1]}
    read newport <&\${CAT[0]}
    kill \$CAT_PID

    error=0
    {
        python ${baseDir}/import_genes.py \
            --gtf ${params.genes_gtf} \
            --port \${newport} \
            --mysqlyaml ${baseDir}/mysql.yml \
            --output ${outfile_base}
        bgzip ${outfile_base} && tabix -s 1 -b 2 -e 3 ${outfile_base}.gz
        cp ${outfile_base}.gz ${params.genes_gtf_converted}
        cp ${outfile_base}.gz.tbi ${params.genes_gtf_converted}.tbi
    } || {
        echo "\$? script error"
        error=1
    }
    ssh -S ${ssh_socket} -O exit ${params.hostname}

    exit \$error
    """
}

Patient_file_ch = Channel.fromPath(params.patient_table)

process 'import_patients' {
    clusterOptions "${params.clusterOptionsBasic}"
    input:
        file(patient_file) from Patient_file_ch
        val(create_db_done) from create_tables_ch
    output:
        val(import_patients_done) into Import_patients_ch
    script:
        import_patients_done = true
        ssh_socket = "ssh.socket"
    """
    set +eu
    source ~/.bashrc
    set -eu
    source ${baseDir}/bin/utils.sh

    # establish ssh tunnel, and get port
    coproc CAT { cat; }
    ssh_tunnel ${params.mysql_tunnel_port} ${params.ssh_key} ${params.hostname} ${params.hostport} ${ssh_socket} >&\${CAT[1]}
    read newport <&\${CAT[0]}
    kill \$CAT_PID

    error=0
    {
        python ${baseDir}/import_patients.py --input ${params.patient_table} --port \${newport} --mysqlyaml ${baseDir}/mysql.yml
    } || {
        echo "\$? script error"
        error=1
    }
    ssh -S ${ssh_socket} -O exit ${params.hostname}

    exit \$error

    """
}

vcf_to_tsv_input_ch = Channel.fromPath(params.patient_table)
    .splitCsv(header: true, sep: '\t')
process 'vcf_to_tsv' {
    tag "${row.participant_id}"
    queue 'short'
    //memory "10G"
    clusterOptions "${params.clusterOptionsBasic}"
    publishDir params.vcf_to_tsv_path, mode: "copy"
    input:
        val(row) from vcf_to_tsv_input_ch
        val(import_gene_done) from import_gene_ch
    output:
        set file("${outfile}.gz"), file("${outfile}.gz.tbi") into vcf_to_tsv_ch
    when:
        row.genome == 'Available'
    script:
        ssh_socket = "ssh.socket"
        outfile = "${row.participant_id}_${row.plate_key}.tsv"
        // don't change outfile definition, since it will be used to figure out participant_id and plate_key downstream
    """
    set +eu
    source ~/.bashrc
    set -eu
    source ${baseDir}/bin/utils.sh

    input=${row.path}/${row.plate_key}/Variations/${row.plate_key}.SV.vcf.gz
    python ${baseDir}/vcf_to_tsv.py \
        --input \$input \
        --gtf ${params.genes_gtf_converted} \
        --build ${row.genome_build} \
        --padding ${params.bnd_genes_padding} \
        --liftOver_chainfile ${params.liftOver_chainfile} >temp
    grep '^#' temp >${outfile}
    grep -v '^#' temp | sort -k1,1 -k2,2n -k3,3n >>${outfile} && bgzip ${outfile} && tabix -s 1 -b 2 -e 3 ${outfile}.gz
    rm temp

    """
}

Chromosome_ch = Channel.from(
    'chr1',
    'chr2',
    'chr3',
    'chr4',
    'chr5',
    'chr6',
    'chr7',
    'chr8',
    'chr9',
    'chr10',
    'chr11',
    'chr12',
    'chr13',
    'chr14',
    'chr15',
    'chr16',
    'chr17',
    'chr18',
    'chr19',
    'chr20',
    'chr21',
    'chr22',
    'chrX',
    'chrY',
    'chrM'
).combine(Channel.from('INV','GAIN','LOSS'))

process walk_intervals {
    tag "${chrom}-${type}"
    queue "short"
    clusterOptions "${params.clusterOptionsBasic}"
    publishDir "${baseDir}/walk_intervals", mode: "copy"
    input:
        set val(chrom), val(type) from Chromosome_ch
    output:
        file(output) into Walk_intervals_ch
    //when:
    //    chrom == 'chr1' && type == 'LOSS'
    script:
        output = "${chrom}.${type}.count.tsv"
    """
    set   +eu
    source ~/.bashrc
    set -eu
    source ${baseDir}/bin/utils.sh

    echo ${HOSTNAME}

    python ${baseDir}/walk_intervals.py \
        --chrom ${chrom} \
        --input_path ${params.vcf_to_tsv_path} \
        --type ${type} \
        --output ${output}
    """

}

N = 0
Chromosome_ch = Channel.fromPath("${baseDir}/walk_intervals/*.tsv").map { it -> [N++, it] }
//Chromosome_ch = Channel.fromPath("${baseDir}/test/*.tsv").map { it -> [N++, it] }

process 'make_files' {
    tag "${chrom} ${interval_type}"
    queue "short"
    memory "8G"
    errorStrategy "finish"
    //publishDir "${baseDir}/make_files", mode: "copy"
    clusterOptions "${params.clusterOptionsBasic}"
    input:
        set val(ind), file(input_file) from Chromosome_ch
    output:
        set file(interval_file),
            file(interval_gene_file),
            file(patient_interval_file),
            file(group_file),
            file(interval_group_file) into Make_files_ch
    script:
        input_file_tokens = input_file.name.toString().tokenize('.')
        chrom = input_file_tokens.get(0)
        interval_type = input_file_tokens.get(1)
        interval_file = "${chrom}_${interval_type}.interval.tsv"
        interval_gene_file = "${chrom}_${interval_type}.interval_gene.tsv"
        patient_interval_file = "${chrom}_${interval_type}.patient_interval.tsv"
        group_file = "${chrom}_${interval_type}.group.tsv"
        interval_group_file = "${chrom}_${interval_type}.interval_group.tsv"
    """
    set   +eu
    source ~/.bashrc
    set -eu
    echo ${HOSTNAME}

    python ${baseDir}/import_intervals_make_files.py \
        --region_file ${input_file} \
        --interval_output ${interval_file} \
        --interval_gene_output ${interval_gene_file} \
        --patient_interval_output ${patient_interval_file} \
        --group_output ${group_file} \
        --interval_group_output ${interval_group_file} \

    """
}