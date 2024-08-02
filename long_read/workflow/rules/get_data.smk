# Rule rsync downloads the data from the corresponding remote server.
rule rsync:
    output:
        directory(str(SAMPLES_DIR) + "/{sample}/runs")
    log:
        str(SAMPLES_DIR) + "/{sample}/log/rsync.log"
    params:
        rsync_path=lambda wilds: RSYNC_PATH,
        date = ('[0-9]'*8),
    resources:
        runtime=2*24*60
    shell:
        """
        rsync -a --info=progress2 \
                {params.rsync_path}/{params.date}_{wildcards.sample}/runs/* \
                {output}/ \
                &> {log}
        """