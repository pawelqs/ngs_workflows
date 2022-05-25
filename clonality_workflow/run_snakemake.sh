mkdir -p outputs


JOB_NAME=$(grep name config/config.yml | cut -f2 -d" ")

snakemake -j 100 \
    -k --use-conda \
    --configfile config/config.yml \
    --cluster-config workflow/cluster.yml \
    --cluster="sbatch -A {cluster.account} --mem={cluster.mem} --time={cluster.time} --cpus-per-task={cluster.cpus} --partition={cluster.partition} --job-name={cluster.job_name} --output={cluster.output}"'
