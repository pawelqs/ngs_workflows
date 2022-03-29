# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/pkus/programs/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/pkus/programs/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/pkus/programs/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/pkus/programs/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate snakemake
mkdir -p outputs

#cluster_command="sbatch\
# -A {cluster.account}\
# --mem={cluster.mem}\
# --time={cluster.time}\
# --cpus-per-task={cluster.cpus}\
# --partition={cluster.partition}\
# --job-name={cluster.job_name}\
# --output={cluster.output}"

JOB_NAME=$(grep name config/config.yml | cut -f2 -d" ")

sbatch \
    --output=outputs/snakemake.out \
    --partition=ibm_small \
    --job-name=$JOB_NAME \
    --wrap='snakemake -j 50 -k --use-conda --cluster-config config/cluster.yml --cluster="sbatch -A {cluster.account} --mem={cluster.mem} --time={cluster.time} --cpus-per-task={cluster.cpus} --partition={cluster.partition} --job-name={cluster.job_name} --output={cluster.output}"'
