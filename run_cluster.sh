snakemake --latency-wait 90 \
          --drmaa ' -N Preprocessor -pe BWA {cluster.threads} -l h_vmem={cluster.vmem} -cwd -V' \
          --drmaa-log-dir cluster_logs \
          --jobs 100 \
          --max-jobs-per-second 10 \
          --cluster-config cluster_settings.yaml \
          --directory pipeline_output


