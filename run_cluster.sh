snakemake --latency-wait 90 \
          --drmaa ' -N preprocessor -pe BWA {cluster.threads} -l h_vmem={cluster.vmem} -q all.q -cwd -V -j Y' \
          --drmaa-log-dir cluster_logs \
          --jobs 100 \
          --max-jobs-per-second 10 \
          --cluster-config cluster_settings.yaml \
          --directory pipeline_output


