snakemake --latency-wait 90 \
          --drmaa ' -N Preprocessor -pe BWA {cluster.threads} -l h_vmem={cluster.vmem} -cwd -V -e cluster_logs -o cluster_logs' \
          --jobs 100 \
          --max-jobs-per-second 10 \
          --cluster-config cluster_settings.yaml \
          --directory pipeline_output


