#! /bin/bash

OUTPUT_DIR=pipeline_output

while getopts ":d:" opt; do
    case $opt in
        d)
            OUTPUT_DIR=$OPTARG
            echo "Writing pipeline output to $OUTPUT_DIR" >&2
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

snakemake --latency-wait 90 \
          --drmaa ' -N preprocessor -pe BWA {cluster.threads} -l h_vmem={cluster.vmem} -q all.q -cwd -V -j Y' \
          --drmaa-log-dir cluster_logs \
          --jobs 100 \
          --max-jobs-per-second 10 \
          --cluster-config cluster_settings.yaml \
          --directory $OUTPUT_DIR


