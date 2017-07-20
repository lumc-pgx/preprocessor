#! /bin/bash

OUTPUT_DIR=""
EXTRA_CONFIG=""
CONFIG_FILE=""

while getopts ":d:c:e:" opt; do
    case $opt in
        d)
            OUTPUT_DIR="--directory $OPTARG"
            echo "Writing pipeline output to $OPTARG" >&2
            ;;
        c)
            CONFIG_FILE="--configfile $OPTARG"
            echo "Using additional config file $OPTARG"
            ;;
        e)
            EXTRA_CONFIG="--config ${OPTARG[@]}"
            echo "Using config options ${OPTARG[@]}" >&2
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo "Valid options are:" >&2
            echo "  -d : output directory" >&2
            echo "  -c : config file" >&2
            echo "  -e : extra config values" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

snakemake $OUTPUT_DIR \
          $EXTRA_CONFIG \
          $CONFIG_FILE \
          --latency-wait 90 \
          --drmaa ' -N preprocessor -pe BWA {cluster.threads} -l h_vmem={cluster.vmem} -q all.q -cwd -V -j Y' \
          --drmaa-log-dir cluster_logs \
          --jobs 100 \
          --max-jobs-per-second 10 \
          --cluster-config cluster_settings.yaml

