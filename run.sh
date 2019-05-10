#!/bin/bash

set -e

nextflow \
    run \
    main.nf \
    -with-docker ubuntu:16.04 \
    -process.executor 'local' \
    -work-dir work/ \
    -resume
