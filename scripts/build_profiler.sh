#!/bin/bash

sudo docker build -t rnaseq-umi-cpp -f Dockerfile.build  ${PWD}
sudo docker run --rm -v ${PWD}:/local rnaseq-umi-cpp /bin/sh -c "cp -r source/w* /local/. && cp /usr/local/bin/bwa /local/."
sudo docker build -t biodepot/rnaseq-umi-cpp:profiler  -f Dockerfile-profiler  ${PWD}
sudo rm -rf w96 w384
