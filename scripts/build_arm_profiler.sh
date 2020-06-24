#!/bin/bash

echo  "sudo docker build -t rnaseq-umi-cpp -f Dockerfile.build_ARM64  ${PWD}"
sudo docker build -t rnaseq-umi-cpp -f Dockerfile.build_ARM64  ${PWD}

echo sudo docker run --rm -v ${PWD}:/local rnaseq-umi-cpp /bin/sh -c "cp -r source/w* /local/. "
sudo docker run --rm -v ${PWD}:/local rnaseq-umi-cpp /bin/sh -c "cp -r source/w* /local/. "

echo "sudo docker build -t bwa:build -f Dockerfile.build_bwa ."
sudo docker build -t bwa:build -f Dockerfile.build_bwa .

echo "sudo docker run --rm -v ${PWD}:/local bwa:build /bin/sh -c '"' cp -r /usr/local/bin/bwa /local/.'"'"
sudo docker run --rm -v ${PWD}:/local bwa:build /bin/sh -c "cp ./bwa /local/."

echo "sudo docker build -t biodepot/rnaseq-umi-cpp:profiler  -f Dockerfile-arm-profiler  ${PWD}"
sudo docker build -t biodepot/rnaseq-umi-cpp:profiler  -f Dockerfile-arm-profiler  ${PWD}

echo "sudo rm -rf w96 w384"
sudo rm -rf w96 w384 bwa
sudo docker rmi bwa:build
