#!/bin/bash

docker_build_file=Dockerfile.build_ubuntu
docker_build_bwa_file=Dockerfile.build_bwa_ubuntu
temp_rna_seq_name=rnaseq-umi-cpp
temp_bwa="bwa:temp"

#build rna_seq executables
echo "sudo docker build -t $temp_rna_seq_name -f $docker_build_file  ${PWD}"
sudo docker build -t $temp_rna_seq_name -f $docker_build_file  ${PWD}
#copy rna_seq executables
echo 'sudo docker run --rm -v ${PWD}:/local  $temp_rna_seq_name /bin/sh -c "cp -r source/w* /local/."'
sudo docker run --rm -v ${PWD}:/local  $temp_rna_seq_name /bin/sh -c "cp -r source/w* /local/."
#build correct version of bwa
echo "sudo docker build -t $temp_bwa -f $docker_build_bwa_file  ${PWD}"
sudo docker build -t $temp_bwa -f $docker_build_bwa_file  ${PWD}
echo 'sudo docker run --rm -v ${PWD}:/local  $temp_bwa /bin/sh -c "cp -r bwa /local/bwa"'
sudo docker run --rm -v ${PWD}:/local  $temp_bwa /bin/sh -c "cp -r bwa /local/bwa"

#build final container
sudo docker build -t biodepot/rnaseq-umi-cpp:ubuntu_18.04__bwa_0.7.15 -f Dockerfile-ubuntu  ${PWD}
#build profiler container
sudo docker build -t biodepot/rnaseq-umi-cpp:profiler__ubuntu_18.04__bwa_0.7.15 -f Dockerfile-ubuntu-profiler  ${PWD}
sudo rm -rf w96 w384
