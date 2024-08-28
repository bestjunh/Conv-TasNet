#!/bin/bash

# 특정 디렉토리로 이동
cd  /usr/local/MATLAB/R2023b/bin  

# MATLAB 스크립트를 여러 번 실행
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 
./matlab -nodisplay -nodesktop -r "cd /home/albert/Albert-Conv-TasNet/create-mixture/;rng('shuffle');run('create_wav_2speakers_rir.m'); exit;" & 

# 모든 백그라운드 작업이 완료될 때까지 기다림
wait

# 완료 메시지 출력
echo "All MATLAB scripts have been executed."
