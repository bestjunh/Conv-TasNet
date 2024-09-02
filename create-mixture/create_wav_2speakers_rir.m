% create_wav_2_speakers.m
%
% Create 2-speaker mixtures
% 
% This script assumes that WSJ0's wv1 sphere files have already
% been converted to wav files, using the original folder structure
% under wsj0/, e.g., 
% 11-1.1/wsj0/si_tr_s/01t/01to030v.wv1 is converted to wav and 
% stored in YOUR_PATH/wsj0/si_tr_s/01t/01to030v.wav, and
% 11-6.1/wsj0/si_dt_05/050/050a0501.wv1 is converted to wav and
% stored in YOUR_PATH/wsj0/si_dt_05/050/050a0501.wav.
% Relevant data from all disks are assumed merged under YOUR_PATH/wsj0/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2016 Mitsubishi Electric Research Labs 
%                          (Jonathan Le Roux, John R. Hershey, Zhuo Chen)
%   Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("functions/")

if ispc
    nasPath = 'Z:/nas1_data/';
    devDataPath = 'T:/dev60-data-mount/';
elseif isunix
    nasPath = '/home/nas/';
    devDataPath = '/home/dev60-data-mount/';
else
    disp('Unknown operating system.');
end


data_type = {'tr','cv','tt'};
wsj0root = [nasPath 'DB/wsj_wav/']; % YOUR_PATH/, the folder containing wsj0/
noiseroot = [nasPath 'user/Uihyeop/DB/wham_noise/'];
rirroot = [devDataPath 'albert/DB/CONV-TASNet-RIR-v2/'];
output_dir16k=[devDataPath 'albert/DB/wsj0-mix_20240903_rir/2speakers_noise/wav16k/'];

min_max = {'min'};%{'min','max'};
upFs = 3*16000;
early_idx = floor(0.05 * upFs); % 50ms

nch = 1;
NNoise = 18;
RT60 = [0.2 0.3 0.4 0.5 0.6 0.7];

ROOMDIM = {[4.01, 3.65, 2.81]
           %[6.02, 3.82, 2.53]
           %[3.35, 4.94, 2.83]
           %[4.55, 5.76, 2.61]
           %[5.33, 3.51, 2.72]
           [4.34, 6.02, 2.55]
           [4.55, 3.54, 2.81]
           [4.76, 3.86, 2.74]
           [4.87, 3.78, 2.57]
           [5.08, 4.02, 2.64]
           [5.94, 7.45, 2.82]
           [7.40, 5.55, 2.87]
           [7.13, 5.78, 2.82]
           [7.23, 5.22, 2.91]
           %[8.12, 5.45, 3.01]
           [8.20, 5.55, 2.89]
           [6.63, 4.75, 2.84]
           [6.30, 6.15, 2.92]
           [8.40, 5.45, 2.95]
           [8.10, 6.42, 3.00]};

ROOMDIM_test = {
    [6.02, 3.82, 2.53]
    [3.35, 4.94, 2.83]
    [4.55, 5.76, 2.61]
    [5.33, 3.51, 2.72]
    [8.12, 5.45, 3.01]

    };
          

LOCDELTA = {[-0.52 -0.15,-0.11]
            [-0.35 -0.41 0.09]
            [-0.16 0.38 -0.10]
            [-0.45 0.20 0.12]
            [0.45 -0.18 -0.14]
            [0.55 -0.33 0.08]
            [0.31 0.22 0.02]
            [0.27 0.43 -0.03]};


for i_mm = 1:length(min_max)
    for i_type = 1:length(data_type)
        if ~exist([output_dir16k '/' min_max{i_mm} '/' data_type{i_type}],'dir')
            mkdir([output_dir16k '/' min_max{i_mm} '/' data_type{i_type}]);
        end

        status = mkdir([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s1/']);
        status = mkdir([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s2/']);
        status = mkdir([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/mix/']);
        status = mkdir([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s1_early/']);
        status = mkdir([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s2_early/']);
        status = mkdir([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s1_dry/']);
        status = mkdir([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s2_dry/']);
                
        TaskFile = ['mix_2_spk_' data_type{i_type} '.txt'];
        fid=fopen(TaskFile,'r');
        C=textscan(fid,'%s %f %s %f');
        
        Source1File = ['mix_2_spk_' min_max{i_mm} '_' data_type{i_type} '_1'];
        Source2File = ['mix_2_spk_' min_max{i_mm} '_' data_type{i_type} '_2'];
        MixFile     = ['mix_2_spk_' min_max{i_mm} '_' data_type{i_type} '_mix'];
        
        fid_s1 = fopen(Source1File,'w');
        fid_s2 = fopen(Source2File,'w');
        fid_m  = fopen(MixFile,'w');
        
        num_files = length(C{1});

        scaling_16k = zeros(num_files,2);        
        scaling16bit_16k = zeros(num_files,1);
        
        
        fprintf(1,'%s\n',[min_max{i_mm} '_' data_type{i_type}]);

        noiseWAVList = dir([noiseroot data_type{i_type}]);noiseWAVList = noiseWAVList(3:end);



        tic
        for i = randperm(num_files)%1:num_files
            
            [inwav1_dir,invwav1_name,inwav1_ext] = fileparts(C{1}{i});
            [inwav2_dir,invwav2_name,inwav2_ext] = fileparts(C{3}{i});
            fprintf(fid_s1,'%s\n',C{1}{i});
            fprintf(fid_s2,'%s\n',C{3}{i});
            inwav1_snr = C{2}(i);
            inwav2_snr = C{4}(i);
            mix_name = [invwav1_name,'_',num2str(inwav1_snr),'_',invwav2_name,'_',num2str(inwav2_snr)];
            fprintf(fid_m,'%s\n',mix_name);

            if exist([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s2_dry/' mix_name '.wav'],'file') == 2
                disp(['Already exist: ' num2str(i_mm) '/' num2str(length(min_max)) ' ' num2str(i_type) '/' num2str(length(data_type)) ' ' num2str(i) '/' num2str(num_files)]);
                continue
            else
                disp(['Creating: ' num2str(i_mm) '/' num2str(length(min_max)) ' ' num2str(i_type) '/' num2str(length(data_type)) ' ' num2str(i) '/' num2str(num_files)]);
            end

            % load RIRs
            
            rtIdx = randi(length(RT60));            
            micIdx = randi(length(LOCDELTA));
            if data_type{i_type} == 'tt'
                roomIdx = randi(length(ROOMDIM_test));
            else
                roomIdx = randi(length(ROOMDIM));
            end
            roomDim = cell2mat(ROOMDIM(roomIdx));
            locDelta = cell2mat(LOCDELTA(micIdx));
            targetRIRPath = [rirroot 'target/RT' num2str(RT60(rtIdx)) '/ROOM' num2str(roomDim(1)) 'x' num2str(roomDim(2)) 'x' num2str(roomDim(3)) '/LOC' num2str(locDelta(1)) 'x' num2str(locDelta(2)) 'x' num2str(locDelta(3)) '/'];
            noiseRIRPath = [rirroot 'noise/RT' num2str(RT60(rtIdx)) '/ROOM' num2str(roomDim(1)) 'x' num2str(roomDim(2)) 'x' num2str(roomDim(3)) '/LOC' num2str(locDelta(1)) 'x' num2str(locDelta(2)) 'x' num2str(locDelta(3)) '/'];
            target_rirList = dir(targetRIRPath);target_rirList = target_rirList(3:end);
            noiseRIR = dir(noiseRIRPath);noiseRIR = noiseRIR(3:end);
            

            rir1Idx = randi(length(target_rirList));
            while(1)
                rir2Idx = randi(length(target_rirList));
                if rir1Idx ~= rir2Idx
                    break
                end
            end
            
            h = load([targetRIRPath target_rirList(rir1Idx).name]);
            h1 = h.h;
            h = load([targetRIRPath target_rirList(rir2Idx).name]);
            h2 = h.h;

            h_noise = cell(1,18);
            minLength = Inf;
            maxVal = 0;
            for src = 1:NNoise
                htmp = load([noiseRIR(src).folder '/' noiseRIR(src).name]);
                htmp = htmp.h;
                h_noise{src} = htmp;
                if minLength > length(htmp)
                    minLength = length(htmp);
                end
                if maxVal < max(htmp)
                    maxval = max(htmp);
                end
            end

            minLength = min([length(h1),length(h2),minLength]);
            h1 = h1(1:minLength);
            h2 = h2(1:minLength);
            maxVal = max([max(h1),max(h2),maxval]);
            h1 = h1/maxVal;
            h2 = h2/maxVal;
            for src = 1:NNoise
                h_noise{src} = h_noise{src}(1:minLength);
                h_noise{src} = h_noise{src}/maxVal;
            end

            [~,idx] = max(h1);
            h1_early = h1;
            h1_early(idx+early_idx:end) = 0;

            [~,idx] = max(h2);
            h2_early = h2;
            h2_early(idx+early_idx:end) = 0;
                        
            % get input wavs
            [s1, fs] = audioread([wsj0root C{1}{i}]);
            s2       = audioread([wsj0root C{3}{i}]);

            cropWavLen = min(length(s1),length(s2));
            noise = zeros(cropWavLen,NNoise);
            for src = 1:NNoise
                while(1)
                    idx = randi(num_files);
                    [tmp, fsNoise] = audioread([noiseroot data_type{i_type} '/' noiseWAVList(idx).name]);
                    if length(tmp) < cropWavLen
                        continue
                    end
                    startSample = randi(length(tmp)-cropWavLen)+1;
                    noise(:,src) = tmp(startSample:startSample+cropWavLen-1,1);
                    break
                end
                
            end
            
            
            % Convolve with RIRs            
            s1_up = resample(s1,upFs,fs);
            s2_up = resample(s2,upFs,fs);
            s1_up_rir = fft_filter(h1,1,s1_up);
            s2_up_rir = fft_filter(h2,1,s2_up);
            s1_early_up_rir = fft_filter(h1_early,1,s1_up);
            s2_early_up_rir = fft_filter(h2_early,1,s2_up);
            
            % resample, normalize 16 kHz file, save scaling factor
            s1_early_16k_rir=resample(s1_early_up_rir,fs,upFs);
            s2_early_16k_rir=resample(s2_early_up_rir,fs,upFs);
            s1_16k_rir=resample(s1_up_rir,fs,upFs);
            s2_16k_rir=resample(s2_up_rir,fs,upFs);
            s1_16k = resample(s1,16000,fs);
            s2_16k = resample(s2,16000,fs);

            noise_up_rir = zeros(min(length(s1_up),length(s2_up)),NNoise);
            noise_16k_rir = zeros(min(length(s1_16k),length(s2_16k)),NNoise);
            for src = 1:NNoise
                tmp_up = resample(noise(:,src),upFs,fs);                
                noise_up_rir(:,src) = fft_filter(h_noise{src},1,tmp_up);
                noise_16k_rir(:,src) = resample(noise_up_rir(:,src),fs,upFs);
            end

            
            % [s1_16k_rir,lev1]=activlev(s1_16k_rir,16000,'n'); % y_norm = y /sqrt(lev);            
            % [s2_16k_rir,lev2]=activlev(s2_16k_rir,16000,'n');


            
            energyInterference = sqrt(sum(sum(abs(noise_16k_rir).^2)));
            
            energyTarget = sqrt(sum(abs(s1_16k_rir(1:length(noise_16k_rir))).^2));            
            energyNormal = energyTarget/energyInterference;
            scale1 = 1/(energyNormal/sqrt(10^(inwav1_snr/10)));
            

            energyTarget = sqrt(sum(abs(s2_16k_rir(1:length(noise_16k_rir))).^2));            
            energyNormal = energyTarget/energyInterference;
            scale2 = 1/(energyNormal/sqrt(10^(inwav2_snr/10)));
            
            s1_16k_rir = s1_16k_rir.*scale1;
            s2_16k_rir = s2_16k_rir.*scale2;
            s1_early_16k_rir = s1_early_16k_rir.*scale1;
            s2_early_16k_rir = s2_early_16k_rir.*scale2;
            s1_16k = s1_16k.*scale1;
            s2_16k = s2_16k.*scale2;

            % energyTarget = sqrt(sum(abs(s2_16k_rir(1:length(noise_16k_rir))).^2));
            % 20*log10(energyTarget/energyInterference)
            % debug = 1;
            
            
            switch min_max{i_mm}
                case 'max'
                    % mix_16k_length = max(length(s1_16k_rir),length(s2_16k_rir));
                    % s1_16k_rir = cat(1,s1_16k_rir,zeros(mix_16k_length - length(s1_16k_rir),1));
                    % s2_16k_rir = cat(1,s2_16k_rir,zeros(mix_16k_length - length(s2_16k_rir),1));
                    % s1_16k = cat(1,s1_16k,zeros(mix_16k_length - length(s1_16k),1));
                    % s2_16k = cat(1,s2_16k,zeros(mix_16k_length - length(s2_16k),1));
                case 'min'
                    mix_16k_length = min(length(s1_16k_rir),length(s2_16k_rir));
                    s1_16k_rir = s1_16k_rir(1:mix_16k_length);
                    s2_16k_rir = s2_16k_rir(1:mix_16k_length);
                    s1_early_16k_rir = s1_early_16k_rir(1:mix_16k_length);
                    s2_early_16k_rir = s2_early_16k_rir(1:mix_16k_length);
                    tmp_length = min(length(s1_16k),length(s2_16k));
                    s1_16k = s1_16k(1:tmp_length);
                    s2_16k = s2_16k(1:tmp_length);
                    s1_16k = cat(1,s1_16k,zeros(mix_16k_length - tmp_length,1));
                    s2_16k = cat(1,s2_16k,zeros(mix_16k_length - tmp_length,1));
            end
            mix_16k = s1_16k_rir + s2_16k_rir + sum(noise_16k_rir,2);
            
            max_amp_16k = max(cat(1,abs(mix_16k(:)),abs(s1_16k_rir(:)),abs(s2_16k_rir(:)),abs(s1_early_16k_rir(:)),abs(s2_early_16k_rir(:)),abs(s1_16k(:)),abs(s2_16k(:))));
            mix_scaling_16k = 1/max_amp_16k*0.9;
            s1_16k_rir = mix_scaling_16k * s1_16k_rir;
            s2_16k_rir = mix_scaling_16k * s2_16k_rir;
            s1_early_16k_rir = mix_scaling_16k * s1_early_16k_rir;
            s2_early_16k_rir = mix_scaling_16k * s2_early_16k_rir;
            mix_16k = mix_scaling_16k * mix_16k;      
            s1_16k = mix_scaling_16k * s1_16k;
            s2_16k = mix_scaling_16k * s2_16k;      
            
            % save 8 kHz and 16 kHz mixtures, as well as
            % necessary scaling factors
            
            scaling_16k(i,1) = mix_scaling_16k*scale1;
            scaling_16k(i,2) = mix_scaling_16k*scale2;
            scaling16bit_16k(i) = mix_scaling_16k;

 
            
            audiowrite([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s1/' mix_name '.wav'],s1_16k_rir,fs);            
            audiowrite([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s2/' mix_name '.wav'],s2_16k_rir,fs);            
            audiowrite([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/mix/' mix_name '.wav'],mix_16k,fs);
            audiowrite([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s1_early/' mix_name '.wav'],s1_early_16k_rir,fs);
            audiowrite([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s2_early/' mix_name '.wav'],s2_early_16k_rir,fs);
            audiowrite([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s1_dry/' mix_name '.wav'],s1_16k,fs);
            audiowrite([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s2_dry/' mix_name '.wav'],s2_16k,fs);
      
            
        end        
        save([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/scaling.mat'],'scaling_16k','scaling16bit_16k');
        
        fclose(fid);
        fclose(fid_s1);
        fclose(fid_s2);
        fclose(fid_m);
        toc
    end
end
