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

if ispc
    nasPath = 'Z:/nas1_data/';
    devPath = 'W:/';
elseif isunix
    nasPath = '/home/nas/';
    devPath = '/home/';
else
    disp('Unknown operating system.');
end


data_type = {'tr','cv','tt'};
wsj0root = [nasPath 'DB/wsj_wav/']; % YOUR_PATH/, the folder containing wsj0/
output_dir16k=[devPath 'data/albert/DB/wsj0-mix_20240717/2speakers/wav16k'];
output_dir8k=[devPath 'data/albert/DB/wsj0-mix_20240717/2speakers/wav8k'];

min_max = {'min','max'};

useaudioread = 0;
if exist('audioread','file')
    useaudioread = 1;
end

% 기존 병렬 풀 종료
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    delete(poolobj);
end

% 새로운 병렬 풀을 워커 수와 함께 시작
numWorkers = 32;
parpool('local', numWorkers);

disp(['Parallel pool restarted with ', num2str(numWorkers), ' workers.']);


for i_mm = 1:length(min_max)
    for i_type = 1:length(data_type)
        if ~exist([output_dir16k '/' min_max{i_mm} '/' data_type{i_type}],'dir')
            mkdir([output_dir16k '/' min_max{i_mm} '/' data_type{i_type}]);
        end
        if ~exist([output_dir8k '/' min_max{i_mm} '/' data_type{i_type}],'dir')
            mkdir([output_dir8k '/' min_max{i_mm} '/' data_type{i_type}]);
        end
        status = mkdir([output_dir8k  '/' min_max{i_mm} '/' data_type{i_type} '/s1/']); %#ok<NASGU>
        status = mkdir([output_dir8k  '/' min_max{i_mm} '/' data_type{i_type} '/s2/']); %#ok<NASGU>
        status = mkdir([output_dir8k  '/' min_max{i_mm} '/' data_type{i_type} '/mix/']); %#ok<NASGU>
        status = mkdir([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s1/']); %#ok<NASGU>
        status = mkdir([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/s2/']); %#ok<NASGU>
        status = mkdir([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/mix/']);
                
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
        fs8k=8000;
        
        scaling_16k = zeros(num_files,2);
        scaling_8k = zeros(num_files,2);
        scaling16bit_16k = zeros(num_files,1);
        scaling16bit_8k = zeros(num_files,1);
        fprintf(1,'%s\n',[min_max{i_mm} '_' data_type{i_type}]);

        var.C = C;
        var.fid_s1 = fid_s1;
        var.fid_s2 = fid_s2;
        var.fid_m = fid_m;
        var.wsj0root = wsj0root;
        var.output_dir8k = output_dir8k;
        var.output_dir16k = output_dir16k;
        var.min_max = min_max{i_mm};
        var.data_type = data_type{i_type};
        var.useaudioread = useaudioread;
        var.fs8k = fs8k;

        parfor i = 1:num_files
            disp(['Creating ' num2str(i_mm) '/' num2str(length(min_max)) ' ' num2str(i_type) '/' num2str(length(data_type)) ' ' num2str(i) '/' num2str(num_files)]);
            [scaling_8k(i,:), scaling16bit_8k(i), scaling_16k(i,:), scaling16bit_16k(i)] = mix_wav(i,var);
            
        end
        save([output_dir8k  '/' min_max{i_mm} '/' data_type{i_type} '/scaling.mat'],'scaling_8k','scaling16bit_8k');
        save([output_dir16k '/' min_max{i_mm} '/' data_type{i_type} '/scaling.mat'],'scaling_16k','scaling16bit_16k');
        
        fclose(fid);
        fclose(fid_s1);
        fclose(fid_s2);
        fclose(fid_m);
    end
end
