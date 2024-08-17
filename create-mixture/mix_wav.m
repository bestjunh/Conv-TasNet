function [scaling_8k, scaling16bit_8k, scaling_16k,scaling16bit_16k] = mix_wav(i,var)
C = var.C;
fid_s1 = var.fid_s1;
fid_s2 = var.fid_s2;
fid_m = var.fid_m;
wsj0root = var.wsj0root;
output_dir8k = var.output_dir8k;
output_dir16k = var.output_dir16k;
useaudioread = var.useaudioread;
fs8k = var.fs8k;


[inwav1_dir,invwav1_name,inwav1_ext] = fileparts(C{1}{i});
[inwav2_dir,invwav2_name,inwav2_ext] = fileparts(C{3}{i});
fprintf(fid_s1,'%s\n',C{1}{i});
fprintf(fid_s2,'%s\n',C{3}{i});
inwav1_snr = C{2}(i);
inwav2_snr = C{4}(i);
mix_name = [invwav1_name,'_',num2str(inwav1_snr),'_',invwav2_name,'_',num2str(inwav2_snr)];
fprintf(fid_m,'%s\n',mix_name);

% get input wavs
if useaudioread
    [s1, fs] = audioread([wsj0root C{1}{i}]);
    s2       = audioread([wsj0root C{3}{i}]);
else
    [s1, fs] = wavread([wsj0root C{1}{i}]); %#ok<*DWVRD>
    s2       = wavread([wsj0root C{3}{i}]);
end

% resample, normalize 8 kHz file, save scaling factor
s1_8k=resample(s1,fs8k,fs);
[s1_8k,lev1]=activlev(s1_8k,fs8k,'n'); % y_norm = y /sqrt(lev);
s2_8k=resample(s2,fs8k,fs);
[s2_8k,lev2]=activlev(s2_8k,fs8k,'n');

weight_1=10^(inwav1_snr/20);
weight_2=10^(inwav2_snr/20);

s1_8k = weight_1 * s1_8k;
s2_8k = weight_2 * s2_8k;

switch var.min_max
    case 'max'
        mix_8k_length = max(length(s1_8k),length(s2_8k));
        s1_8k = cat(1,s1_8k,zeros(mix_8k_length - length(s1_8k),1));
        s2_8k = cat(1,s2_8k,zeros(mix_8k_length - length(s2_8k),1));
    case 'min'
        mix_8k_length = min(length(s1_8k),length(s2_8k));
        s1_8k = s1_8k(1:mix_8k_length);
        s2_8k = s2_8k(1:mix_8k_length);
end
mix_8k = s1_8k + s2_8k;

max_amp_8k = max(cat(1,abs(mix_8k(:)),abs(s1_8k(:)),abs(s2_8k(:))));
mix_scaling_8k = 1/max_amp_8k*0.9;
s1_8k = mix_scaling_8k * s1_8k;
s2_8k = mix_scaling_8k * s2_8k;
mix_8k = mix_scaling_8k * mix_8k;

% apply same gain to 16 kHz file
s1_16k = weight_1 * s1 / sqrt(lev1);
s2_16k = weight_2 * s2 / sqrt(lev2);

switch var.min_max
    case 'max'
        mix_16k_length = max(length(s1_16k),length(s2_16k));
        s1_16k = cat(1,s1_16k,zeros(mix_16k_length - length(s1_16k),1));
        s2_16k = cat(1,s2_16k,zeros(mix_16k_length - length(s2_16k),1));
    case 'min'
        mix_16k_length = min(length(s1_16k),length(s2_16k));
        s1_16k = s1_16k(1:mix_16k_length);
        s2_16k = s2_16k(1:mix_16k_length);
end
mix_16k = s1_16k + s2_16k;

max_amp_16k = max(cat(1,abs(mix_16k(:)),abs(s1_16k(:)),abs(s2_16k(:))));
mix_scaling_16k = 1/max_amp_16k*0.9;
s1_16k = mix_scaling_16k * s1_16k;
s2_16k = mix_scaling_16k * s2_16k;
mix_16k = mix_scaling_16k * mix_16k;

% save 8 kHz and 16 kHz mixtures, as well as
% necessary scaling factors
scaling_16k = zeros(1,2);
scaling_8k = zeros(1,2);

scaling_16k(1) = weight_1 * mix_scaling_16k/ sqrt(lev1);
scaling_16k(2) = weight_2 * mix_scaling_16k/ sqrt(lev2);
scaling_8k(1) = weight_1 * mix_scaling_8k/ sqrt(lev1);
scaling_8k(2) = weight_2 * mix_scaling_8k/ sqrt(lev2);

scaling16bit_16k = mix_scaling_16k;
scaling16bit_8k  = mix_scaling_8k;

if useaudioread
    s1_8k = int16(round((2^15)*s1_8k));
    s2_8k = int16(round((2^15)*s2_8k));
    mix_8k = int16(round((2^15)*mix_8k));
    s1_16k = int16(round((2^15)*s1_16k));
    s2_16k = int16(round((2^15)*s2_16k));
    mix_16k = int16(round((2^15)*mix_16k));
    audiowrite([output_dir8k '/' var.min_max '/' var.data_type '/s1/' mix_name '.wav'],s1_8k,fs8k);
    audiowrite([output_dir16k '/' var.min_max '/' var.data_type '/s1/' mix_name '.wav'],s1_16k,fs);
    audiowrite([output_dir8k '/' var.min_max '/' var.data_type '/s2/' mix_name '.wav'],s2_8k,fs8k);
    audiowrite([output_dir16k '/' var.min_max '/' var.data_type '/s2/' mix_name '.wav'],s2_16k,fs);
    audiowrite([output_dir8k '/' var.min_max '/' var.data_type '/mix/' mix_name '.wav'],mix_8k,fs8k);
    audiowrite([output_dir16k '/' var.min_max '/' var.data_type '/mix/' mix_name '.wav'],mix_16k,fs);
else
    wavwrite(s1_8k,fs8k,[output_dir8k '/' var.min_max '/' var.data_type '/s1/' mix_name '.wav']); %#ok<*DWVWR>
    wavwrite(s1_16k,fs,[output_dir16k '/' var.min_max '/' var.data_type '/s1/' mix_name '.wav']);
    wavwrite(s2_8k,fs8k,[output_dir8k '/' var.min_max '/' var.data_type '/s2/' mix_name '.wav']);
    wavwrite(s2_16k,fs,[output_dir16k '/' var.min_max '/' var.data_type '/s2/' mix_name '.wav']);
    wavwrite(mix_8k,fs8k,[output_dir8k '/' var.min_max '/' var.data_type '/mix/' mix_name '.wav']);
    wavwrite(mix_16k,fs,[output_dir16k '/' var.min_max '/' var.data_type '/mix/' mix_name '.wav']);
end

% if mod(i,10)==0
%     fprintf(1,'.');
%     if mod(i,200)==0
%         fprintf(1,'\n');
%     end
% end
end