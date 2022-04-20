% Main script to run card_interpolate.mat
% Instructions: 
% 1. Replace ['...'] with the directory in which the waveform files are
% stored.
% 2. If the files are compressed in .gz format, uncomment the according
% section.
% 3. The preprocessed files and fixed files will be in two separate folders
% under the inputted directory, with '_physOUT' and '_physOUT_fixed'
% suffix.

% Input directory of PPG waveforms here
user_dir = ['...'];
file_list = dir([user_dir, '/*.tsv']);

if ~exist([user_dir, '/fixed'],'dir')
    mkdir([user_dir, '/fixed']);
end

% Pure text format (csv, tsv, etc.) required
for i = 1:length(file_list)
    file_name = [file_list(i).folder, '/', file_list(i).name];
    
    % Decompress a .gz file if necessary
    % extrt_file = gunzip(file_name);
    % file_name = extrt_file{1,1};
    
    % Preprocess the physio waveforms to prepare for the fix
    preproc_physio(file_name, [file_list(i).folder, '/preprocessed'], 500, 700);
    close all;
    clean_name = extractBefore(file_list(i).name, '.');
    mat_name = [clean_name, '_physOUT.mat'];

    % Interpolate the waveform and display difference in figure
    card_interpolate([file_list(i).folder, '/preprocessed/', mat_name]);
    card_view(file_list(i).folder, mat_name);
end

function card_view(file_folder, file_name)
% Helper Function: card_view
%
% Given file_folder and file_name, display in new figure an overall
% contrast between the original(preprocessed) and the fixed waveform.
% Also marks the corrected peak locations on top of the corrected waveform.

% Load the preprocessed and the fixed files
[~,fn,ext] = fileparts(file_name);
old_file = load([file_folder, '/preprocessed/', file_name]);
new_file = load([file_folder, '/fixed/', fn, '_fixed', ext]);

% Plot the two waveforms and the peak locations
figure, plot(old_file.OUT_p.card_bpf, 'r');
hold on, plot(new_file.OUT_p.card_bpf, 'g');
plot(new_file.OUT_p.card_trig_samples, ...
     new_file.OUT_p.card_bpf(new_file.OUT_p.card_trig_samples), '*b');
title(fn);

% Continue to next file
input('Press Enter to Continue');
close all;
end