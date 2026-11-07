
%% RID-Based Phase-Locking Value (PLV) Computation
% Gamma-band EEG functional connectivity (30-50 Hz)
% Society for Neuroscience 2026
%
% Computes participant-averaged PLV matrices from preprocessed EEG
% using RID-Rihaczek time-frequency representations.
%
% Conditions:
%   PRE resting state, PRE meditation
%   POST resting state, POST meditation
%
% Parameters:
%   Gamma band:       30-50 Hz
%   RID frequency bins: 256
%   Window length:    400 ms
%   Window step:      200 ms (50% overlap)
%   Meditation period: minutes 2-6
%   Electrode subset: 12 frontoparietal channels
%
% Outputs:
%   Participant-averaged 12 x 12 PLV matrices
%   Per-participant checkpoint files
%
% Dependencies:
%   EEGLAB
%   rid_rihaczek4.m
%   MATLAB Parallel Computing Toolbox

clear; clc; close all;
restoredefaultpath;
rehash toolboxcache;

%% 1. Paths and dependencies

projectDir = '/mnt/home/ozbekogl/PLV_project';
rootDir = fullfile(projectDir,'data');

eeglabDir = fullfile(projectDir,'sccn-eeglab-6d3c273');
bctDir = fullfile(projectDir,'BCT');

outputDir = fullfile(projectDir,'output', ...
    'gamma_30to50_256bins_12ch_checkpoints');

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

addpath(projectDir);
addpath(fullfile(projectDir,'code'));

assert(exist(eeglabDir,'dir')==7, ...
    'EEGLAB directory not found.');

addpath(genpath(eeglabDir));

if exist(bctDir,'dir')
    addpath(genpath(bctDir));
end

assert(exist('pop_loadset','file')==2, ...
    'pop_loadset not found.');

assert(exist('rid_rihaczek4','file')==2, ...
    'rid_rihaczek4 not found.');

%% 2. Analysis parameters

excludeSubjects = ...
    [31 32 35 37 42 44 45 50 52 54 55 56 57 58 59];

sessions = {'ses-premedita','ses-posmedita'};

restTask = 'restCE01';
medTask = 'slMedita';

subNetLabels = { ...
    'F3','F4','Fz', ...
    'FC3','FC4','FCz', ...
    'P3','P4','Pz', ...
    'CP3','CP4','CPz'};

n_freqs = 256;

gamma_low = 30;
gamma_high = 50;

win_ms = 400;
step_ms = 200;

med_start_min = 2;
med_end_min = 6;

saveFile = fullfile(projectDir, ...
    'RID_PLV_frontoparietal_12ch_gamma_30to50_PRE_POST_2to6min_256BINS.mat');

fprintf('\n============================================\n');
fprintf('CORRECTED GAMMA RID-PLV EXTRACTION\n');
fprintf('============================================\n');

%% 3. Identify eligible participants

allSubs = dir(fullfile(rootDir,'sub-*'));
allSubs = allSubs([allSubs.isdir]);

subNames = {};
subIDs = [];

for k = 1:numel(allSubs)

    tok = regexp(allSubs(k).name, ...
        '^sub-(\d+)$','tokens','once');

    if isempty(tok)
        continue;
    end

    sid = str2double(tok{1});

    if ~ismember(sid,excludeSubjects)

        subNames{end+1} = allSubs(k).name;
        subIDs(end+1) = sid;

    end
end

[subIDs,sortIdx] = sort(subIDs);
subNames = subNames(sortIdx);

%% 4. Verify required EEG files

validSubs = {};
validIDs = [];

for s = 1:numel(subNames)

    subName = subNames{s};
    valid = true;

    for ss = 1:numel(sessions)

        sess = sessions{ss};

        for task = {restTask,medTask}

            taskName = task{1};

            filename = fullfile(rootDir,subName,sess, ...
                sprintf('%s_%s_task-%s_eeg_preproc_clean.set', ...
                subName,sess,taskName));

            if ~exist(filename,'file')

                fprintf('Missing: %s\n',filename);
                valid = false;

            end
        end
    end

    if valid

        validSubs{end+1} = subName;
        validIDs(end+1) = subIDs(s);

    end
end

nSubjects = numel(validSubs);

assert(nSubjects>0,'No valid subjects found.');

fprintf('Valid participants: %d\n',nSubjects);
disp(validSubs(:));

%% 5. Sampling rate and electrode selection

firstSub = validSubs{1};

firstFile = fullfile(rootDir,firstSub, ...
    'ses-premedita', ...
    sprintf('%s_ses-premedita_task-%s_eeg_preproc_clean.set', ...
    firstSub,restTask));

EEG_tmp = pop_loadset(firstFile);

fs = EEG_tmp.srate;

assert(abs(fs-250)<1e-6, ...
    'Expected 250 Hz sampling frequency.');

allLabels = {EEG_tmp.chanlocs.labels};

subIdx = zeros(1,numel(subNetLabels));

for c = 1:numel(subNetLabels)

    matches = find(strcmpi(strtrim(allLabels), ...
        subNetLabels{c}));

    assert(numel(matches)==1, ...
        'Missing or duplicate electrode: %s', ...
        subNetLabels{c});

    subIdx(c) = matches;

end

nChan = numel(subIdx);
chanlocs = EEG_tmp.chanlocs(subIdx);

win_samples = round(win_ms*fs/1000);
step_samples = round(step_ms*fs/1000);

med_start_sample = round(med_start_min*60*fs)+1;
med_end_sample = round(med_end_min*60*fs);

%% 6. RID frequency-axis correction

% RID uses unshifted FFT frequency indexing.
% Frequency spacing is fs/n_freqs, not fs/(2*n_freqs).
% This ensures that the selected bins correspond to 30-50 Hz.

freq_axis = (0:n_freqs-1)*(fs/n_freqs);

gamma_bins = find( ...
    freq_axis>=gamma_low & ...
    freq_axis<=gamma_high);

assert(~isempty(gamma_bins), ...
    'No gamma frequency bins selected.');

fprintf('\nSampling frequency: %.2f Hz\n',fs);
fprintf('Frequency bins: %d\n',n_freqs);
fprintf('Frequency resolution: %.6f Hz\n',fs/n_freqs);

fprintf('Gamma bin indices:\n');
disp(gamma_bins);

fprintf('Actual gamma frequencies:\n');
disp(freq_axis(gamma_bins));

fprintf('Gamma range: %.4f - %.4f Hz\n', ...
    freq_axis(gamma_bins(1)), ...
    freq_axis(gamma_bins(end)));

fprintf('\nSelected electrodes:\n');
disp(subNetLabels(:));

fprintf('Number of unique edges: %d\n', ...
    nChan*(nChan-1)/2);

fprintf('Window: %d ms\n',win_ms);
fprintf('Step: %d ms\n',step_ms);

%% 7. Initialize output storage

PLV_pre_rest = cell(1,nSubjects);
PLV_pre_med = cell(1,nSubjects);

PLV_post_rest = cell(1,nSubjects);
PLV_post_med = cell(1,nSubjects);

t_pre_rest = cell(1,nSubjects);
t_pre_med = cell(1,nSubjects);

t_post_rest = cell(1,nSubjects);
t_post_med = cell(1,nSubjects);

nWin_pre_rest = nan(nSubjects,1);
nWin_pre_med = nan(nSubjects,1);

nWin_post_rest = nan(nSubjects,1);
nWin_post_med = nan(nSubjects,1);

valid_pre_pair = false(nSubjects,1);
valid_post_pair = false(nSubjects,1);

%% 8. Initialize parallel processing

nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));

if isnan(nWorkers) || nWorkers<1
    nWorkers = 4;
end

pool = gcp('nocreate');

if isempty(pool)
    pool = parpool('local',nWorkers);
end

pctRunOnAll addpath(projectDir);
pctRunOnAll addpath(fullfile(projectDir,'code'));
pctRunOnAll addpath(genpath(eeglabDir));

fprintf('\nStarting parallel gamma extraction...\n');

%% 9. Compute gamma-band RID-PLV for each participant

parfor s = 1:nSubjects

    subName = validSubs{s};

    checkpointFile = fullfile(outputDir, ...
        sprintf('%s_gamma_256bins_12ch_checkpoint.mat', ...
        subName));

    if exist(checkpointFile,'file')

        C = load(checkpointFile);

        % Reuse only checkpoints matching the analysis configuration.
        assert(C.n_freqs==n_freqs && ...
            isequal(C.gamma_bins,gamma_bins) && ...
            isequal(C.subNetLabels,subNetLabels), ...
            'Incompatible checkpoint for %s.',subName);

        PLV_pre_rest{s} = C.local_PLV_pre_rest;
        PLV_pre_med{s} = C.local_PLV_pre_med;

        PLV_post_rest{s} = C.local_PLV_post_rest;
        PLV_post_med{s} = C.local_PLV_post_med;

        t_pre_rest{s} = C.local_t_pre_rest;
        t_pre_med{s} = C.local_t_pre_med;

        t_post_rest{s} = C.local_t_post_rest;
        t_post_med{s} = C.local_t_post_med;

        nWin_pre_rest(s) = C.local_nWin_pre_rest;
        nWin_pre_med(s) = C.local_nWin_pre_med;

        nWin_post_rest(s) = C.local_nWin_post_rest;
        nWin_post_med(s) = C.local_nWin_post_med;

        valid_pre_pair(s) = C.local_valid_pre_pair;
        valid_post_pair(s) = C.local_valid_post_pair;

        fprintf('Loaded checkpoint: %s\n',subName);

        continue;
    end

    fprintf('Processing %s\n',subName);

    local_PLV_pre_rest = [];
    local_PLV_pre_med = [];

    local_PLV_post_rest = [];
    local_PLV_post_med = [];

    local_t_pre_rest = [];
    local_t_pre_med = [];

    local_t_post_rest = [];
    local_t_post_med = [];

    local_nWin_pre_rest = NaN;
    local_nWin_pre_med = NaN;

    local_nWin_post_rest = NaN;
    local_nWin_post_med = NaN;

    local_valid_pre_pair = false;
    local_valid_post_pair = false;

    %% PRE session

    preSess = 'ses-premedita';

    preRestFile = fullfile(rootDir,subName,preSess, ...
        sprintf('%s_%s_task-%s_eeg_preproc_clean.set', ...
        subName,preSess,restTask));

    preMedFile = fullfile(rootDir,subName,preSess, ...
        sprintf('%s_%s_task-%s_eeg_preproc_clean.set', ...
        subName,preSess,medTask));

    EEG_pre_rest = pop_loadset(preRestFile);
    EEG_pre_med = pop_loadset(preMedFile);

    assert(EEG_pre_rest.srate==fs && ...
        EEG_pre_med.srate==fs, ...
        'PRE sampling rate mismatch.');

    data_pre_rest = select_subnetwork( ...
        EEG_pre_rest,subNetLabels);

    data_pre_med_full = select_subnetwork( ...
        EEG_pre_med,subNetLabels);

    if size(data_pre_med_full,2)>=med_end_sample

        data_pre_med = data_pre_med_full( ...
            :,med_start_sample:med_end_sample);

        [local_PLV_pre_rest, ...
         local_t_pre_rest, ...
         local_nWin_pre_rest] = compute_rid_plv( ...
            data_pre_rest,fs,win_samples, ...
            step_samples,n_freqs,gamma_bins);

        [local_PLV_pre_med, ...
         local_t_pre_med, ...
         local_nWin_pre_med] = compute_rid_plv( ...
            data_pre_med,fs,win_samples, ...
            step_samples,n_freqs,gamma_bins);

        local_valid_pre_pair = true;

    end

    %% POST session

    postSess = 'ses-posmedita';

    postRestFile = fullfile(rootDir,subName,postSess, ...
        sprintf('%s_%s_task-%s_eeg_preproc_clean.set', ...
        subName,postSess,restTask));

    postMedFile = fullfile(rootDir,subName,postSess, ...
        sprintf('%s_%s_task-%s_eeg_preproc_clean.set', ...
        subName,postSess,medTask));

    EEG_post_rest = pop_loadset(postRestFile);
    EEG_post_med = pop_loadset(postMedFile);

    assert(EEG_post_rest.srate==fs && ...
        EEG_post_med.srate==fs, ...
        'POST sampling rate mismatch.');

    data_post_rest = select_subnetwork( ...
        EEG_post_rest,subNetLabels);

    data_post_med_full = select_subnetwork( ...
        EEG_post_med,subNetLabels);

    if size(data_post_med_full,2)>=med_end_sample

        data_post_med = data_post_med_full( ...
            :,med_start_sample:med_end_sample);

        [local_PLV_post_rest, ...
         local_t_post_rest, ...
         local_nWin_post_rest] = compute_rid_plv( ...
            data_post_rest,fs,win_samples, ...
            step_samples,n_freqs,gamma_bins);

        [local_PLV_post_med, ...
         local_t_post_med, ...
         local_nWin_post_med] = compute_rid_plv( ...
            data_post_med,fs,win_samples, ...
            step_samples,n_freqs,gamma_bins);

        local_valid_post_pair = true;

    end

    %% Assign participant results

    PLV_pre_rest{s} = local_PLV_pre_rest;
    PLV_pre_med{s} = local_PLV_pre_med;

    PLV_post_rest{s} = local_PLV_post_rest;
    PLV_post_med{s} = local_PLV_post_med;

    t_pre_rest{s} = local_t_pre_rest;
    t_pre_med{s} = local_t_pre_med;

    t_post_rest{s} = local_t_post_rest;
    t_post_med{s} = local_t_post_med;

    nWin_pre_rest(s) = local_nWin_pre_rest;
    nWin_pre_med(s) = local_nWin_pre_med;

    nWin_post_rest(s) = local_nWin_post_rest;
    nWin_post_med(s) = local_nWin_post_med;

    valid_pre_pair(s) = local_valid_pre_pair;
    valid_post_pair(s) = local_valid_post_pair;

    %% Save participant checkpoint

    checkpoint = struct();

    checkpoint.subName = subName;
    checkpoint.subID = validIDs(s);

    checkpoint.local_PLV_pre_rest = local_PLV_pre_rest;
    checkpoint.local_PLV_pre_med = local_PLV_pre_med;

    checkpoint.local_PLV_post_rest = local_PLV_post_rest;
    checkpoint.local_PLV_post_med = local_PLV_post_med;

    checkpoint.local_t_pre_rest = local_t_pre_rest;
    checkpoint.local_t_pre_med = local_t_pre_med;

    checkpoint.local_t_post_rest = local_t_post_rest;
    checkpoint.local_t_post_med = local_t_post_med;

    checkpoint.local_nWin_pre_rest = local_nWin_pre_rest;
    checkpoint.local_nWin_pre_med = local_nWin_pre_med;

    checkpoint.local_nWin_post_rest = local_nWin_post_rest;
    checkpoint.local_nWin_post_med = local_nWin_post_med;

    checkpoint.local_valid_pre_pair = local_valid_pre_pair;
    checkpoint.local_valid_post_pair = local_valid_post_pair;

    checkpoint.n_freqs = n_freqs;
    checkpoint.gamma_bins = gamma_bins;
    checkpoint.freq_axis = freq_axis;

    checkpoint.subNetLabels = subNetLabels;

    checkpoint.fs = fs;
    checkpoint.win_ms = win_ms;
    checkpoint.step_ms = step_ms;

    checkpoint.med_start_min = med_start_min;
    checkpoint.med_end_min = med_end_min;

    save_checkpoint(checkpointFile,checkpoint);

    fprintf('Finished %s\n',subName);

end

%% 10. Processing summary

fprintf('\n============================================\n');
fprintf('GAMMA EXTRACTION SUMMARY\n');
fprintf('============================================\n');

for s = 1:nSubjects

    fprintf('%s | PRE=%d | POST=%d\n', ...
        validSubs{s}, ...
        valid_pre_pair(s), ...
        valid_post_pair(s));

end

fprintf('Valid PRE: %d/%d\n', ...
    sum(valid_pre_pair),nSubjects);

fprintf('Valid POST: %d/%d\n', ...
    sum(valid_post_pair),nSubjects);

%% 11. Save combined participant-level PLV matrices

save(saveFile, ...
    'validSubs','validIDs','excludeSubjects', ...
    'sessions','restTask','medTask', ...
    'subNetLabels','subIdx', ...
    'chanlocs','nChan', ...
    'fs','n_freqs','freq_axis','gamma_bins', ...
    'win_ms','step_ms', ...
    'win_samples','step_samples', ...
    'med_start_min','med_end_min', ...
    'med_start_sample','med_end_sample', ...
    'PLV_pre_rest','PLV_pre_med', ...
    'PLV_post_rest','PLV_post_med', ...
    't_pre_rest','t_pre_med', ...
    't_post_rest','t_post_med', ...
    'nWin_pre_rest','nWin_pre_med', ...
    'nWin_post_rest','nWin_post_med', ...
    'valid_pre_pair','valid_post_pair', ...
    'outputDir','-v7.3');

fprintf('\nSaved combined output:\n%s\n',saveFile);
fprintf('\nGAMMA EXTRACTION COMPLETE\n');

%% Local function: RID-PLV computation

function [PLV_avg,win_centers_sec,nWin] = ...
    compute_rid_plv(data,fs,win_samples, ...
    step_samples,n_freqs,gamma_bins)

    nChan = size(data,1);
    nSamples = size(data,2);

    start_idx = 1:step_samples: ...
        (nSamples-win_samples+1);

    nWin = numel(start_idx);

    win_centers_sec = nan(nWin,1);

    if nWin<1
        PLV_avg = nan(nChan,nChan);
        return;
    end

    PLV_sum = zeros(nChan,nChan);
    PLV_count = zeros(nChan,nChan);

    for w = 1:nWin

        idx = start_idx(w): ...
            (start_idx(w)+win_samples-1);

        xw = data(:,idx);

        W = nan(nChan,nChan);

        for i = 1:nChan
            for j = i+1:nChan

                data1 = xw(i,:);
                data2 = xw(j,:);

                % RID-Rihaczek time-frequency representations
                Rx = rid_rihaczek4(data1,n_freqs);
                Ry = rid_rihaczek4(data2,n_freqs);

                % Normalize to obtain phase
                Rx_phase = Rx./(abs(Rx)+eps);
                Ry_phase = Ry./(abs(Ry)+eps);

                % Compute phase difference
                phdif = Rx_phase.*conj(Ry_phase);

                % PLV averaged across time for each frequency
                plv_f = abs(mean(phdif,2));

                % Average PLV across the selected gamma frequencies
                gamma_plv = mean( ...
                    plv_f(gamma_bins),'omitnan');

                W(i,j) = gamma_plv;
                W(j,i) = gamma_plv;

            end
        end

        W(1:nChan+1:end) = NaN;

        validMask = isfinite(W);

        PLV_sum(validMask) = ...
            PLV_sum(validMask)+W(validMask);

        PLV_count(validMask) = ...
            PLV_count(validMask)+1;

        win_centers_sec(w) = mean(idx)/fs;

    end

    PLV_avg = PLV_sum./PLV_count;

    PLV_avg(PLV_count==0) = NaN;
    PLV_avg(1:nChan+1:end) = NaN;

end

%% Local function: select frontoparietal electrodes

function data = select_subnetwork(EEG,subNetLabels)

    labels = {EEG.chanlocs.labels};

    idx = zeros(1,numel(subNetLabels));

    for c = 1:numel(subNetLabels)

        matches = find(strcmpi( ...
            strtrim(labels),subNetLabels{c}));

        assert(numel(matches)==1, ...
            'Missing or duplicate electrode: %s', ...
            subNetLabels{c});

        idx(c) = matches;

    end

    data = double(EEG.data(idx,:));

end

%% Local function: save participant checkpoint

function save_checkpoint(filename,checkpoint)

    save(filename,'-struct','checkpoint','-v7.3');

end
