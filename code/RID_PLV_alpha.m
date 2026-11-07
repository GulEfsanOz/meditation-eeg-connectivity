
%% Split-alpha EEG connectivity | RID–PLV extraction
% Supplementary analysis code — Society for Neuroscience 2026
% Fronto-parietal 12-channel subnetwork; 256 RID frequency bins.
% Requires EEGLAB (pop_loadset), rid_rihaczek4, and MATLAB Parallel Computing Toolbox.
% The input EEG data and external dependencies are not included in this file.
%
% Conditions:
%   PRE rest, PRE meditation, POST rest, POST meditation
%
% Frequency bands:
%   Low alpha:  8-10 Hz
%   High alpha: 10-12 Hz
%
% Outputs:
%   Subject average: 12 x 12 x 2
%   Window level:    12 x 12 x 2 x nWindows
%% ================================================================

clear;
clc;
close all;

restoredefaultpath;
rehash toolboxcache;

fprintf('\n===== RID-PLV SPLIT ALPHA: 256 BINS =====\n');

%% ================================================================
% PATHS
%% ================================================================

projectDir = fullfile(getenv('HOME'), 'PLV_project');
rootDir = fullfile(projectDir, 'data');
outputDir = fullfile(projectDir, 'data');

eeglabDir = fullfile(projectDir, 'sccn-eeglab-6d3c273');
bctDir = fullfile(projectDir, 'BCT');

addpath(projectDir);
addpath(genpath(eeglabDir));

if exist(bctDir, 'dir')
    addpath(genpath(bctDir));
end

assert(exist('pop_loadset', 'file') == 2, ...
    'pop_loadset not found.');

assert(exist('rid_rihaczek4', 'file') == 2, ...
    'rid_rihaczek4 not found.');

%% ================================================================
% PARAMETERS
%% ================================================================

excludeSubjects = ...
    [31 32 35 37 42 44 45 50 52 54 55 56 57 58 59];

sessions = {'ses-premedita', 'ses-posmedita'};

restTask = 'restCE01';
medTask = 'slMedita';

% Number of RID frequency bins used in the analysis
n_freqs = 256;

% 5-second windows, 50% overlap
win_ms = 5000;
step_ms = 2500;

% Meditation interval
med_start_min = 2;
med_end_min = 6;

% Split alpha bands
bandNames = {'alpha_8_10', 'alpha_10_12'};

bandEdges = [8 10;
             10 12];

% 12-electrode subnetwork
subNetLabels = { ...
    'F3','F4','Fz', ...
    'FC3','FC4','FCz', ...
    'P3','P4','Pz', ...
    'CP3','CP4','CPz'};

% Separate output filename to protect previous results
outputFile = fullfile(outputDir, ...
    ['RID_PLV_frontoparietal_12ch_alpha_' ...
     '8to10_10to12_PRE_POST_2to6min_' ...
     'WITH_WINDOWS_256BINS.mat']);

%% ================================================================
% BUILD SUBJECT LIST
%% ================================================================

allSubs = dir(fullfile(rootDir, 'sub-*'));
allSubs = allSubs([allSubs.isdir]);

subNames = {};
subIDs = [];

for k = 1:numel(allSubs)

    tok = regexp(allSubs(k).name, ...
        '^sub-(\d+)$', 'tokens', 'once');

    if isempty(tok)
        continue;
    end

    sid = str2double(tok{1});

    if ismember(sid, excludeSubjects)
        continue;
    end

    subNames{end+1} = allSubs(k).name;
    subIDs(end+1) = sid;

end

[subIDs, sortIdx] = sort(subIDs);
subNames = subNames(sortIdx);

%% ================================================================
% KEEP SUBJECTS WITH ALL FOUR REQUIRED EEG FILES
%% ================================================================

validSubs = {};
validIDs = [];

for s = 1:numel(subNames)

    subName = subNames{s};
    ok = true;

    for ss = 1:numel(sessions)

        sess = sessions{ss};

        restFile = fullfile(rootDir, subName, sess, ...
            sprintf('%s_%s_task-%s_eeg_preproc_clean.set', ...
            subName, sess, restTask));

        medFile = fullfile(rootDir, subName, sess, ...
            sprintf('%s_%s_task-%s_eeg_preproc_clean.set', ...
            subName, sess, medTask));

        if ~exist(restFile, 'file') || ...
           ~exist(medFile, 'file')

            ok = false;
            break;
        end

    end

    if ok
        validSubs{end+1} = subName;
        validIDs(end+1) = subIDs(s);
    end

end

nSubjects = numel(validSubs);

fprintf('Using %d valid subjects.\n', nSubjects);

assert(nSubjects > 0, ...
    'No valid subjects found. Check data paths.');

%% ================================================================
% GET SAMPLING RATE AND CHANNEL INFORMATION
%% ================================================================

firstSub = validSubs{1};

tmpFile = fullfile(rootDir, firstSub, ...
    'ses-premedita', ...
    sprintf('%s_ses-premedita_task-%s_eeg_preproc_clean.set', ...
    firstSub, restTask));

EEG_tmp = pop_loadset(tmpFile);

% Read actual sampling frequency from EEG
fs = EEG_tmp.srate;

chanlocs_full = EEG_tmp.chanlocs;
allLabels = string({chanlocs_full.labels});

win_samples = round(win_ms * fs / 1000);
step_samples = round(step_ms * fs / 1000);

med_start_sample = round(med_start_min * 60 * fs) + 1;
med_end_sample = round(med_end_min * 60 * fs);

%% ================================================================
% CORRECT RID FREQUENCY AXIS
%% ================================================================

% Match the unshifted FFT indexing used by rid_rihaczek4.
% The spacing is fs/n_freqs, not fs/(2*n_freqs).

freq_axis = (0:n_freqs-1) * fs / n_freqs;

bandBins = cell(1, numel(bandNames));

for b = 1:numel(bandNames)

    bandBins{b} = find( ...
        freq_axis >= bandEdges(b,1) & ...
        freq_axis <= bandEdges(b,2));

    assert(~isempty(bandBins{b}), ...
        'No frequency bins found for %s.', bandNames{b});

    fprintf('\n%s\n', bandNames{b});

    fprintf('Selected frequencies: %s Hz\n', ...
        mat2str(freq_axis(bandBins{b}), 8));

end

% At fs = 250 Hz and n_freqs = 256:
%
% 8-10 Hz:
%   8.7890625, 9.765625
%
% 10-12 Hz:
%   10.7421875, 11.71875

%% ================================================================
% IDENTIFY THE 12 ELECTRODES
%% ================================================================

subIdx = zeros(1, numel(subNetLabels));

for k = 1:numel(subNetLabels)

    idx = find(strcmpi(allLabels, subNetLabels{k}), 1);

    if isempty(idx)
        error('Could not find electrode %s.', ...
            subNetLabels{k});
    end

    subIdx(k) = idx;

end

subnetLabels = string(subNetLabels);
chanlocs = chanlocs_full(subIdx);
nNodes = numel(subIdx);

fprintf('\nSampling rate: %.2f Hz\n', fs);
fprintf('Frequency bins: %d\n', n_freqs);
fprintf('Window length: %d samples\n', win_samples);
fprintf('Step length: %d samples\n', step_samples);

fprintf('Meditation crop: %.1f-%.1f minutes\n', ...
    med_start_min, med_end_min);

disp(table(subnetLabels(:), subIdx(:), ...
    'VariableNames', {'Electrode','OriginalIndex'}));

%% ================================================================
% STORAGE
%% ================================================================

% Subject-level averages
PLV_pre_rest = cell(1, nSubjects);
PLV_pre_med = cell(1, nSubjects);

PLV_post_rest = cell(1, nSubjects);
PLV_post_med = cell(1, nSubjects);

% Individual-window PLVs
PLVwin_pre_rest = cell(1, nSubjects);
PLVwin_pre_med = cell(1, nSubjects);

PLVwin_post_rest = cell(1, nSubjects);
PLVwin_post_med = cell(1, nSubjects);

% Window counts
nWin_pre_rest = nan(nSubjects,1);
nWin_pre_med = nan(nSubjects,1);

nWin_post_rest = nan(nSubjects,1);
nWin_post_med = nan(nSubjects,1);

% Validity indicators
valid_pre_pair = false(nSubjects,1);
valid_post_pair = false(nSubjects,1);

%% ================================================================
% START PARALLEL POOL
%% ================================================================

pool = gcp('nocreate');

if isempty(pool)

    nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));

    if isnan(nWorkers) || nWorkers < 1
        nWorkers = 1;
    end

    fprintf('Starting %d parallel workers.\n', nWorkers);

    pool = parpool('local', nWorkers);

end

pctRunOnAll addpath(projectDir);
pctRunOnAll addpath(genpath(eeglabDir));

if exist(bctDir, 'dir')
    pctRunOnAll addpath(genpath(bctDir));
end

dq = parallel.pool.DataQueue;
afterEach(dq, @printProgress);

fprintf('\nStarting parallel PLV computation...\n');

%% ================================================================
% MAIN PARALLEL SUBJECT LOOP
%% ================================================================

parfor s = 1:nSubjects

    subName = validSubs{s};

    send(dq, sprintf('START subject %d/%d: %s', ...
        s, nSubjects, subName));

    % Local storage
    local_PLV_pre_rest = [];
    local_PLV_pre_med = [];
    local_PLV_post_rest = [];
    local_PLV_post_med = [];

    local_PLVwin_pre_rest = [];
    local_PLVwin_pre_med = [];
    local_PLVwin_post_rest = [];
    local_PLVwin_post_med = [];

    local_nWin_pre_rest = NaN;
    local_nWin_pre_med = NaN;
    local_nWin_post_rest = NaN;
    local_nWin_post_med = NaN;

    local_valid_pre_pair = false;
    local_valid_post_pair = false;

    %% ------------------------------------------------------------
    % PRE SESSION
    %% ------------------------------------------------------------

    % Build the session-specific EEG filenames.

    pre_rest_file = fullfile(rootDir, subName, ...
        'ses-premedita', ...
        sprintf('%s_ses-premedita_task-%s_eeg_preproc_clean.set', ...
        subName, restTask));

    pre_med_file = fullfile(rootDir, subName, ...
        'ses-premedita', ...
        sprintf('%s_ses-premedita_task-%s_eeg_preproc_clean.set', ...
        subName, medTask));

    EEG_pre_rest = pop_loadset(pre_rest_file);
    EEG_pre_med = pop_loadset(pre_med_file);

    data_pre_rest = double(EEG_pre_rest.data(subIdx,:));
    data_pre_med_full = double(EEG_pre_med.data(subIdx,:));

    pre_med_ok = size(data_pre_med_full,2) >= med_end_sample;

    local_valid_pre_pair = pre_med_ok;

    if pre_med_ok

        data_pre_med = data_pre_med_full( ...
            :, med_start_sample:med_end_sample);

        [local_PLV_pre_rest, ...
         local_PLVwin_pre_rest, ...
         local_nWin_pre_rest] = ...
            compute_split_alpha_subjectavg_plv( ...
                data_pre_rest, ...
                win_samples, step_samples, ...
                n_freqs, bandBins);

        [local_PLV_pre_med, ...
         local_PLVwin_pre_med, ...
         local_nWin_pre_med] = ...
            compute_split_alpha_subjectavg_plv( ...
                data_pre_med, ...
                win_samples, step_samples, ...
                n_freqs, bandBins);

    end

    send(dq, sprintf('PRE DONE subject %d/%d: %s', ...
        s, nSubjects, subName));

    %% ------------------------------------------------------------
    % POST SESSION
    %% ------------------------------------------------------------

    % Build the session-specific EEG filenames.

    post_rest_file = fullfile(rootDir, subName, ...
        'ses-posmedita', ...
        sprintf('%s_ses-posmedita_task-%s_eeg_preproc_clean.set', ...
        subName, restTask));

    post_med_file = fullfile(rootDir, subName, ...
        'ses-posmedita', ...
        sprintf('%s_ses-posmedita_task-%s_eeg_preproc_clean.set', ...
        subName, medTask));

    EEG_post_rest = pop_loadset(post_rest_file);
    EEG_post_med = pop_loadset(post_med_file);

    data_post_rest = double(EEG_post_rest.data(subIdx,:));
    data_post_med_full = double(EEG_post_med.data(subIdx,:));

    post_med_ok = size(data_post_med_full,2) >= med_end_sample;

    local_valid_post_pair = post_med_ok;

    if post_med_ok

        data_post_med = data_post_med_full( ...
            :, med_start_sample:med_end_sample);

        [local_PLV_post_rest, ...
         local_PLVwin_post_rest, ...
         local_nWin_post_rest] = ...
            compute_split_alpha_subjectavg_plv( ...
                data_post_rest, ...
                win_samples, step_samples, ...
                n_freqs, bandBins);

        [local_PLV_post_med, ...
         local_PLVwin_post_med, ...
         local_nWin_post_med] = ...
            compute_split_alpha_subjectavg_plv( ...
                data_post_med, ...
                win_samples, step_samples, ...
                n_freqs, bandBins);

    end

    send(dq, sprintf('POST DONE subject %d/%d: %s', ...
        s, nSubjects, subName));

    %% ------------------------------------------------------------
    % WRITE BACK TO SUBJECT CELLS
    %% ------------------------------------------------------------

    PLV_pre_rest{s} = local_PLV_pre_rest;
    PLV_pre_med{s} = local_PLV_pre_med;

    PLV_post_rest{s} = local_PLV_post_rest;
    PLV_post_med{s} = local_PLV_post_med;

    PLVwin_pre_rest{s} = local_PLVwin_pre_rest;
    PLVwin_pre_med{s} = local_PLVwin_pre_med;

    PLVwin_post_rest{s} = local_PLVwin_post_rest;
    PLVwin_post_med{s} = local_PLVwin_post_med;

    nWin_pre_rest(s) = local_nWin_pre_rest;
    nWin_pre_med(s) = local_nWin_pre_med;

    nWin_post_rest(s) = local_nWin_post_rest;
    nWin_post_med(s) = local_nWin_post_med;

    valid_pre_pair(s) = local_valid_pre_pair;
    valid_post_pair(s) = local_valid_post_pair;

    send(dq, sprintf('FINISHED subject %d/%d: %s', ...
        s, nSubjects, subName));

end

%% ================================================================
% SAVE RESULTS
%% ================================================================

save(outputFile, ...
    'PLV_pre_rest', ...
    'PLV_pre_med', ...
    'PLV_post_rest', ...
    'PLV_post_med', ...
    'PLVwin_pre_rest', ...
    'PLVwin_pre_med', ...
    'PLVwin_post_rest', ...
    'PLVwin_post_med', ...
    'validSubs', ...
    'validIDs', ...
    'excludeSubjects', ...
    'sessions', ...
    'restTask', ...
    'medTask', ...
    'med_start_min', ...
    'med_end_min', ...
    'med_start_sample', ...
    'med_end_sample', ...
    'fs', ...
    'win_ms', ...
    'step_ms', ...
    'win_samples', ...
    'step_samples', ...
    'n_freqs', ...
    'freq_axis', ...
    'bandEdges', ...
    'bandNames', ...
    'bandBins', ...
    'subNetLabels', ...
    'subnetLabels', ...
    'subIdx', ...
    'chanlocs', ...
    'nWin_pre_rest', ...
    'nWin_pre_med', ...
    'nWin_post_rest', ...
    'nWin_post_med', ...
    'valid_pre_pair', ...
    'valid_post_pair', ...
    '-v7.3');

fprintf('\nSaved output to:\n%s\n', outputFile);

%% ================================================================
% VERIFY SAVED WINDOW DIMENSIONS
%% ================================================================

firstValid = find( ...
    ~cellfun(@isempty, PLVwin_pre_med), ...
    1, 'first');

if ~isempty(firstValid)

    fprintf('\nSubject %s PRE meditation dimensions:\n', ...
        validSubs{firstValid});

    disp(size(PLVwin_pre_med{firstValid}));

    fprintf('Expected: 12 x 12 x 2 x nWindows\n');

end

fprintf('\n===== EXTRACTION COMPLETE =====\n');

%% ================================================================
% FUNCTION: COMPUTE RID-PLV FOR ALL WINDOWS
%% ================================================================

function [PLVmat, plv_windows, nWin] = ...
    compute_split_alpha_subjectavg_plv( ...
        data, win_samples, step_samples, n_freqs, bandBins)

    nChan = size(data,1);
    nBands = numel(bandBins);

    start_idx = 1:step_samples: ...
        (size(data,2) - win_samples + 1);

    nWin = numel(start_idx);

    if nWin < 1

        PLVmat = nan(nChan, nChan, nBands);
        plv_windows = nan(nChan, nChan, nBands, 0);

        return;
    end

    plv_windows = nan(nChan, nChan, nBands, nWin);

    for w = 1:nWin

        idx = start_idx(w): ...
            (start_idx(w) + win_samples - 1);

        xw = data(:,idx);

        W = nan(nChan, nChan, nBands);

        for i = 1:nChan

            for j = i+1:nChan

                data1 = xw(i,:);
                data2 = xw(j,:);

                % RID-Rihaczek distributions
                Rx = rid_rihaczek4(data1, n_freqs);
                Ry = rid_rihaczek4(data2, n_freqs);

                % Normalize to obtain phase
                Rx_phase = Rx ./ (abs(Rx) + eps);
                Ry_phase = Ry ./ (abs(Ry) + eps);

                % Phase difference
                phdif = Rx_phase .* conj(Ry_phase);

                % PLV at each frequency
                % Average phase difference over time
                plv_f = abs(mean(phdif, 2));

                % Average frequency-specific PLVs
                % within each alpha sub-band
                for b = 1:nBands

                    bins = bandBins{b};

                    band_plv = mean( ...
                        plv_f(bins), 'omitnan');

                    W(i,j,b) = band_plv;
                    W(j,i,b) = band_plv;

                end

            end
        end

        % Diagonal remains NaN.
        % Only off-diagonal electrode pairs are assigned.

        % Store each individual window
        plv_windows(:,:,:,w) = W;

    end

    % Average across windows for each participant
    PLVmat = mean(plv_windows, 4, 'omitnan');

end

%% ================================================================
% FUNCTION: PROGRESS DISPLAY
%% ================================================================

function printProgress(msg)

    persistent t0 finishedCount

    if isempty(t0)
        t0 = tic;
        finishedCount = 0;
    end

    elapsedMin = toc(t0) / 60;

    if startsWith(msg, 'FINISHED')

        finishedCount = finishedCount + 1;

        fprintf( ...
            '[%.1f min] %s | completed subjects: %d\n', ...
            elapsedMin, msg, finishedCount);

    else

        fprintf('[%.1f min] %s\n', elapsedMin, msg);

    end

end
