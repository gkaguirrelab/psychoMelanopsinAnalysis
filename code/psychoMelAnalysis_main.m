
% psychMelAnalysis_main
%
% This routine loads the results of psychophysical measurement of the
% perceptual properties of pulses of spectral modulation


%% Housekeeping
clear variables
close all

[~, userName] = system('whoami');
userName = strtrim(userName);
dropboxDir = ...
    fullfile('/Users', userName, '/Dropbox (Aguirre-Brainard Lab)');

%% Set paths to surveys and output
dataDir = '/MELA_data/MaxPulsePsychophysics/';
analysisDir = '/MELA_analysis/psychoMelanopsinAnalysis/';

% Set the output filenames
outputFileName=fullfile(dropboxDir, analysisDir, 'FileNameHere');


gitInfo=GetGITInfo(getpref('surveyMelanopsinAnalysis', 'projectDir'));
notesText{4}=['Local code path: ' gitInfo.Path];
notesText{5}=['Remote code path: ' gitInfo.RemoteRepository{1}];
notesText{6}=['Revision: ' gitInfo.Revision];


