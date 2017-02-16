
% psychMelAnalysis_main
%
% This routine loads the results of psychophysical measurement of the
% perceptual properties of pulses of spectral modulation


%% Housekeeping
clear variables
close all
clc

%% Set directory and filename paths
[~, userName] = system('whoami');
userName = strtrim(userName);
dropboxDir = ...
    fullfile('/Users', userName, '/Dropbox (Aguirre-Brainard Lab)');

dataDir = '/MELA_data/MaxPulsePsychophysics/';
analysisDir = '/MELA_analysis/psychoMelanopsinAnalysis/';

%% Subject list
subjectIDs={'MELA_0094','MELA_0049','MELA_0050','MELA_0075','MELA_0077','MELA_0081'};

%% Load the data
% Turn off the warnings that arise from loading the saved files without
% having instantiated an OL object
warning('off','MATLAB:load:cannotInstantiateLoadedVariable');
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off','MATLAB:class:EnumerableClassNotFound');

for ss=1:length(subjectIDs)
    subjectDir=fullfile(dropboxDir,dataDir,subjectIDs{ss});
    fileList = getAllFiles(subjectDir);
    matFileIdx=find(~cellfun(@isempty, strfind(fileList,'.mat')));
    if isempty(matFileIdx)
        warning(['Subject ' subjectsIDs{ss} ' does not have a .mat file']);
    else
        if length(matFileIdx)>1
        warning(['Subject ' subjectsIDs{ss} ' has more than one .mat file']);
        else
            dataCellArray{ss}=load(fileList{matFileIdx});
        end
    end    
end % Loop over subjects
warning('on','MATLAB:load:cannotInstantiateLoadedVariable');
warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('on','MATLAB:class:EnumerableClassNotFound');

%% Assemble the responses into a table with categorical variables
for ss=1:length(subjectIDs)
    subjectTable=struct2table(dataCellArray{ss}.data);
    subjectTable.stimLabel=categorical(subjectTable.stimLabel);
    subjectTable.perceptualDimension=categorical(subjectTable.perceptualDimension);
    repTag=table(categorical([repmat(1,1,27) repmat(2,1,27)],[1 2],{'rep1','rep2'})');
    repTag.Properties.VariableNames={'repetitionTag'};
    dataTable{ss}=[subjectTable repTag];
end

%% Calculate for each subject the correlation between rep1 and rep2 responses
for ss=1:length(subjectIDs)
    subjectTable=dataTable{ss};
    rep1Idx=find(subjectTable.repetitionTag=='rep1');
    rep2Idx=find(subjectTable.repetitionTag=='rep2');
    withinSubReliability(ss)=corr(subjectTable.response(rep1Idx),subjectTable.response(rep2Idx),'type','Spearman');
end

%% Average the first and second responses for each subject
for ss=1:length(subjectIDs)
    subjectTable=dataTable{ss};
    rep1Idx=find(subjectTable.repetitionTag=='rep1');
    rep2Idx=find(subjectTable.repetitionTag=='rep2');
    averageScores=mean([subjectTable.response(rep1Idx) subjectTable.response(rep2Idx)],2);
    subjectTable(28:end,:)=[];
    subjectTable.response=averageScores;
    subjectTable.repetitionTag=[];
    foldedDataTable{ss}=subjectTable;
end

%% Calculate the correlation of each subject to the average of the other subjects
for ss=1:length(subjectIDs)
    subjectIdx=find([1:1:length(subjectIDs)]~=ss);
    for jj=1:length(subjectIdx)
       responseMatrix(jj,:)=foldedDataTable{subjectIdx(jj)}.response;
    end
    averageResponses=median(responseMatrix);
    betweenSubConsistency(ss)=corr(foldedDataTable{ss}.response,averageResponses','type','Spearman');
end

%% Obtain the mean and SEM across subjects for each perceptual dimension crossed with each direction
for ss=1:length(subjectIDs)
    responseMatrix(ss,:)=foldedDataTable{ss}.response;
end
tmpTable=foldedDataTable{1};
resultTableByStimulus=tmpTable(:,[2 3]);
resultTableByStimulus.medianResponse = median(responseMatrix)';
resultTableByStimulus.iqrResponse = iqr(responseMatrix)';

%% Assemble the result table by subject
resultTableBySubject=cell2table(subjectIDs');
resultTableBySubject.Properties.VariableNames{1}='subjectID';
resultTableBySubject.withinSubReliability=withinSubReliability';
resultTableBySubject.betweenSubConsistency=betweenSubConsistency';

%% Dump the tables to the console
resultTableBySubject
resultTableByStimulus

%% Write tables to excel file
outputFileName=fullfile(dropboxDir, analysisDir, 'resultTableBySubject.csv');
writetable(resultTableBySubject,outputFileName);

outputFileName=fullfile(dropboxDir, analysisDir, 'resultTableByStimulus.csv');
writetable(resultTableByStimulus,outputFileName);
