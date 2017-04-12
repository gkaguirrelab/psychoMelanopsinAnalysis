% psychMelAnalysis_main
%
% This routine loads the results of psychophysical measurement of the
% perceptual properties of pulses of spectral modulation


%% Housekeeping
clear variables
close all
clc

%% Set directory and filename paths, machine dependent.
%
% Could move these to a preference set in a local hook, but
% for now just dealing with it.
[~, localHostName] = system('scutil --get LocalHostName');
[~, userName] = system('whoami');
localHostName = strtrim(localHostName);
userName = strtrim(userName);
switch (localHostName)
    case 'eagleray'
        % DHB's desktop
        dropboxDir = fullfile(filesep,'Volumes','Users1','Dropbox (Aguirre-Brainard Lab)');
        
    otherwise
        % Some unspecified machine, go with more typical default.
        dropboxDir = fullfile('/Users', userName, '/Dropbox (Aguirre-Brainard Lab)');
end
dataDir = '/MELA_data/MaxPulsePsychophysics/';
analysisDir = '/MELA_analysis/psychoMelanopsinAnalysis/';
figureDir = fullfile(dropboxDir,analysisDir,'figures');
if (~exist(figureDir,'dir'))
    mkdir(figureDir);
end

%% Subject list
subjectIDs={...
    'MELA_0074',...
    'MELA_0087',...
    'MELA_0089',...
    'MELA_0026',...
    'MELA_0082',...
    'MELA_0038',...
    'MELA_0094',...
    'MELA_0096',...
    'MELA_0088',...
    'MELA_0079',...
    'MELA_0043',...
    'MELA_0073',...
    'MELA_0080',...
    'MELA_0090',...
    'MELA_0037',...
    'MELA_0049',...
    'MELA_0050',...
    'MELA_0075',...
    'MELA_0077',...
    'MELA_0081',...
    };

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
        warning(['Subject ' subjectIDs{ss} ' does not have a .mat file']);
    else
        if length(matFileIdx)>1
            warning(['Subject ' subjectIDs{ss} ' has more than one .mat file. Taking the last in the list.']);
            matFileIdx=matFileIdx(end);
            dataCellArray{ss}=load(fileList{matFileIdx});
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

%% Sort the rows in order by repTag, then stimLabel, then trialNumber
for ss=1:length(subjectIDs)
    subjectTable=dataTable{ss};
    subjectTable=sortrows(subjectTable,{'repetitionTag','stimLabel','trialNum'});
    dataTable{ss}=subjectTable;
end


%% Calculate for each subject the correlation between rep1 and rep2 responses
for ss=1:length(subjectIDs)
    subjectTable=dataTable{ss};
    rep1Idx=find(subjectTable.repetitionTag=='rep1');
    rep2Idx=find(subjectTable.repetitionTag=='rep2');
    withinSubReliability(ss)=corr(subjectTable.response(rep1Idx),subjectTable.response(rep2Idx),'type','Spearman');
end

%% Calculate the difference between first and second responses for each subject
for ss=1:length(subjectIDs)
    subjectTable=dataTable{ss};
    rep1Idx=find(subjectTable.repetitionTag=='rep1');
    rep2Idx=find(subjectTable.repetitionTag=='rep2');
    differenceValues{ss}=(subjectTable.response(rep1Idx) - subjectTable.response(rep2Idx));
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
% resultTableBySubject
% resultTableByStimulus

%% Write summary tables to excel files
outputFileName=fullfile(dropboxDir, analysisDir, 'resultTableBySubject.csv');
writetable(resultTableBySubject,outputFileName);
outputFileName=fullfile(dropboxDir, analysisDir, 'resultTableByStimulus.csv');
writetable(resultTableByStimulus,outputFileName);

%% Get data into format for pca, svm
perceptualDimensions = unique(dataTable{1}.perceptualDimension,'stable');
nPerceptualDimensions = length(perceptualDimensions);
nSubjects = length(subjectIDs);
dataVectorsLightFlux1 = zeros(nSubjects,nPerceptualDimensions);
dataVectorsLightFlux2 = zeros(nSubjects,nPerceptualDimensions);
dataVectorsLMS1 = zeros(nSubjects,nPerceptualDimensions);
dataVectorsLMS2 = zeros(nSubjects,nPerceptualDimensions);
dataVectorsMel1 = zeros(nSubjects,nPerceptualDimensions);
dataVectorsMel2 = zeros(nSubjects,nPerceptualDimensions);
for ss=1:length(subjectIDs)
    subjectTable = dataTable{ss};
    rep1Sel = subjectTable.repetitionTag=='rep1';
    rep2Sel = subjectTable.repetitionTag=='rep2';
    lightFluxSel = subjectTable.stimLabel=='Light Flux';
    LMSSel = subjectTable.stimLabel=='MaxLMS';
    MelSel = subjectTable.stimLabel=='MaxMel';
    for dd = 1:nPerceptualDimensions
        dimSel = subjectTable.perceptualDimension == perceptualDimensions(dd);
        dataVectorsLightFlux1(ss,dd) = subjectTable.response(rep1Sel & lightFluxSel & dimSel);
        dataVectorsLightFlux2(ss,dd) = subjectTable.response(rep2Sel & lightFluxSel & dimSel);
        dataVectorsLMS1(ss,dd) = subjectTable.response(rep1Sel & LMSSel & dimSel);
        dataVectorsLMS2(ss,dd) = subjectTable.response(rep2Sel & LMSSel & dimSel);
        dataVectorsMel1(ss,dd) = subjectTable.response(rep1Sel & MelSel & dimSel);
        dataVectorsMel2(ss,dd) = subjectTable.response(rep2Sel & MelSel & dimSel);
    end
end
dataVectorsLightFluxMean = (dataVectorsLightFlux1 + dataVectorsLightFlux2)/2;
dataVectorsLMSMean = (dataVectorsLMS1 + dataVectorsLMS2)/2;
dataVectorsMelMean = (dataVectorsMel1 + dataVectorsMel2)/2;


%% Classification analysis
%
% Set up data
theResponses = [dataVectorsLightFluxMean ; dataVectorsLMSMean ; dataVectorsMelMean];
labelsLightFluxClassify = ones(size(dataVectorsLightFluxMean,1),1);
labelsLMSClassify = ones(size(dataVectorsLMSMean,1),1);
labelsMelClassify = 2*ones(size(dataVectorsLMSMean,1),1);
labelsLightFluxSep = ones(size(dataVectorsLightFluxMean,1),1);
labelsLMSSep = 2*ones(size(dataVectorsLMSMean,1),1);
labelsMelSep = 3*ones(size(dataVectorsLMSMean,1),1);
theLabelsClassify = [labelsLightFluxClassify; labelsLMSClassify ; labelsMelClassify];
theLabelsSep = [labelsLightFluxSep; labelsLMSSep ; labelsMelSep];

% Cross-validated linear SVM performance on PCA representation
nFolds = 10;
CVO = cvpartition(theLabelsClassify,'kfold',nFolds);
testCheckIndex = [];
thePrediction = NaN*ones(size(theLabelsClassify));
keepDim = 2;
for cc = 1:CVO.NumTestSets
  % Split set
    trainingIndex = CVO.training(cc);
    testIndex = CVO.test(cc);
    testCheckIndex = [testCheckIndex ; find(testIndex)];
    
    % Do PCA
    thePCABasis{cc} = pca(theResponses(trainingIndex,:));
    pcaResponses{cc} = (thePCABasis{cc}\theResponses')';
    
    % Do the classify
    classifyInfo{cc} = fitcsvm(pcaResponses{cc}(trainingIndex,1:keepDim),theLabelsClassify(trainingIndex),'KernelFunction','linear','Solver','SMO');
    thePrediction(testIndex) = predict(classifyInfo{cc},pcaResponses{cc}(testIndex,1:keepDim));
end
testCheckIndex = sort(testCheckIndex);
if (any(testCheckIndex ~= (1:length(theLabelsClassify))'))
    error('We do not understand the cvpartition object');
end
crossValCorrectMel = thePrediction == theLabelsClassify;
percentCorrectLMS_Mel = 100*sum(crossValCorrectMel)/length(crossValCorrectMel);

%% Classifier on all the data, for looking at boundary
%
% Do PCA
thePCABasisAll = pca(theResponses);
pcaResponsesAll = (thePCABasisAll\theResponses')';
classifyInfoAll = fitcsvm(pcaResponsesAll(:,1:keepDim),theLabelsClassify,'KernelFunction','linear','Solver','SMO');
w1 = classifyInfoAll.Beta(1);
w2 = classifyInfoAll.Beta(2);
b = classifyInfoAll.Bias;
fitX = linspace(min(pcaResponsesAll(:,1)),max(pcaResponsesAll(:,1)),100);
fitY = (-w1/w2)*fitX - b/w2;
discrim = w1*thePCABasisAll(:,1) + w2*thePCABasisAll(:,2);
discrim = discrim/norm(discrim);

%% Boostrap the classifier on all data, to get error bars on the dimensions
nSubjects = size(dataVectorsLMSMean,1);
nBootstraps = 1000;
for bb = 1:nBootstraps
    bindex = randi(nSubjects,nSubjects,1);
    theResponsesBoot = [dataVectorsLightFluxMean(bindex,:) ; dataVectorsLMSMean(bindex,:) ; dataVectorsMelMean(bindex,:)];
    labelsLightFluxClassifyBoot = ones(size(dataVectorsLightFluxMean,1),1);
    labelsLMSClassifyBoot = ones(size(dataVectorsLMSMean,1),1);
    labelsMelClassifyBoot = 2*ones(size(dataVectorsLMSMean,1),1);
    theLabelsClassifyBoot = [labelsLightFluxClassifyBoot; labelsLMSClassifyBoot ; labelsMelClassifyBoot];
    
    thePCABasisBoot = pca(theResponsesBoot);
    pcaResponsesBoot = (thePCABasisBoot\theResponsesBoot')';
    classifyInfoBoot = fitcsvm(pcaResponsesBoot(:,1:keepDim),theLabelsClassifyBoot,'KernelFunction','linear','Solver','SMO');
    w1Boot = classifyInfoBoot.Beta(1);
    w2Boot = classifyInfoBoot.Beta(2);
    bBoot = classifyInfoBoot.Bias;
    fitXBoot = linspace(min(pcaResponsesBoot(:,1)),max(pcaResponsesBoot(:,1)),100);
    fitYBoot = (-w1Boot/w2Boot)*fitXBoot - bBoot/w2Boot;
    discrimBoot(:,bb) = w1Boot*thePCABasisBoot(:,1) + w2Boot*thePCABasisBoot(:,2);
    discrimBoot(:,bb) = discrimBoot(:,bb)/norm(discrimBoot(:,bb));
end
meanDiscrimBoot = mean(discrimBoot,2);
stderrDiscrimBoot = std(discrimBoot,[],2);

% Plot of classification
classifyFigure = figure; clf; hold on
set(gca,'FontName','Helvetica','FontSize',14);
set(gcf,'Position',[100 100 700 700]);
plot(pcaResponsesAll(theLabelsSep == labelsLightFluxSep(1),1),pcaResponsesAll(theLabelsSep == labelsLightFluxSep(1),2),'ko','MarkerSize',12,'MarkerFaceColor','k');
plot(pcaResponsesAll(theLabelsSep == labelsLMSSep(1),1),pcaResponsesAll(theLabelsSep == labelsLMSSep(1),2),'ro','MarkerSize',12,'MarkerFaceColor','r');
plot(pcaResponsesAll(theLabelsSep == labelsMelSep(1),1),pcaResponsesAll(theLabelsSep == labelsMelSep(1),2),'go','MarkerSize',12,'MarkerFaceColor','g');
plot(fitX,fitY,'k','LineWidth',2);
xlabel('PCA Dimension 1','FontSize',18);
ylabel('PCA Dimension 2','FontSize',18);
title({'LightFlux/LMS versus Mel' ; sprintf('Classification Accuracy %d%%',round(percentCorrectLMS_Mel)) ; ''},'FontSize',18);
legend({'LightFlux', 'LMS', 'Mel'},'Location','NorthWest','FontSize',14);
curdir = pwd;
cd(figureDir);
FigureSave('ClassifyAll.pdf',classifyFigure,'pdf');
cd(curdir);

% Plot discriminant weights
barColor = [0.5 0.5 0.5];
barEdgeColor = [0 0 0];
errorBarColor = [0.5 0.5 0.5];
xColor = [0 0 0];
[~,index] = sort(meanDiscrimBoot);
discrimFigure = figure; clf; hold on
set(gcf,'Position',[100 100 1750 750]);
h = bar(1:length(perceptualDimensions),discrimBoot(index));
set(h,'FaceColor',barColor,'EdgeColor',barEdgeColor);
errorbar(1:length(perceptualDimensions),discrimBoot(index),stderrDiscrimBoot(index),'b.','Color',errorBarColor);
%plot(1:length(perceptualDimensions),discrim(index),'x','MarkerSize',16,'MarkerFaceColor',xColor,'MarkerEdgeColor',xColor);
perceptualDimensionsCell{1} = '';
for ll = 2:length(perceptualDimensions)+1
    perceptualDimensionsCell{ll} = char(perceptualDimensions(index(ll-1)));
end
perceptualDimensionsCell{length(perceptualDimensions)+2} = '';
set(gca,'XTickLabel',perceptualDimensionsCell,'FontSize',14);
xlabel('Perceptual Dimension','FontSize',16);
ylabel('Melanopishness','FontSize',16);
ylim([-0.75 0.75]);
set(gca,'YTick',[-0.75 -0.5 -0.25 0 0.25 0.5 0.75]);
title({'LMS versus Mel' ; 'Dimensional Interpreation from Classifier' ; ''},'FontSize',18);
curdir = pwd;
cd(figureDir);
FigureSave('DimInterpretAll.pdf',discrimFigure,'pdf');
cd(curdir);



