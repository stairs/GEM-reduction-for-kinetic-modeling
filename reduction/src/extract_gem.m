% Context-Specific GEM Extraction Script
%
% This script extracts a HT-29 cell line-specific GEM from a generic
% human GEM using DepMap RNASeq data.
%
% Usage:
%   1. Ensure a .env file exists in the project root directory.
%   2. Specify the path to the Human-GEM repo root in the .env file:
%       HUMAN_GEM_PATH=path/to/HumanGEM
%   3. Run the script.
%
% Output:
% HT-29 cell line-specific GEM saved in SBML format.
%

rootDir = fileparts(fileparts(which(mfilename)));
humanGEMPath = dotenv(fullfile(rootDir,'.env')).env.HUMAN_GEM_PATH;

modelPath = fullfile(humanGEMPath,'model','Human-GEM.mat');
humanGEM = load(modelPath).ihuman;

convertGenes = false;
essentialTasksFilePath = fullfile(humanGEMPath,'data','metabolicTasks','metabolicTasks_Essential.txt');
rxnsFilePath = fullfile(humanGEMPath,'model','reactions.tsv');
prepData = prepHumanModelForftINIT(humanGEM, convertGenes, essentialTasksFilePath, rxnsFilePath);
prepData.essentialRxns = [prepData.essentialRxns; {'MAR03964'}];

ht29rnaSeq = readtable(fullfile(rootDir,'data','ht29-depmap-rnaseq-tpm.csv'));
transcriptionData.genes = ht29rnaSeq.gene;
transcriptionData.tissues = ht29rnaSeq.Properties.VariableNames(2:end);
transcriptionData.levels = ht29rnaSeq.HT29;
transcriptionData.threshold = 1;

condition = 'HT29';
initSteps = getHumanGEMINITSteps('1+0');
removeGenes = false;
useScoresForTasks = true;
ht29GEM = ftINIT(prepData, condition, [], [], ...
    transcriptionData, {}, initSteps, removeGenes, useScoresForTasks);

ht29GEM.id = 'HT29 GEM';
ht29GEM.name = ht29GEM.id;
ht29GEM.description = 'Genome-scale metabolic model for HT-29 cell line';
exportModel(ht29GEM, fullfile(rootDir,'out','HT29-GEM.xml'));

closedModelPath = fullfile(rootDir,'out','HT29-GEM-closed.xml');
exportModel(closeModel(ht29GEM), closedModelPath);
