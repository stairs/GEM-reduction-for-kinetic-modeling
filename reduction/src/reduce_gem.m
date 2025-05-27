% HT29-GEM Reduction Script
%
% This script prunes and compresses ht29 genome-scale metabolic network
% while preserving the following functions:
%   • Growth on defined medium
%   • Anaerobic (fermentation) ATP yield
%   • Aerobic (respiration) ATP yield
%   • Growth without lactate excretion
%

rootDir = fileparts(fileparts(which(mfilename)));
closedModelPath = fullfile(rootDir,'model','HT29-GEM-closed.xml');
cnaModel = CNAsbmlModel2MFNetwork(closedModelPath);

% Block all exchange reactions
exchanges = getExchanges(cnaModel);
for i = 1:length(exchanges)
    rxnIndex = mfindstr(cnaModel.reacID, exchanges(i));
    cnaModel.lb(rxnIndex) = 0;
end

medium = readtable(fullfile(rootDir,'data','medium.csv'));
[mediumConstraintMatrix, mediumConstraintBoundary] = generateMediumConstraint(cnaModel, medium.rxn_id, medium.lb);

% Protected function #1: growth on the defined medium with a maximal growth 
% rate that matches the full model at a high L/G ratio
protectedFunctions(1).D = mediumConstraintMatrix;
protectedFunctions(1).d = mediumConstraintBoundary;

Bio_exc = mfindstr(cnaModel.reacID,'R_MAR10024');
protectedFunctions(1).D(end + 1, Bio_exc) = -1;
protectedFunctions(1).d(end + 1) = -0.06482 * 0.99;

O2T = mfindstr(cnaModel.reacID,'R_MAR09048');
protectedFunctions(1).D(end + 1, O2T) = -1;
protectedFunctions(1).d(end + 1) = 0.030071;

% Protected function #2: regenerate ATP from glucose under anaerobic 
% conditions (fermentation) with a yield that matches the full model 
protectedFunctions(2).D = zeros(0, cnaModel.numr);

ATPM = mfindstr(cnaModel.reacID,'R_MAR03964');
GLUT = mfindstr(cnaModel.reacID,'R_MAR09034');
protectedFunctions(2).D(1, ATPM)  = -1;
protectedFunctions(2).D(2, GLUT)  = -1;

protectedFunctions(2).d = [-2; 1];

% Protected function #3: regenerate ATP from glucose under aerobic 
% conditions (respiration) with a yield that matches the full model 
protectedFunctions(3).D = zeros(0, cnaModel.numr);
protectedFunctions(3).D(1, ATPM) = -1;
protectedFunctions(3).D(2, GLUT) = -1;
protectedFunctions(3).D(3, O2T) = -1;
protectedFunctions(3).d = [-157.45; 5; 1000];

% Protected function #4: growth on the defined medium with a maximal 
% growth rate that matches the full model and no lactate excretion
[mediumConstraintMatrix, mediumConstraintBoundary] = generateMediumConstraint(cnaModel, medium.rxn_id, medium.lb);

protectedFunctions(4).D = mediumConstraintMatrix;
protectedFunctions(4).d = mediumConstraintBoundary;

Bio_exc = mfindstr(cnaModel.reacID,'R_MAR10024');
protectedFunctions(4).D(end + 1, Bio_exc) = -1;
protectedFunctions(4).d(end + 1) = -0.06482 * 0.99;

LACT = mfindstr(cnaModel.reacID,'R_MAR09135');
protectedFunctions(4).D(end + 1, LACT) = 1;
protectedFunctions(4).d(end + 1) = 0;

% Generate a list of protected reaction indices
protected_rxn_ids = ["R_MAR04394" "R_MAR04379" "R_MAR04391" "R_MAR04368" "R_MAR04363" "R_MAR04358" "R_MAR04381"...
    "R_MAR04375" "R_MAR04373" "R_MAR04365" "R_MAR09135" "R_MAR09034" "R_MAR10024"];
protected_rxn_count = length(protected_rxn_ids);
protected_rxn_indices = zeros(1, protected_rxn_count);

for i = 1:protected_rxn_count
   rxn_index = mfindstr(cnaModel.reacID, protected_rxn_ids(i));
   protected_rxn_indices(i) = rxn_index;
end

% Define the remaining parameters for the model reduction
solver = 3; % GUROBI
min_degrees_of_freedom = 1;
protected_mets = [];

reactions_feasible = 1;
min_reactions = 1;
compression = 1;  % change this parameter to 0 to generate a pruned but not compressed model
rational = 0;
verbose = 1;

reducedCnaModel = CNAreduceMFNetwork(cnaModel, min_degrees_of_freedom, protectedFunctions, protected_mets, ...
    protected_rxn_indices, reactions_feasible, min_reactions, solver, compression, rational);

Bio_exc = mfindstr(reducedCnaModel.reacID,'R_MAR10024');
reducedCnaModel.objFunc = zeros(reducedCnaModel.numr, 1);
reducedCnaModel.objFunc(Bio_exc) = -1;

if compression    
    % Assign metabolites names in the reduced models as NetworkReducer
    % does not preserve after compression
    reducedCnaModel.specLongName = getSpeciesNamesByIds(cnaModel, reducedCnaMode.specID);

    % Make sure that reaction ids are not too long (maximum id length in Cobra is 255 characters)
    reducedCnaModel.reacID = normalizeReactionIds(reducedCnaModel.reacID);
    
    reducedCnaModel = CNAremoveConsRel(reducedCnaModel, false, true, true);
end

reducedCnaModel.reacNotes = repmat({''}, 1, reducedCnaModel.numr);

CNAMFNetwork2sbml(reducedCnaModel, fullfile(rootDir,'out','HT29-GEM-reduced.xml'));

function speciesNames = getSpeciesNamesByIds(model, speciesIds)    
    indices = zeros(length(speciesIds), 1);
    
    for i = 1:length(speciesIds)
        id = deblank(speciesIds(i,:));
        indices(i) = mfindstr(model.specID, id);
    end
    speciesNames = model.specLongName(indices, :);
end

function rxnIds = normalizeReactionIds(rxnIds)    
    for i = 1:size(rxnIds, 1)
        rxnId = deblank(rxnIds(i,:));
        
        if length(rxnId) > 255
            rxnIds(i,:) = blanks(size(rxnIds, 2));
            newRxnId = strjoin({'Reaction', num2str(i)});
            rxnIds(i, 1:numel(newRxnId)) = newRxnId;
        end        
    end
end

function [constraintMatrix, constraintBoundary] = generateMediumConstraint(cnaModel, rxnIds, lb)
    constraintMatrix = zeros(0, cnaModel.numr);
    constraintBoundary = zeros(length(rxnIds), 1);
    
    for i = 1:length(rxnIds)
        rxnIndex = mfindstr(cnaModel.reacID, rxnIds(i));
        
        if rxnIndex ~= 0
            constraintMatrix(i, rxnIndex) = -1;
            constraintBoundary(i) = lb(i);
        end
    end
end

function exchangeRxns = getExchanges(cnaModel)
    S = cnaModel.stoichMat;
    exchIdx = find(any(S(cnaModel.specExternal, :), 1));
    ids = string(cnaModel.reacID);
    exchangeRxns = ids(exchIdx);
end
