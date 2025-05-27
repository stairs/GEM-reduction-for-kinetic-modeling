function [Fluxes, Conc, N, N_red, L, S, R, ParametersID, ParameterValue, MO, REQ, ntwkreactions, intracellular_index, Cnc, vDependentIdx, nDependent] = ReadSBML(modelFile, concFile, fluxFile, SD)
% ReadSBML  Generic SBML importer for SimBiology models
%   [Fluxes, Conc, N, N_red, L, S, R, ParametersID, ParameterValue, MO, REQ,
%    ntwkreactions, intracellular_index, Cnc, vDependentIdx, nDependent] = 
%    ReadSBML(modelFile, concFile, fluxFile)
%   Optional input:
%     SD           Cell array of conserved moiety species names (with or without
%                  compartment prefixes). If omitted or empty, sbioconsmoiety
%                  will be called to generate SD internally.

% Import SBML model
MO = sbmlimport(modelFile);

% Determine conserved species list (SD)
if nargin < 4 || isempty(SD)
    [~, SD, ~, ~, ~] = sbioconsmoiety(MO, 'link');
end

% Strip compartments and find dependent species indices
vctDependent = zeros(1, numel(SD));
for i = 1:numel(SD)
    name = SD{i};
    % remove known prefixes
    stripped = erase(name, {'CYT.', 'EXT.'});
    SD{i} = stripped;
    % locate in MO.Species
    idx = find(strcmp({MO.Species.Name}, stripped), 1);
    if ~isempty(idx)
        vctDependent(i) = idx;
    end
end
% Clean dependent index vector
vDependentIdx = unique(vctDependent(vctDependent>0));
nDependent   = numel(vDependentIdx);

% Classify species by compartment
intracellular_index = [];
extracellular_index = [];
for i = 1:numel(MO.Species)
    compName = MO.Species(i).Parent.Name;
    if contains(compName, 'EXT')
        extracellular_index(end+1) = i;
    else
        intracellular_index(end+1) = i;
    end
end

% Define independent species: exclude external + dependent
mask = unique([extracellular_index, vDependentIdx]);
indepSpecies = setdiff(intracellular_index, mask);

% Define network reactions indices
nr = numel(MO.Reactions);
ntwkreactions = 1:nr;

% Stoichiometric matrix and reduction
[Nmat, ~, ~] = getstoichmatrix(MO);
Nf = full(Nmat);
N = Nf(intracellular_index, ntwkreactions);
N_red = Nf(indepSpecies, ntwkreactions);
L = N / N_red;

% Read steady-state data
Cnc = readtable(concFile);
Flx = readtable(fluxFile);

% Metabolite and reaction names
S = {MO.Species(intracellular_index).Name};
R = {MO.Reactions(ntwkreactions).Name};

% Assign concentrations
ssc = NaN(1, numel(intracellular_index));
for j = 1:height(Cnc)
    specName = char(Cnc.Species(j));
    idx = find(strcmp({MO.Species.Name}, specName), 1);
    if ~isempty(idx)
        MO.Species(idx).Value = Cnc.Concentration(j);
    end
end
for k = 1:numel(intracellular_index)
    ssc(k) = MO.Species(intracellular_index(k)).Value;
end
Conc = ssc;

% Assign fluxes
ssf = NaN(1, numel(ntwkreactions));
for j = 1:height(Flx)
    rxnName = Flx.Reaction{j};
    idx = find(strcmp({MO.Reactions.Name}, rxnName), 1);
    if ~isempty(idx)
        MO.Reactions(idx).UserData = Flx.Flux(j);
    end
end
for k = 1:numel(ntwkreactions)
    ssf(k) = MO.Reactions(ntwkreactions(k)).UserData;
end
Fluxes = ssf;

% Clean kinetic laws and collect reaction expressions
for i = 1:numel(MO.Reactions)
    rate = MO.Reactions(i).ReactionRate;
    rate = erase(rate, {'EXT.', 'CYT.'});
    MO.Reactions(i).ReactionRate = rate;
end
REQ = cell(1, numel(ntwkreactions));
for i = 1:numel(ntwkreactions)
    REQ{i} = MO.Reactions(ntwkreactions(i)).ReactionRate;
end

% Retrieve parameters
ParametersID = {};
ParameterValue = [];
for p = 1:numel(MO.Parameters)
    ParametersID{end+1}   = MO.Parameters(p).Name;
    ParameterValue(end+1) = MO.Parameters(p).Value;
end
for i = ntwkreactions
    law = MO.Reactions(i).KineticLaw;
    if ~isempty(law)
        for param = law.Parameters
            ParametersID{end+1}   = param.Name;
            ParameterValue(end+1) = param.Value;
        end
    end
end

end
