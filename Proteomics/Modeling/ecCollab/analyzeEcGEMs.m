% Define adapter location and set default model adapter
adapterLocation = fullfile(findGECKOroot, 'ecCollab', 'ecCollabAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();
fluxData = loadFluxData();
model = loadConventionalGEM();
allSolutions = struct();
nSamples = 5000;

% Loop through each condition
for i = 1:length(fluxData.conds)
    modelName = ['ec', fluxData.conds{i}];
    fprintf('Working on: %s \n', modelName);

    % Load and constrain the model
    ecModel = loadEcModel([modelName, '_prot.yml']);
    sol = solveLP(ecModel);
    ecModel = setParam(ecModel, 'var', 'xBIOMASS', -sol.f, 10);
    
    % Set protein pool
    usageRxnsIdx = startsWith(ecModel.rxns, 'usage_prot_');
    protPoolIdx = find(ismember(ecModel.mets, 'prot_pool'));
    ecModel.S(protPoolIdx, usageRxnsIdx) = 1;
    ecModel = setProtPoolSize(ecModel, fluxData.Ptot(i));
    ecModel = setParam(ecModel, 'obj', 'prot_pool_exchange', 1);
    sol = solveLP(ecModel);
    ecModel = setParam(ecModel, 'var', 'prot_pool_exchange', -sol.f, 10);

    % Perform random sampling
    [~, goodRxns] = randomSampling(ecModel, 1, true, true, true);
    solutions = randomSampling(ecModel, nSamples, true, true, true, goodRxns);
    fluxes = mean(full(solutions), 2);
    usageData.(modelName) = enzymeUsage(ecModel, fluxes);

    % Normalize solutions
    idx = getIndexes(ecModel, params.c_source, 'rxns');
    solutionsNorm = solutions ./ abs(solutions(idx, :));

    % Store solutions and statistics
    allSolutions.(modelName) = solutionsNorm;
    allSol4Plot.(modelName) = mapRxnsToConv(ecModel, model, solutionsNorm);
    mean_Fluxes.(modelName) = mean(full(allSol4Plot.(modelName)), 2);
    standardDev.(modelName) = std(full(allSol4Plot.(modelName)), 0, 2);

    % Data treatment
    fluxMean = full(mean(solutionsNorm, 2));
    fluxSD = full(std(solutionsNorm, 0, 2));
    logVec = abs(fluxMean) > abs(fluxSD);
    fluxMeanSD = logVec .* fluxMean;

    % Export data
    fluxesForEscher(ecModel.rxns, fluxMean, [modelName, '_prot_FBA_SD.json']);
    fluxesForEscher(model.rxns, full(mean(mapRxnsToConv(ecModel, model, solutionsNorm), 2)), [modelName, '_FBA_SD.json']);
end

% Output to TSV
outputFileName = 'Supplementary_Table_fluxes_data.tsv';
outputFilePath = 'output'; % You can set this to the desired directory path
modelAdapter = ModelAdapterManager.getDefault(); % Assuming you have already set the default adapter

writeTSV(outputFileName, model.rxnNames, mean_Fluxes, standardDev, fluxData.conds, outputFilePath, modelAdapter);

% Define pairs of models to compare
modelPairs = {
    'ecOKYL029', 'ecST6512';
    'ecJFYL07', 'ecOKYL029';
    'ecJFYL14', 'ecJFYL07';
    'ecJFYL18', 'ecJFYL07';
};

% Calculate fold changes, p-values, and adjusted p-values
results = computeStatistics(modelPairs, mean_Fluxes, allSol4Plot, model);

% Filter, sort, and update results
results = filterAndSortResults(results, standardDev, modelPairs, mean_Fluxes, model);

% Analysis and output
analyzeResults(results, modelPairs, 'SupTable_foldChangeAnalysis.tsv', outputFilePath, modelAdapter);

%% Heatmap with Reduced Font Size in Cells
% Define specific reactions and their corresponding names
specificReactions = {'958', '714', '713', '300', '302', '280', '658', '216', '471', '217', '1889'}; % Reactions IDs
specificReactionNames = {'PC', 'MDHc', 'MDHm','CSm', 'ACONTa', 'ACONTb', ...
                         'ICDHx', 'ASPTAc', 'GLUDy', 'ASPTAm', 'EX_Glu'}; % Names for the reactions
numReactions = length(specificReactions);

% Get the model names from the mean_Fluxes struct
modelNames = fieldnames(mean_Fluxes);
numModels = length(modelNames);

% Initialize the matrix to store flux values
specificFluxMatrix = nan(numReactions, numModels);

% Loop through each specific reaction and fetch the flux values from mean_Fluxes
for i = 1:numReactions
    reaction = specificReactions{i};
    
    for j = 1:numModels
        modelName = modelNames{j};
        
        % Find the index of the specific reaction in the model
        rxnIndex = find(strcmp(model.rxns, reaction));
        
        if ~isempty(rxnIndex)
            % Store the mean flux value for the specific reaction and model
            specificFluxMatrix(i, j) = abs(mean_Fluxes.(modelName)(rxnIndex));
        else
            warning('Reaction %s not found in the model.', reaction);
        end
    end
end

% Create the heatmap
h = heatmap(modelNames, specificReactionNames, specificFluxMatrix); % Using the custom names for reactions

% Customize the heatmap
h.Title = 'Specific Metabolic Fluxes Heatmap with Reduced Font Size';
h.XLabel = 'Models';
h.YLabel = 'Specific Reactions';
h.Colormap = parula;  % You can choose other colormaps like jet, hot, cool, etc.
h.ColorbarVisible = 'on';

% Adjust color scaling if needed
clim([min(specificFluxMatrix(:)), max(specificFluxMatrix(:))]);  % Set the color axis scaling

% Reduce font size of heatmap labels
h.YDisplayLabels = specificReactionNames;
h.YDisplayLabels = specificReactionNames;
h.FontSize = 10; % Set font size

% Display the heatmap
disp('Heatmap with reduced font size in cells generated successfully.');


%% Functions
% Function to write TSV file
function writeTSV(fileName, rxnNames, mean_Fluxes, standardDev, conds, filePath, modelAdapter)
    if nargin < 7 || isempty(modelAdapter)
        modelAdapter = ModelAdapterManager.getDefault();
        if isempty(modelAdapter)
            error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
        end
    end
    params = modelAdapter.getParameters();

    if nargin < 6 || isempty(filePath)
        filePath = fullfile(params.path, 'output');
    end

    % Ensure the directory exists, create it if needed
    if ~exist(filePath, 'dir')
        mkdir(filePath);
    end

    % Full file path including directory
    fullFilePath = fullfile(filePath, fileName);

    % Open the file for writing
    fileID = fopen(fullFilePath, 'w');
    if fileID > 0
        fieldNames = fieldnames(mean_Fluxes);
        fprintf(fileID, 'model.rxnNames');
        for i = 1:length(fieldNames)
            fprintf(fileID, '\t%s_avg\t%s_std', fieldNames{i}, fieldNames{i});
        end
        fprintf(fileID, '\n');
        for j = 1:length(rxnNames)
            fprintf(fileID, '%s', rxnNames{j});
            for i = 1:length(conds)
                modelName = ['ec', conds{i}];
                fprintf(fileID, '\t%.8f\t%.8f', mean_Fluxes.(modelName)(j), standardDev.(modelName)(j));
            end
            fprintf(fileID, '\n');
        end
        fclose(fileID);
        disp(['TSV file created: ', fullFilePath]);
    else
        error('Unable to open the file for writing.');
    end
end

% Function to compute statistics (fold changes, p-values, adjusted p-values)
function results = computeStatistics(modelPairs, mean_Fluxes, allSol4Plot, model)
    results = struct();
    for k = 1:size(modelPairs, 1)
        model1 = modelPairs{k, 1};
        model2 = modelPairs{k, 2};
        pairName = [model1, '_vs_', model2];
        fold_changes = mean_Fluxes.(model1) ./ mean_Fluxes.(model2);
        [~, p_vals] = ttest2(full(allSol4Plot.(model1)), full(allSol4Plot.(model2)), 'Dim', 2);
        adj_p_vals = mafdr(p_vals, 'BHFDR', true);
        results.(pairName) = struct('rxns', {model.rxns}, 'rxnNames', {model.rxnNames}, ...
                                    'fold_changes', fold_changes, 'p_values', p_vals, ...
                                    'adjusted_p_values', adj_p_vals);
    end
end

% Function to filter, sort, and update results
function results = filterAndSortResults(results, standardDev, modelPairs, mean_Fluxes, model)
    pValCutoff = 0.001;
    xchangeRxnsIDs = getExchangeRxns(model);
    transpRxnsIDs = model.rxns(getTransportRxns(model));
    for k = 1:size(modelPairs, 1)
        model1 = modelPairs{k, 1};
        model2 = modelPairs{k, 2};
        pairName = [model1, '_vs_', model2];
        data = results.(pairName);
        fold_changes = data.fold_changes;
        adj_p_values = data.adjusted_p_values;
        flux_sd_model1 = standardDev.(model1);
        avg_flux_model1 = mean_Fluxes.(model1);

        % Filter fold changes
        fold_changes(adj_p_values > pValCutoff | flux_sd_model1 > avg_flux_model1) = NaN;
        isExchangeTransport = ismember(data.rxns, [xchangeRxnsIDs; transpRxnsIDs]);
        isPseudoRxn = contains(data.rxnNames, 'pseudoreaction');
        filtered_idx = ~isnan(fold_changes) & ~isExchangeTransport & ~isPseudoRxn;

        % Update filtered results
        filtered_results = structfun(@(field) field(filtered_idx), data, 'UniformOutput', false);
        [~, sort_idx] = sort(filtered_results.fold_changes, 'descend');
        results.(pairName) = structfun(@(field) field(sort_idx), filtered_results, 'UniformOutput', false);
    end
end

% Function to analyze results and output to a file
function analyzeResults(results, modelPairs, fileName, filePath, modelAdapter)
    if nargin < 5 || isempty(modelAdapter)
        modelAdapter = ModelAdapterManager.getDefault();
        if isempty(modelAdapter)
            error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
        end
    end
    params = modelAdapter.getParameters();

    if nargin < 4 || isempty(filePath)
        filePath = fullfile(params.path, 'output');
    end

    % Ensure the directory exists, create it if needed
    if ~exist(filePath, 'dir')
        mkdir(filePath);
    end

    % Full file path including directory
    fullFilePath = fullfile(filePath, fileName);

    header_row = 'Pair\tNum Increases\tAvg Increase\tMedian Increase\tNum Decreases\tAvg Decrease\tMedian Decrease\tNum Rev Changes\tTotal Sig Changes\n';
    fid = fopen(fullFilePath, 'w');  % Open in write mode
    if fid > 0
        fprintf(fid, header_row);
        for k = 1:size(modelPairs, 1)
            pairName = [modelPairs{k, 1}, '_vs_', modelPairs{k, 2}];
            data = results.(pairName);
            fold_changes = data.fold_changes;
            num_increases = sum(fold_changes > 1 & ~isinf(fold_changes));
            num_decreases = sum(fold_changes > 0 & fold_changes <= 1 & ~isinf(fold_changes));
            num_reversible_changes = sum(fold_changes < 0 & ~isinf(fold_changes));
            fold_changes_finite = fold_changes(~isinf(fold_changes));
            avg_increase = nanmean(fold_changes_finite(fold_changes_finite > 1));
            median_increase = nanmedian(fold_changes_finite(fold_changes_finite > 1));
            avg_decrease = nanmean(fold_changes_finite(fold_changes_finite > 0 & fold_changes_finite <= 1));
            median_decrease = nanmedian(fold_changes_finite(fold_changes_finite > 0 & fold_changes_finite <= 1));
            num_significant_changes = length(fold_changes);
            output_string = sprintf('%s\t%d\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t%d\t%d\n', ...
                                    pairName, num_increases, avg_increase, median_increase, ...
                                    num_decreases, avg_decrease, median_decrease, ...
                                    num_reversible_changes, num_significant_changes);
            fprintf(fid, output_string);
        end
        fclose(fid);
        disp(['TSV file created: ', fullFilePath]);
    else
        error('Unable to open the file for writing.');
    end
end