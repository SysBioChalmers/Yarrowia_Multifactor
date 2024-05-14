adapterLocation = fullfile(findGECKOroot,'ecCollab','ecCollabAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();
fluxData = loadFluxData();
model = loadConventionalGEM();
allSolutions = struct();  % Initialize allSolutions as a structure

for i = 1:length(fluxData.conds)

    modelName = ['ec',fluxData.conds{i}];
    fprintf('Working on: %s \n', modelName)

    fileName = [modelName,'_prot.yml'];

    % Load model: already constrained with 10% variance around chemostat fluxes
    ecModel = loadEcModel(fileName);

    % Set bounds for biomass too
    sol = solveLP(ecModel);
    ecModel = setParam(ecModel,'var', 'xBIOMASS', -sol.f, 10);

    % Minimize all protein usage
    % All draw from pool
    usageRxnsIdx = startsWith(ecModel.rxns,'usage_prot_');
    protPoolIdx = find(ismember(ecModel.mets,'prot_pool'));
    ecModel.S(protPoolIdx, usageRxnsIdx) = 1;

    % Reset pool to all enzymes
    ecModel = setProtPoolSize(ecModel,fluxData.Ptot(i));

    % Minimize total enzyme usage
    ecModel = setParam(ecModel,'obj','prot_pool_exchange',1);
    sol = solveLP(ecModel);
    ecModel = setParam(ecModel,'var', 'prot_pool_exchange', -sol.f, 10);

    % Get good reactions
    [~, goodRxns] = randomSampling(ecModel,1,true,true,true);

    % Run random sampling
    solutions = randomSampling(ecModel,5000,true,true,true,goodRxns);

    % Get index of C-source rxn
    idx = getIndexes(ecModel,params.c_source,'rxns');

    % Compute absolute value of the element in the specified row
    abs_value = abs(solutions(idx, :));
    
    % Divide each column by the absolute value
    solutionsNorm = solutions ./ abs_value;

    % Store solutions in a field named after the model name
    allSolutions.(modelName) = solutionsNorm;
    allSol4Plot.(modelName) = mapRxnsToConv(ecModel,model,solutionsNorm);

    % Data treatment
    fluxMean = full(mean(solutionsNorm,2));
    fluxSD = full(std(solutionsNorm,0,2));
    logVec = abs(fluxMean) > abs(fluxSD);
    fluxMeanSD = logVec.*fluxMean;

    % Export data
    fluxesForEscher(ecModel.rxns,fluxMean,[modelName,'_prot_FBA_SD.json']);
    fluxesForEscher(model.rxns,full(mean(mapRxnsToConv(ecModel,model,solutionsNorm),2)),[modelName,'_FBA_SD.json']);
end

%%
% Calculate means of the fluxes for each reaction across all conditions
mean_fluxes.ecST6512 = mean(full(allSol4Plot.ecST6512), 2);
mean_fluxes.ecOKYL029 = mean(full(allSol4Plot.ecOKYL029), 2);
mean_fluxes.ecJFYL07 = mean(full(allSol4Plot.ecJFYL07), 2);
mean_fluxes.ecJFYL14 = mean(full(allSol4Plot.ecJFYL14), 2);
mean_fluxes.ecJFYL18 = mean(full(allSol4Plot.ecJFYL18), 2);

% Calculate stddevs for all conditions
standardDev.ecST6512 = std(full(allSol4Plot.ecST6512), 0, 2);
standardDev.ecOKYL029 = std(full(allSol4Plot.ecOKYL029), 0, 2);
standardDev.ecJFYL07 = std(full(allSol4Plot.ecJFYL07), 0, 2);
standardDev.ecJFYL14 = std(full(allSol4Plot.ecJFYL14), 0, 2);
standardDev.ecJFYL18 = std(full(allSol4Plot.ecJFYL18), 0, 2);

% Calculate fold changes (or differences in means)
fold_changes = mean_fluxes.ecJFYL18 ./ mean_fluxes.ecJFYL07;

% Perform two-sample t-test to compute p-values
[~, p_values] = ttest2(full(allSol4Plot.ecJFYL18), full(allSol4Plot.ecJFYL07), 'Dim', 2);

% Adjust p-values using Benjamini and Hochberg method
adjusted_p_values = mafdr(p_values,'BHFDR',true);

fluxCutOff = 0.00001;

fold_changes = fold_changes.*(abs(mean_fluxes.ecJFYL18) > fluxCutOff).*(abs(mean_fluxes.ecJFYL18) > fluxCutOff);

% Plot volcano plot
figure;
scatter(fold_changes, adjusted_p_values, 'ko'); % Plot p-values
set(gca, 'XScale', 'log'); % Set x-axis to log scale
set(gca, 'YScale', 'log'); % Set y-axis to log scale
set(gca, 'YDir', 'reverse'); % Reverse the y-axis direction
hold on;

% Highlight significant points (if desired)
alpha = 0.001; % significance level
significant_points = adjusted_p_values < alpha;
scatter(fold_changes(significant_points), adjusted_p_values(significant_points), 'r*'); % Highlight significant points

% Customize plot
xlabel('Fold Change');
ylabel('Adjusted p-value');
title('Volcano Plot of Flux Differences');
legend('All points', 'Significant points (adjusted p < 0.01)', 'Location', 'best');
grid on;
%%
% TO-DO
% - Get 10 biggest fold changes (p < 0.01)
% - Get 10 smalles fold changes (p < 0.01)
%   - Match them with model.rxnNames
% - Get 0, Inf and negative flux changes 
%   - match them with model.rxnNames

% Find indices of significant fold changes
significant_fold_changes = find(adjusted_p_values < 0.001);
% Get the corresponding fold changes
significant_fold_changes_values = fold_changes(significant_fold_changes);
% Sort fold changes in descending order to get the biggest ones
[sorted_fold_changes, sorted_indices] = sort(significant_fold_changes_values, 'descend');
% Select the 10 biggest fold changes
top_10_biggest_fold_changes = sorted_fold_changes(1:min(10, length(sorted_indices)));
% Get the names of the reactions corresponding to the top 10 biggest fold changes
top_10_biggest_fold_changes_names = model.rxns(significant_fold_changes(sorted_indices));

pValues_sorted = adjusted_p_values(significant_fold_changes(sorted_indices));

% Sort fold changes in ascending order to get the smallest ones
[sorted_fold_changes, sorted_indices] = sort(significant_fold_changes_values, 'ascend');
% Select the 10 smallest fold changes
top_10_smallest_fold_changes = sorted_fold_changes(1:min(10, length(sorted_indices)));
% Get the names of the reactions corresponding to the top 10 smallest fold changes
top_10_smallest_fold_changes_names = model.rxnNames(sorted_indices);

