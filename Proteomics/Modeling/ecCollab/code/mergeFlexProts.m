function newStruct = mergeFlexProts(flexProt, flexProt2)
    % Check if flexProt2 is provided and non-empty
    if nargin < 2 || isempty(flexProt2)
        % If flexProt2 is not provided or empty, return flexProt
        newStruct = flexProt;
        return;
    end
    
    % Combine the uniprotIDs, oldConcs, and flexConcs from both structures
    all_uniprotIDs = [flexProt.uniprotIDs; flexProt2.uniprotIDs];
    all_oldConcs = [flexProt.oldConcs; flexProt2.oldConcs];
    all_flexConcs = [flexProt.flexConcs; flexProt2.flexConcs];
    all_frequence = [flexProt.frequence; flexProt2.frequence];

    % Get unique Uniprot IDs
    unique_uniprotIDs = unique(all_uniprotIDs);

    % Initialize newStruct fields as empty cell arrays or arrays
    newStruct.uniprotIDs = {};
    newStruct.oldConcs = [];
    newStruct.flexConcs = [];
    newStruct.frequence = [];
    newStruct.ratioIncr = [];

    % Iterate over each unique Uniprot ID
    for i = 1:numel(unique_uniprotIDs)
        % Current Uniprot ID
        ID = unique_uniprotIDs{i};

        % Indices of the current Uniprot ID in both structures
        idx_flexProt = strcmp(flexProt.uniprotIDs, ID);
        idx_flexProt2 = strcmp(flexProt2.uniprotIDs, ID);

        % Extract relevant data for the current Uniprot ID
        oldConcs = [flexProt.oldConcs(idx_flexProt); flexProt2.oldConcs(idx_flexProt2)];
        flexConcs = [flexProt.flexConcs(idx_flexProt); flexProt2.flexConcs(idx_flexProt2)];
        frequence = [flexProt.frequence(idx_flexProt); flexProt2.frequence(idx_flexProt2)];

        % Calculate properties for the current Uniprot ID
        minOldConcs = min(oldConcs);
        maxFlexConcs = max(flexConcs);
        totalFrequence = sum(frequence);
        recalculatedRatioIncr = maxFlexConcs / minOldConcs;

        % Append the properties to newStruct
        newStruct.uniprotIDs = [newStruct.uniprotIDs; {ID}];
        newStruct.oldConcs = [newStruct.oldConcs; minOldConcs];
        newStruct.flexConcs = [newStruct.flexConcs; maxFlexConcs];
        newStruct.frequence = [newStruct.frequence; totalFrequence];
        newStruct.ratioIncr = [newStruct.ratioIncr; recalculatedRatioIncr];
    end

    % Sort newStruct by descending order of ratioIncr
    [~, idx_sorted] = sort([newStruct.ratioIncr], 'descend');
    newStruct.uniprotIDs = newStruct.uniprotIDs(idx_sorted);
    newStruct.oldConcs = newStruct.oldConcs(idx_sorted);
    newStruct.flexConcs = newStruct.flexConcs(idx_sorted);
    newStruct.frequence = newStruct.frequence(idx_sorted);
    newStruct.ratioIncr = newStruct.ratioIncr(idx_sorted);
end