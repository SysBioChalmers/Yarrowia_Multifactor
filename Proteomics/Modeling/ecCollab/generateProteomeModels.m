%% TO-DO:
% Combine the flexProt variables of each model and save them

verbose = false;
saveModels = true;

% Get model adapter
adapterLocation = fullfile(findGECKOroot,'ecCollab','ecCollabAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

% Load pooled model
ecModel = loadEcModel('eciYali.yml');

fluxData = loadFluxData();
fluxDataIT = loadFluxData(fullfile(params.path,'data','fluxData_itaconate.tsv'));
lipidNchainData = loadLipidNchainData();
protData = loadProtData([3, 3, 3, 3, 3]);
sharedDeletions = {'336u_EXP_1','336u_EXP_2','336u_1_EXP_1','336u_1_EXP_2','2884g','2884g_REV','2948g','2948g_REV'};

% Calculate GANnonPol for later
GAMnonPol = calculateGAM(ecModel);

for i = 1:length(fluxData.conds)
    % Update the biomass equation for each condition
    modelName = ['ec',fluxData.conds{i}];
    fprintf('Working on: %s \n', modelName)
    ecModel_new = ecModel;

    % TO-DO: update models genetic background
    if ismember(fluxData.conds{i}, {'JFYL07', 'JFYL14', 'JFYL18'})

        % Block flux through the deleted genes
        ecModel_new = setParam(ecModel_new,'eq',sharedDeletions,0);

        % Update lipid equation
        equations.mets = {'m1000','m1579','m1581','m1640','m1641','m1648','m1651','m1705','m1700','m1701','m1631','m1727'};
        equations.stoichCoeffs = {[-0.031129914,-0.000147006,-0.000448367,0,-0.004633618,-0.005292203,-0.001237787,-0.004996722,-0.014849039,-0.020933603,-0.002891601,1]};
        ecModel_new = changeRxns(ecModel_new, {'xLIPID'}, equations, 1);

        if ismember(fluxData.conds{i}, {'JFYL14'})
            % Add cytosolic itaconate
            metsToAdd.mets = {'m1841','m1842'};
            metsToAdd.metNames = {'itaconate','itaconate'};
            metsToAdd.compartments = {'c','e'};
            metsToAdd.unconstrained = [0, 1];
            ecModel_new = addMets(ecModel_new,metsToAdd);

            % Add heterologous genes
            genesToAdd.genes = {'ATEG09970','ATEG09971'};
            genesToAdd.geneShortNames = {'mttA','cad1'};
            ecModel_new = addGenesRaven(ecModel_new,genesToAdd);

            % Add non enzymetic reactions
            rxnsToAdd.rxns = {'AtMTT','ITACNce','EX_ITACN_OUT'};
            rxnsToAdd.equations = {'cis-aconitate[m] + (S)-malate[c] => cis-aconitate[c] + (S)-malate[m]',...
                'itaconate[c] => itaconate[e]',...
                'itaconate[e] =>'};
            rxnsToAdd.rxnNames = {'tricarboxylic acid transporter',...
                'itaconate transport',...
                'itaconate exchange'};
            rxnsToAdd.grRules = {'ATEG09970','',''};

            ecModel_new = addRxns(ecModel_new, rxnsToAdd, 3);

            % Add enzymetic reactions
            newRxns.rxns = {'AtCAD'};
            newRxns.rxnNames = {'cis-aconitate decarboxylase'};
            newRxns.equations = {'cis-aconitate[c] + H+_p+1[c] => itaconate[c] + carbon dioxide[c]'};
            newRxns.grRules = {'ATEG09971'};

            newEnzymes.enzymes = {'B3IUN8'};
            newEnzymes.genes = {'ATEG09971'};
            newEnzymes.mw = 52754;

            ecModel_new = addNewRxnsToEC(ecModel_new, newRxns, newEnzymes, ModelAdapter);

            % Constrain the kcat
            % Kcat estimated using DLKcat
            ecModel_new = setKcatForReactions(ecModel_new,'AtCAD',35.8357);
            ecModel_new = applyKcatConstraints(ecModel_new);

        elseif ismember(fluxData.conds{i}, {'JFYL18'})
            ecModel_new = setParam(ecModel_new,'eq',{'662','734_EXP_1','734_EXP_2'},0);
        end
    end

    % Change protein content in the biomass equation
    ecModel_new = scaleBioMassYali4(ecModel_new, 'protein', fluxData.Ptot(i));

    % Change lipid content in the biomass equation
    ecModel_new = scaleBioMassYali4(ecModel_new, 'lipid', lipidNchainData.Ltot(i));

    % Adjust fatty acid distribution
    ecModel_new = updateAcylPool(ecModel_new, lipidNchainData, i, ModelAdapter);

    % Balance out mass with carbohydrate content
    [X,~,C,~,~,~,~] = sumBioMassYali4(ecModel_new, false);

    delta = X - 1;  % difference to balance
    fC = (C - delta)/C;
    ecModel_new = rescalePseudoReaction(ecModel_new, 'carbohydrate', fC);

    % Recalculate GAM
    [~,~,ecModel_new] = calculateGAM(ecModel_new, GAMnonPol,true);

    % Rename model
    ecModel_new.id = [modelName,'_pooled'];

    % Save pooled model
    if saveModels == true
        saveEcModel(ecModel_new,[ecModel_new.id,'.yml']);
        saveEcModel(ecModel_new,[ecModel_new.id,'.xml']);
    end

    % Constrain model with proteomics data
    ecModel_new = fillEnzConcs(ecModel_new, protData, i);
    ecModel_new = constrainEnzConcs(ecModel_new);

    % Update protein pool
    ecModel_new = updateProtPool(ecModel_new,fluxData.Ptot(i));

    feasability = -1;

    flexProtMerged = [];

    while feasability < 0

        % Load flux data
        if ismember(fluxData.conds{i}, {'JFYL14'})
            ecModel_new = constrainFluxData(ecModel_new, fluxDataIT, i,'max','loose');
        else
            ecModel_new = constrainFluxData(ecModel_new, fluxData, i,'max','loose');
        end

        sol = solveLP(ecModel_new); % To observe if growth was reached.
        fprintf('Growth rate that is reached: %f /hour.\n', abs(sol.f))

        % This strain needs flexibilizeEnzConcs or it will not reach the
        % growth rate
        if ismember(fluxData.conds{i}, {'OKYL029'})
            [ecModel_new, flexProt] = flexibilizeEnzConcs(ecModel_new, fluxData.grRate(i), 10);
            flexProtMerged = mergeFlexProts(flexProt,flexProtMerged);
        end

        idx = getIndexes(ecModel_new,'xMAINTENANCE','rxns');
        NGAM = ecModel_new.lb(idx);
        feasability = 1;

        while feasability > 0
            NGAM = NGAM + 0.1;
            ecModel_new = setParam(ecModel_new,'lb','xMAINTENANCE',NGAM);
            sol = solveLP(ecModel_new,1); % To observe if growth was reached.
            if verbose == true
                fprintf('Growth rate that is reached: %f /hour. NGAM equals %.1f \n', abs(sol.f), NGAM)
            end
            feasability = sol.stat;
        end

        ecModel_new = setParam(ecModel_new,'lb','xMAINTENANCE',NGAM-0.1);

        % Flexibilize protein concentrations
        [ecModel_new, flexProt] = flexibilizeEnzConcs(ecModel_new, fluxData.grRate(i), 10);
        flexProtMerged = mergeFlexProts(flexProt,flexProtMerged);

        % Sometime NGAM is too high and stops growth rate from being
        % achievable. This should fix that.
        sol = solveLP(ecModel_new);
        if -sol.f < 0.975*fluxData.grRate(i)
            ecModel_new = setParam(ecModel_new,'lb','xMAINTENANCE',0);
            [ecModel_new, flexProt] = flexibilizeEnzConcs(ecModel_new, fluxData.grRate(i), 10);
            flexProtMerged = mergeFlexProts(flexProt,flexProtMerged);
        end

        % Check if model works with tight contraints
        % Load flux data
        if ismember(fluxData.conds{i}, {'JFYL14'})
            ecModel_new = constrainFluxData(ecModel_new, fluxDataIT, i,'max',10);
        else
            ecModel_new = constrainFluxData(ecModel_new, fluxData, i,'max',10);
        end
        ecModel_new = setParam(ecModel_new,'ub','1992',0);
        sol = solveLP(ecModel_new);


        if sol.stat == 1

            NGAM = 0;
            ecModel_new = setParam(ecModel_new,'lb','xMAINTENANCE',NGAM);
            sol = solveLP(ecModel_new);
            maxGrowth = round(-sol.f,4);

            while round(-sol.f,4) == maxGrowth
                NGAM = NGAM + 0.1;
                ecModel_new = setParam(ecModel_new,'lb','xMAINTENANCE',NGAM);
                sol = solveLP(ecModel_new); % To observe if growth was reached.
                if verbose == true
                    fprintf('Growth rate that is reached: %f /hour. NGAM equals %.1f \n', abs(sol.f), NGAM)
                end
            end

            ecModel_new = setParam(ecModel_new,'lb','xMAINTENANCE',NGAM-0.1);
            break
        end
        feasability = sol.stat;
    end

    % Rename model
    ecModel_new.id = [modelName,'_prot'];

    % Save proteome model
    if saveModels == true
        saveEcModel(ecModel_new,[ecModel_new.id,'.yml']);
    end
end