% Go to the models folder

% Get model adapter
adapterLocation = fullfile(findGECKOroot,'ecCollab','ecCollabAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

% Load model
model = loadConventionalGEM();

% Load flux data
fluxData = loadFluxData();

i = 1;
% j = 200;
% quitCondition = fluxData.grRate(i);
%
% while quitCondition == fluxData.grRate(i)
%     model = constrainFluxData(model, fluxData, i,'max',j);
%     model = setParam(model,'ub','1992',0);
%     model = setParam(model,'eq','xMAINTENANCE',9.4287);
%     sol = solveLP(model,1); % To observe if growth was reached.
%     %fprintf('Growth rate that is reached: %f /hour. Variance equals %.f percent \n', abs(sol.f), j)
%     j = j-1;
%     quitCondition = -sol.f;
% end

% Constrain upper bound of experimental data (except O2)
model = constrainFluxData(model, fluxData, i,'max','loose');
model = setParam(model,'ub','1992',0);

% Maximize and fix max biomass
sol = solveLP(model);
model = setParam(model,'lb','xBIOMASS',-sol.f);

% Calculate and fix max NGAM
model = setParam(model,'obj','xMAINTENANCE',1);
sol = solveLP(model);
model = setParam(model,'lb','xMAINTENANCE',-sol.f);

% Bring biomass back to the normal spot
model = setParam(model,'obj','xBIOMASS',1);
model = setParam(model,'lb','xBIOMASS',0);

k = 1.00;
feasability = 1;
idx = getIndexes(model,'xMAINTENANCE','rxns');
ogNGAM = model.lb(idx);

while feasability > 0
    k = k + 0.01;
    model = setParam(model,'lb','xMAINTENANCE',k*ogNGAM);
    sol = solveLP(model,1); % To observe if growth was reached.
    fprintf('Growth rate that is reached: %f /hour. Variance equals %f \n', abs(sol.f), k)
    feasability = sol.stat;
end

