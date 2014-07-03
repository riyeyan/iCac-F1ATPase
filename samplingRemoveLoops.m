function samplingData = samplingRemoveLoops(model,nPoints)

changeCobraSolver('gurobi5','LP');
changeCobraSolver('gurobi5','MILP');
changeCobraSolver('gurobi5','MIQP');


options.nFiles = 10;
options.nWarmupPoints=5000;

ratio = 1.2;
options.nPointsReturned = nPoints*ratio;
options.removeLoopsFlag = true;      
options.removeLoopSamplesFlag = true;

[modelSampling,samples] = sampleCbModel(model,'CacModelSamples','ACHR',options);
% sampleFluxVecs = zeros(length(modelSampling.rxns),options.nPointsReturned);
% difference = zeros(options.nPointsReturned,1);
% for i = 1:options.nPointsReturned
%     [sampleFluxVecs(:,i), difference(i)] = nearestLoopLessFlux(modelSampling, samples(:,i), 'MILP');
%     if i > 5
%         i
%     end
% end
[sampleFluxVecs, difference] = nearestLoopLessFlux(modelSampling, samples, 'MILP');
ind = all(sampleFluxVecs==0,1);
sampleFluxVecs(:,ind)=[];difference(ind)=[];
sampleFluxVecs = sampleFluxVecs(:,1:nPoints); difference = difference(1:nPoints);

[liA,locB]=ismember(modelSampling.rxns,model.rxns);

v = zeros(length(model.rxns),nPoints);
v(locB,:)=sampleFluxVecs;





samplingData.options = options;
samplingData.reducedModel = modelSampling;
samplingData.reducedSamples = samples;
samplingData.rxnIndex = locB;
samplingData.samples = v;
samplingData.differences = difference;

% for i = 1:10
%     delete(strcat('KoxModelSamples_',num2str(i)));
% end

% save samplingData samplingData;

end


