% using Spline fitting smoothed specific rates as external constraints for flux
% balance analysis
% Yan Zhu
% 10/1/14
%
clear
clc
close all;

% load model and constraints data
% load '..\GEM models\Ca_iYZ766_13-Mar-2014.mat';
% load '..\GEM models\Ca_iYZ766_16-Jun-2014.mat';
load 'Ca_iYZ766_30-Jun-2014.mat';
load 'f1.mat';
load 'p1.mat';
[mup1,muf1,Qmp1,Qmf1] = specRatesByFitting(f1,p1,0);

% change cobra solver
changeCobraSolver('gurobi5','LP');
% changeCobraSolver('glpk','LP');

% booundary reactions required to be constrained
boundRxns = {'EX_glc(e)','EX_ac(e)','EX_etoh(e)','EX_actn(e)','EX_but(e)',...
    'EX_butoh(e)'};
model = leaveEssentialUptakesCac(model);

% non-growth maintenance
atpm = zeros(5,2);
ngamRxn = {'ATPM'};

% results
data.fluxDist = zeros(length(model.rxns),8);

% variable flux ranges
minFluxDist = zeros(length(model.rxns),8);
maxFluxDist = minFluxDist;

% sampling resultes cell array
sampleStructs = cell(5,2);

for i = 2:5
    % wild type model    
    % index of timepoints in od measurement
    ind_t = find(p1.timepoints==p1.mets.timepoints(i));
    
    modelP1 = model;
    modelF1 = model;
    
    % set exchange flux bounds
    for j=1:6
        modelP1 = changeRxnBounds(modelP1,boundRxns(j),Qmp1(i,j),'b');
        modelF1 = changeRxnBounds(modelF1,boundRxns(j),Qmf1(i,j),'b');
    end
    
    % solve model    
    % 16/6/14 fix growth rate to experimental values
    if mup1(ind_t) >=0
        modelP1 = changeRxnBounds(modelP1,'BIOMASS',mup1(ind_t),'b');
    else
        modelP1 = changeRxnBounds(modelP1,'BIOMASS',0,'b');
    end
    if muf1(ind_t) >=0
        modelF1 = changeRxnBounds(modelF1,'BIOMASS',muf1(ind_t),'b');
    else
        modelF1 = changeRxnBounds(modelF1,'BIOMASS',0,'b');
    end
    solP1 = optimizeCbModel(modelP1); solF1 = optimizeCbModel(modelF1);
    data.fluxDist(:,i-1) = solP1.x;data.fluxDist(:,i+3) = solF1.x;
    atpm(i,:) = [calNGAM(modelP1,ngamRxn) calNGAM(modelF1,ngamRxn)];
    printFluxVector(model,solP1.x,true,true);fprintf('mus = %f\n\n',mup1(ind_t));
    printFluxVector(model,solF1.x,true,true);fprintf('mus = %f\n\n',muf1(ind_t));
    
    % range of flux variability
    [minFluxDist(:,i-1), maxFluxDist(:,i-1)] = fluxVariability(modelP1,100,'max');
    [minFluxDist(:,i+3), maxFluxDist(:,i+3)] = fluxVariability(modelF1,100,'max');
    
    % sampling
    sampleStructs{i-1,1}=samplingRemoveLoops(modelP1,5000);   
    sampleStructs{i-1,2}=samplingRemoveLoops(modelF1,5000);   
end


data.atpm = atpm;
data.minFluxDist = minFluxDist;
data.maxFluxDist = maxFluxDist;
data.samples = sampleStructs;
save simulationData data;
clear all
load simulationData data;
load 'Ca_iYZ766_30-Jun-2014.mat';
fluxDist = data.fluxDist;
fluxDist_ratio = fluxDist(:,5:8)./fluxDist(:,1:4);
maxFlux = abs(fluxDist(findRxnIDs(model,'EX_glc(e)'),:));

