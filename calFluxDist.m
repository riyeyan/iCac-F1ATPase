% using B-spline smoothed specific rates as external constraints for flux
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
load '..\GEM models\Ca_iYZ766_30-Jun-2014.mat';
load '..\f1.mat';
load '..\p1.mat';
load Qrates.mat;
load Mus.mat;

% change cobra solver
changeCobraSolver('gurobi5','LP')

% booundary reactions required to be constrained
boundRxns = {'EX_glc(e)','EX_ac(e)','EX_etoh(e)','EX_actn(e)','EX_but(e)',...
    'EX_butoh(e)'};
model = leaveEssentialUptakesCac(model);

% prediction rerror rate vector
rerr = zeros(5,2);aerr=rerr;

% specific growth rates
expVal = zeros(5,2);
calVal = zeros(5,2);

% non-growth maintenance
atpm = zeros(5,2);
ngamRxn = {'ATPM'};
% flux distribution array
fluxDistributions = zeros(5,2,length(model.rxns));

% variable flux ranges
minfluxDist = zeros(5,2,length(model.rxns));
maxfluxDist = zeros(5,2,length(model.rxns));

% sampling resultes cell array
sampleStructs = cell(5,2);

% figure handle
h = figure;

for i = 2:5
    % wild type model
    
    % index of timepoints in od measurement
    ind_t = find(p1.timepoints==p1.mets.timepoints(i));
    
    % set exchange flux bounds
    for j=1:6
        model = changeRxnBounds(model,boundRxns(j),Qmp1(i,j),'b');
    end
    
    % solve model
    
    % 16/6/14 fix growth rate to experimental values
    if mup1(ind_t) >=0
        model = changeRxnBounds(model,'BIOMASS',mup1(ind_t),'b');
    else
        model = changeRxnBounds(model,'BIOMASS',0,'b');
    end
    sol = optimizeCbModel(model);
    fluxDistributions(i,1,:) = sol.x;
    atpm(i,1) = calNGAM(model,ngamRxn);
    % calcualte lower/upper limits of reaction flux variation
    %     [minfluxDist(i,1,:), maxfluxDist(i,1,:)] = fluxVariability(model,100,'max',model.rxns,true);
    
    % erros
    aerr(i,1)=abs(sol.f-mup1(ind_t));
    rerr(i,1) = aerr(i,1)/mup1(ind_t);
    
    printFluxVector(model,sol.x,true,true);
    
    fprintf('mus = %f\n\n',mup1(ind_t));
    calVal(i,1) = sol.f;
    expVal(i,1)=mup1(ind_t);
    
    % sampling
    %     model_s=changeRxnBounds(model,'BIOMASS',sol.f,'b');
    %     sampleStructs{i-1,1}=samplingRemoveLoops(model_s,5000);
    
    % mutant model
    
    % index of timepoints in od measurement
    ind_t = find(f1.timepoints==f1.mets.timepoints(i));
    
    % set exchange flux bounds
    for j=1:6
        model = changeRxnBounds(model,boundRxns(j),Qmf1(i,j),'b');
    end
    % calculate non-growth maintenance given experiemental derived specific
    % growth rate
    %     atpm(i-1,2) = calNGAM(model,muf1(ind_t));
    %
    %     % fix non-growth maintenance
    %     model = changeRxnBounds(model,'ATPM',atpm(i-1,2),'b');
    
    % solve model
    if muf1(ind_t)>=0
        model = changeRxnBounds(model,'BIOMASS',muf1(ind_t),'b');
    else
        model = changeRxnBounds(model,'BIOMASS',0,'b');
    end
    sol = optimizeCbModel(model);
    fluxDistributions(i,2,:) = sol.x;
    atpm(i,2) = calNGAM(model,ngamRxn);
    % calcualte lower/upper limits of reaction flux variation
    %     [minfluxDist(i,2,:), maxfluxDist(i,2,:)] = fluxVariability(model,100,'max',model.rxns,true);
    
    % erros
    aerr(i,2)=abs(sol.f-muf1(ind_t));
    rerr(i,2) = aerr(i,2)/muf1(ind_t);
    
    printFluxVector(model,sol.x,true,true);
    
    fprintf('mus = %f\n\n',muf1(ind_t));
    calVal(i,2) = sol.f;
    expVal(i,2)=muf1(ind_t);
    
    % sampling
    %     model_s=changeRxnBounds(model,'BIOMASS',sol.f,'b');
    %     sampleStructs{i-1,2}=samplingRemoveLoops(model_s,5000);
end


subplot(1,2,1);
bar([expVal(:,1) calVal(:,1)]);
axis([0 6 0 0.5]);
subplot(1,2,2);
bar([expVal(:,2) calVal(:,2)]);
axis([0 6 0 0.5]);

print(h,'C:\Users\uqyanzhu\Documents\F1-ATPase\F1ATPase manuscript\Figures\Figure_S2.eps','-depsc2');

data.accuracy = aerr;
data.relativeAccu = rerr;
data.fluxDist = zeros(length(model.rxns),8);
% 1-4: pIMP1
data.fluxDist(:,1) = reshape(fluxDistributions(1,1,:),length(model.rxns),1);
data.fluxDist(:,2) = reshape(fluxDistributions(2,1,:),length(model.rxns),1);
data.fluxDist(:,3) = reshape(fluxDistributions(3,1,:),length(model.rxns),1);
data.fluxDist(:,4) = reshape(fluxDistributions(4,1,:),length(model.rxns),1);
% 5-8, pITF1
data.fluxDist(:,5) = reshape(fluxDistributions(1,2,:),length(model.rxns),1);
data.fluxDist(:,6) = reshape(fluxDistributions(2,2,:),length(model.rxns),1);
data.fluxDist(:,7) = reshape(fluxDistributions(3,2,:),length(model.rxns),1);
data.fluxDist(:,8) = reshape(fluxDistributions(4,2,:),length(model.rxns),1);

data.atpm = atpm;
data.minfluxDist = minfluxDist;
data.maxfluxDist = maxfluxDist;
data.samples = sampleStructs;
save simulationData data;
clear all
load simulationData data;
load '..\GEM models\Ca_iYZ766_09-Mar-2014.mat';
fluxDist = data.fluxDist;
fluxDist_ratio = fluxDist(:,1:3)./fluxDist(:,4:6);
maxFlux = abs(fluxDist(findRxnIDs(model,'EX_glyc(e)'),:));
% close all;
