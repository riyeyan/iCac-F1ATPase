% calcualte flux bias and capacities based on flux variability analysis and
% solution space sampling
% Yan Zhu
% 25/12/13


clear all
clc
load simulationData data;
load 'C:\Users\uqyanzhu\Documents\Kox_ORP\kox_model\Kox_iYZ1315.mat';
fluxDist = data.fluxDist;

maxFlux = reshape(abs(fluxDist(findRxnIDs(model,'EX_glyc(e)'),:)),3,2);
close all;

% based on FVA
v_avg = (data.minfluxDist+data.maxfluxDist)/2;
l_sol = abs(data.minfluxDist-data.maxfluxDist);

% based on sampling
vm = zeros(3,2,length(model.rxns));
vs = vm;


for i = 1:3
    for j = 1:2
        %         points = data.samples{i,j}.points;
        points = data.samples{i,j}.samples;
        vm(i,j,:) = mean(points,2);
        vs(i,j,:) = std(points,0,2);
    end
end

vmw = vm/5;
vm1 = vm./repmat(maxFlux,[1,1,length(model.rxns)]);
vs1 = vs./repmat(maxFlux,[1,1,length(model.rxns)]);
fluxDist_ratio = vm(:,1,:)./vm(:,2,:);


% calculate parameters in Table 1
% L-240 mV H-150 mV

% Respiratory quotient
rq = -reshape(vm(:,:,findRxnIDs(model,'EX_co2(e)')),[6,1])./reshape(vm(:,:,findRxnIDs(model,'EX_o2(e)')),[6,1]);
% Biomass yield
Ybiomass = reshape(vm(:,:,findRxnIDs(model,'Biomass')),[6,1])./reshape(vm(:,:,findRxnIDs(model,'EX_glyc(e)')),[6,1]);
% ATP, NADH, NADPH yield
metIDs = {'atp[c]','nadh[c]','nadph[c]'}';
Ycof = zeros(3,6);
for i = 1:3
    
    metRxns = findRxnsFromMets(model,metIDs{i});
    metRxnInds = findRxnIDs(model,metRxns);
    metFluxes = reshape(vm(:,:,metRxnInds),[6,length(metRxnInds)])';
    metTurnover = abs(model.S(findMetIDs(model,metIDs{i}),metRxnInds))*abs(metFluxes)/2;
    Ycof(i,:) = -metTurnover./reshape(vm(:,:,findRxnIDs(model,'EX_glyc(e)')),[1,6]);
end

% redox ratio
rr = (reshape(vm(:,:,findRxnIDs(model,'GLYK')),[6,1])+reshape(vm(:,:,findRxnIDs(model,'GLYCDx')),[6,1]))./reshape(vm(:,:,findRxnIDs(model,'GLYCDH')),[6,1]);


sheet = zeros(length(model.rxns),12);
for i = 1:3
    for j = 1:2
        
        
        sheet(:,2*i-1+6*(2-j)) = reshape(vm(i,j,:),[length(model.rxns),1]);
        sheet(:,2*i+6*(2-j)) = reshape(vs(i,j,:),[length(model.rxns),1]);
    end
end
csvwrite('sheet.csv',sheet);

