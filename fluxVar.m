% calcualte flux bias and capacities based on flux variability analysis and
% solution space sampling
% Yan Zhu
% 3/7/14


% clear all
% clc
load simulationData data;
load 'Ca_iYZ766_30-Jun-2014';
fluxDist = data.fluxDist;

glcFlux = reshape(abs(fluxDist(findRxnIDs(model,'EX_glc(e)'),:)),4,2);
close all;

% based on FVA
v_avg = (data.minFluxDist+data.maxFluxDist)/2;
l_sol = abs(data.minFluxDist-data.maxFluxDist);

% based on sampling
vm = zeros(4,2,length(model.rxns));
vs = vm;


for i = 1:4
    for j = 1:2
        %         points = data.samples{i,j}.points;
        points = data.samples{i,j}.samples;
        vm(i,j,:) = mean(points,2);
        vs(i,j,:) = std(points,0,2);
    end
end

vmw = vm/5;
vm1 = vm./repmat(glcFlux,[1,1,length(model.rxns)]);
vs1 = vs./repmat(glcFlux,[1,1,length(model.rxns)]);
fluxDist_ratio = vm(:,1,:)./vm(:,2,:);


% calculate parameters
% P1 IMP1l; F1 ITF1

% Respiratory quotient
% rq = -reshape(vm(:,:,findRxnIDs(model,'EX_co2(e)')),[6,1])./reshape(vm(:,:,findRxnIDs(model,'EX_o2(e)')),[6,1]);
% Biomass yield
carbonSource = reshape(vm(:,:,findRxnIDs(model,{'EX_glc(e)','EX_ac(e)','EX_but(e)'})),[2*4,3]);
carbonSource_ind = carbonSource<0;
Ybiomass = reshape(vm(:,:,findRxnIDs(model,'BIOMASS')),[2*4,1])./-sum(carbonSource.*carbonSource_ind,2);
% ATP, NADH, NADPH yield
metIDs = {'atp[c]','nadh[c]','nadph[c]'}';
Ycof = zeros(3,2*4);
for i = 1:3
    
    metRxns = findRxnsFromMets(model,metIDs{i});
    metRxnInds = findRxnIDs(model,metRxns);
    metFluxes = reshape(vm(:,:,metRxnInds),[2*4,length(metRxnInds)])';
    metTurnover = abs(model.S(findMetIDs(model,metIDs{i}),metRxnInds))*abs(metFluxes)/2;
    Ycof(i,:) = -metTurnover./reshape(vm(:,:,findRxnIDs(model,'EX_glc(e)')),[1,8]);
end

% redox ratio
% rr = (reshape(vm(:,:,findRxnIDs(model,'GLYK')),[6,1])+reshape(vm(:,:,findRxnIDs(model,'GLYCDx')),[6,1]))./reshape(vm(:,:,findRxnIDs(model,'GLYCDH')),[6,1]);


% sheet = zeros(length(model.rxns),12);
% for i = 1:3
%     for j = 1:2
%         
%         
%         sheet(:,2*i-1+6*(2-j)) = reshape(vm(i,j,:),[length(model.rxns),1]);
%         sheet(:,2*i+6*(2-j)) = reshape(vs(i,j,:),[length(model.rxns),1]);
%     end
% end
% csvwrite('sheet.csv',sheet);

