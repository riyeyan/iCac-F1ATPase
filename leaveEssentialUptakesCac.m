function modelE = leaveEssentialUptakesCac(model)
% useage: function modelE = leaveEssentialUptakes(model)
% change the lower bounds of unessential uptake reactions to zeros,
% specific to Clostridium acetobutylicum
% INPUT:
% model: the original genome-scale metabolic model
% modelE: changed genome-scale metabolic model;

[selExc,selUpt] = findExcRxns(model);
model = changeRxnBounds(model,model.rxns(selUpt),0,'l');
uptakeRxnIDs = {
    'EX_4abz(e)'
    'EX_btn(e)'
    'EX_ca2(e)'
    'EX_cd2(e)'
    'EX_cl(e)'
    'EX_co2(e)'
    'EX_cobalt2(e)'
    'EX_h2o(e)'
    'EX_k(e)'
    'EX_mg2(e)'
    'EX_mn2(e)'
    'EX_na1(e)'
%     'EX_no3(e)'
    'EX_pi(e)'
    'EX_h(e)'
    'EX_h2(e)'
    'EX_so4(e)'
    'EX_zn2(e)'
    'EX_nh4(e)'
    'EX_cu2(e)'
    'EX_fe2(e)'
    'EX_fe3(e)'};
uptakeLb = -1000+zeros(size(uptakeRxnIDs));
% shut the PFL flux
model = changeRxnBounds(model,'PFL',0,'u');
modelE = changeRxnBounds(model,uptakeRxnIDs,uptakeLb,'l');
end
        