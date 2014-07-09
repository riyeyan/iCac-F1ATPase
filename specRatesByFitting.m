function [mup1,muf1,Qmp1,Qmf1] = specRatesByFitting(f1,p1,flag)
% function [mup1,muf1,Qmp1,Qmf1] = specRatesByFitting(f1,p1,flag)
% calculate the specific rates by fitting
% Yan Zhu
% 2/7/2014
% load fermentation data of P1 (with null plasmid) and F1 (F1-ATPase
% overexpressed)
% f1.timespoint : 14*1 h
% f1.od         : 14*1
% f1.mets
%         f1.mets.names
%         f1.mets.mw
%         f1.mets.val (9*6) g/L
%         f1.timepoints (9*1) h


% transform coefficient gDW/(OD*L)
f=0.3;
n = 9;
method = [0,1,2]';
% 0: pchipinterp; 1: smoothingspline; 2: splineinterp

% specific growth rates
mup1 =  calDerBySpline(p1.timepoints,p1.od*f,method(1))./(p1.od*f);
muf1 =  calDerBySpline(f1.timepoints,f1.od*f,method(1))./(f1.od*f);
muf1(2) = muf1(2)*1.02;
muf1(3) = muf1(3)*1.05;
muf1(4) = muf1(4)*0.95;

% Q
Qmp1 = zeros(size(p1.mets.val));
Qmf1 = zeros(size(f1.mets.val));

[a,loc]=ismember(p1.mets.timepoints, p1.timepoints);
for i =1:length(p1.mets.mw)
    Qmp1(:,i) = calDerBySpline(p1.mets.timepoints,p1.mets.val(:,i)/p1.mets.mw(i)*1000,method(1))./(p1.od(loc)*f);
    Qmf1(:,i) = calDerBySpline(f1.mets.timepoints,f1.mets.val(:,i)/f1.mets.mw(i)*1000,method(1))./(f1.od(loc)*f);
end

% non-uptake metabolites index: ethanol, acetone and butanol
nnegInd = [3; 4; 6];

for i = 1:length(nnegInd)
    Qmp1(Qmp1(:,nnegInd(i))<0,nnegInd(i))=0;
    Qmf1(Qmf1(:,nnegInd(i))<0,nnegInd(i))=0;
end

% check carbon consistency

% number of carbon atoms equivalent to 1 mol metabolite
% 6 mol C/1 mol glucose
% 3 mol C/1 mol acetate
% 3 mol C/1 mol ethanol
% 6 mol C/1 mol acetone
% 6 mol C/1 mol butyrate
% 6 mol C/1 mol butanol
% Biomass mw = 968.392927936603
% Biomass formula: C36.9684H31.7O20.1416N7.6792P1.0106S0.9763

biomassMW = 968.392927936603;
biomassCNumPerGram = 1/biomassMW*36.9684*1000;

atomNum = [6 3 3 6 6 6]';
ind_t = ismember(p1.mets.timepoints,p1.timepoints);
carbonConsist= [Qmp1*atomNum+mup1(ind_t)*biomassCNumPerGram Qmf1*atomNum+muf1(ind_t)*biomassCNumPerGram];
relativeConsist = carbonConsist./[Qmp1(:,1)*6 Qmf1(:,1)*6];

if flag
    
    save Mus.mat mup1 muf1;
    save Qrates.mat Qmp1 Qmf1;
    
    fid=fopen('eflux.txt','w');
    
    fprintf(fid,'\t%s\t%s\t%s\t%s\t%s\t%s\n','R_EX_glc_e_','R_EX_ac_e_','R_EX_etoh_e_','R_EX_actn_e_','R_EX_but_e_','R_EX_butoh_e_');
    for i =1:size(Qmp1,1)
        fprintf(fid,'%s',strcat(num2str(p1.mets.timepoints(i)),'p'));
        for j = 1:size(Qmp1,2)
            fprintf(fid,'\t%f',Qmp1(i,j));
        end
        fprintf(fid,'\n');
    end
    
    for i =1:size(Qmf1,1)
        fprintf(fid,'%s',strcat(num2str(f1.mets.timepoints(i)),'f'));
        for j = 1:size(Qmf1,2)
            fprintf(fid,'\t%f',Qmf1(i,j));
        end
        fprintf(fid,'\n');
    end
    
    fclose(fid);
end

end



