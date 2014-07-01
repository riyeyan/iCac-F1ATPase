function ngam = calNGAM(model,ngamRxn)
% calculate non-growth ATP maintenance based on averaging the min/max
% values of ngam
% Yan Zhu
% 1/7/2014

[minNGAM, maxNGAM] = fluxVariability(model, 100, 'max', ngamRxn);
ngam = 0.5*(minNGAM + maxNGAM);

end

