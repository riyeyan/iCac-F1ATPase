function der = calDerBySpline(x,y,methodFlag)
% Yan Zhu
% 10/1/14
% edit on 19/2/14


%{
% k=4; kn=3;
% kn = augknt(x([1 end]),k);
% kn = 1;

figure=1;

% sp = spap2(kn,k, x,y);
% sp = spap2(kn,k,x,y,weights);
% sp = csapi(x,y);

% smoothing spline
sp = csaps(x,y,-1,[],weights);
% sp = spapi(7,x,y);
spd=fnder(sp);
der = fnval(spd,x);
if figure > 0
    plot(x,y,'*r');
    hold on;
    fnplt(sp,'-b');
    hold off;
end

%}

[fitresult, gof] = createFit(x,y,methodFlag);

der = differentiate(fitresult,x);

return

function [fitresult, gof] = createFit(x, y, f)
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
if f == 2
    ft = 'splineinterp';
    opts = fitoptions( ft );
    opts.Normalize = 'on';
elseif f == 1
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( ft );
elseif f == 0
    ft = fittype( 'pchipinterp' );
    opts = fitoptions( ft );
    opts.Normalize = 'on';
    
end

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

