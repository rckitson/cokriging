function [ dmodel, dmc, dmd ] = cokriging2( sc0, yc0, se0, ye0, regr, corr, lb, ub)
% CoKriging

smean = mean([sc0; se0]);
sstd = std([sc0; se0]);
ymean = mean([yc0; ye0]);
ystd = std([yc0; ye0]);

sc = (sc0 - repmat(smean,size(sc0,1),1))./repmat(sstd,size(sc0,1),1);
se = (se0 - repmat(smean,size(se0,1),1))./repmat(sstd,size(se0,1),1);
yc = (yc0 - repmat(ymean,size(yc0,1),1))./repmat(ystd,size(yc0,1),1);
ye = (ye0 - repmat(ymean,size(ye0,1),1))./repmat(ystd,size(ye0,1),1);

dmodel.smean = smean;
dmodel.sstd = sstd;
dmodel.ymean = ymean;
dmodel.ystd = ystd;

opts = optimset('fmincon');
opts.Display = 'off';
opts2 = optimset('fminbnd');
opts2.Display = 'off';

nc = length(sc);
ne = length(se);

th0 = 10^(mean(log10([ub, lb])))*ones(1,size(lb,2));
dmc = dacefit(sc,yc,regr,corr, th0, lb, ub);
fprintf('dmc.theta = % .2g\n', dmc.theta)
yc_e = predictor(se,dmc);

p0 = yc_e\ye;
[u0,s0,v0] = svd(p0);
p = u0*v0';

d = ye - yc_e*p;
dmd = dacefit(se,d,regr,corr, ub, lb, ub);
fprintf('dmd.theta = % .2g\n', dmd.theta)

dmodel.p = p;
dmodel.dmc = dmc;
dmodel.dmd = dmd;


function [cost] = getp(se,yc_e,ye,lb,ub,p,opts)

    d = ye - p*yc_e;
    dmd = dacefit(se, d, @regpoly0, @corrgauss, lb, lb, ub);
    psidee = dmd.C;
    cost = ne/2*log(dmd.sigma2) + 1/2*log(abs(det(psidee)));

    
    
    
    