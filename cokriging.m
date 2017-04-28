function [ dmodel, dmc, dmd ] = cokriging( sc0, yc0, se0, ye0, regr, corr, lb, ub)
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

% th = fminbnd(@(th)gettheta(sc,yc,th),lb,ub);
% th = fmincon(@(th)gettheta(sc,yc,th),ub,[],[],[],[],lb,ub,[],opts);
% fprintf('dmc.theta = %g\n', th);
dmc = dacefit(sc,yc,regr,corr, ub/2, lb, ub);
yc_e = predictor(se,dmc);

psiccc = getpsi(sc,sc,dmc.theta);
psicce = getpsi(sc,se,dmc.theta);
psicee = getpsi(se,se,dmc.theta);

p = fminbnd(@(p)getp(se,yc_e,ye,lb,ub,p,opts),0,1,opts2);
fprintf('p = %g\n', p);

d = ye - p*yc_e;
% th = fminbnd(@(th)gettheta(se,d,th),lb,ub);
th = fmincon(@(th)gettheta(se,d,th),ub,[],[],[],[],lb,ub,[],opts);
fprintf('dmd.theta = %g\n', th);
dmd.corr = corr;
dmd.theta = th;
psidee = getpsi(se,se,dmd.theta);

muc = (ones(1,nc)*(psiccc\yc))/(ones(1,nc)*(psiccc\ones(nc,1)));
mud = (ones(1,ne)*(psidee\d))/(ones(1,ne)*(psidee\ones(ne,1)));

sig2chat = (yc - muc*ones(nc,1))'*(psiccc\(yc - muc*ones(nc,1))/nc);
sig2dhat = (d - mud*ones(ne,1))'*(psidee\(d - mud*ones(ne,1))/ne);

sig2c = var(yc);
sig2d = var(d);

C = [sig2c*psiccc, p*sig2c*psicce;
    p*sig2c*psicce', p^2*sig2c*psicee + sig2d*psidee];

dmodel.mu = p*mean(yc) + mean(d); %(ones(1,nc+ne)*(C\[yc; ye]))/(ones(1,nc+ne)*(C\ones(nc+ne,1)));
dmodel.y_mu = [yc - mean(yc);ye - (p*mean(yc) + mean(d))];
dmodel.C = C;
dmodel.p = p;
dmodel.dmc = dmc;
dmodel.dmd = dmd;
dmodel.sc = sc;
dmodel.se = se;
dmodel.sig2c = sig2c;
dmodel.sig2d = sig2d;
dmodel.sig2dhat = sig2dhat;
dmodel.sig2chat = sig2chat;

function [cost] = getp(se,yc_e,ye,lb,ub,p,opts)

    d = ye - p*yc_e;
%     th = fminbnd(@(th)gettheta(se,d,th),lb,ub);
    th = fmincon(@(th)gettheta(se,d,th),ub,[],[],[],[],lb,ub,[],opts);

    ne = length(se);
    psidee = getpsi(se,se,th);
    mud = (ones(1,ne)*(psidee\d))/(ones(1,ne)*(psidee\ones(ne,1)));
    sig2d = (d - mud*ones(ne,1))'*(psidee\(d - mud*ones(ne,1))/ne);
    cost = ne/2*log(sig2d) + 1/2*log(abs(det(psidee)));

function [cost] = gettheta(sc,yc,th)

nc = length(sc);
psiccc = getpsi(sc,sc,th);
muc = (ones(1,nc)*(psiccc\yc))/(ones(1,nc)*(psiccc\ones(nc,1)));
sig2c = (yc - muc*ones(nc,1))'*(psiccc\(yc - muc*ones(nc,1))/nc);
cost = nc/2*log(sig2c) + 1/2*log(abs(det(psiccc)));
    
function [psi] = getpsi(x1,x2,theta)

n1 = size(x1,1);
n2 = size(x2,1);
temp = zeros(n1,n2);
for ii = 1:n1
    xi = x1(ii,:);
    for jj = 1:n2
        temp(ii,jj) = sum(-theta.*(xi - x2(jj,:)).^2);
    end
end
psi = exp(temp);

% function [ss] = scale(s)
% 
%     ms = mean(s);
%     sts = std(s);
%     ss = (s - repmat(ms,size(s,1),1))./repmat(sts,size(s,1),1);


    
    
    
    