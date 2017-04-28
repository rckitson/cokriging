function [ y, mse ] = predict_cok( x0, dmodel )

y0 = zeros(size(x0,1),1);
mse = zeros(size(x0,1),1);

smean = dmodel.smean;
sstd = dmodel.sstd;

x = (x0 - repmat(smean,size(x0,1),1))./repmat(sstd,size(y0,1),1);

p = dmodel.p;
C = dmodel.C;
mu = dmodel.mu;
y_mu = dmodel.y_mu;
sc = dmodel.sc;
se = dmodel.se;

sig2c = dmodel.sig2c;
sig2d = dmodel.sig2d;

corrc = dmodel.dmc.corr;
corrd = dmodel.dmd.corr;
thetac = dmodel.dmc.theta;
thetad = dmodel.dmd.theta;

for ii = 1:size(x,1);
    nc = length(sc);
    ne = length(se);
    c = zeros(nc+ne,1);
    
    c(1:nc) = p*sig2c*corrc(thetac,sc - repmat(x(ii,:),nc,1));
    c(nc+1:end) = p^2*sig2c*corrc(thetac,se - repmat(x(ii,:),ne,1)) + ...
        sig2d*corrd(thetad,se - repmat(x(ii,:),ne,1));
    
    y0(ii) = mu + c'*(C\y_mu);
    mse(ii) = p^2*sig2c + sig2d - c'*(C\c);
end

ymean = dmodel.ymean;
ystd = dmodel.ystd;

y = y0.*repmat(ystd,size(y0,1),1) + repmat(ymean,size(y0,1),1);



