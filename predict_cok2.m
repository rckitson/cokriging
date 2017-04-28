function [ y, mse ] = predict_cok2( x0, dmodel )

smean = dmodel.smean;
sstd = dmodel.sstd;
ymean = dmodel.ymean;
ystd = dmodel.ystd;

y0 = zeros(size(x0,1),size(ymean,2));

x = (x0 - repmat(smean,size(x0,1),1))./repmat(sstd,size(y0,1),1);

yc = predictor(x, dmodel.dmc);
yd = predictor(x, dmodel.dmd);
y0 = reshape(yd, size(x,1), size(y0,2)) + reshape(yc, size(x,1), size(y0,2))*dmodel.p;

y = y0.*repmat(ystd,size(y0,1),1) + repmat(ymean,size(y0,1),1);



