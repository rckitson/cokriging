clear all
close all
clc

addpath('dace')

% EI, Kriging-Regpoly example

a = 1; b = 100;
lf = @(x) 0.5*asin(1./x);
hf = @(x) asin(1./(x)) - 1;

ncreate = 12;
nsamp = 5*ncreate;
lb = 0;
ub = pi/2;

figure(1)
xx = linspace(lb,ub);

plot(xx,lf(xx),'r')
hold on
plot(xx,hf(xx),'b')
% ylim([0 0.5]);
xlim([lb ub]);
axis off
grid off

set(gcf,'papersize',[4 3],'paperposition',[0 0 4 3])
saveas(gcf,'krigingExample_1.pdf')


clc;
xxlf = linspace(lb,ub,nsamp)';
xxi = randi(size(xxlf,1),ncreate,1);
xx = sort(dsmerge([lb; xxlf(xxi); ub],ones(ncreate+2,1)));
% xx = linspace(lb,ub,ncreate)';


%%
figure(2)
clf;
% Co-Kriging
plot(xxlf,lf(xxlf),'r-');
hold on
scatter(xx,hf(xx),8^2,'bs','filled');
plot(xxlf,hf(xxlf),'b')
% ylim([0 0.5]);
xlim([lb ub]);
axis off 
grid off

ss = sort(lhsamp(nsamp,1)*(ub - lb) + lb);
sk = zeros(size(ss));
for ii = 1:length(ss)
    [ds, ~] = dsmerge([xx; ss(ii)],ones(length(xx)+1,1),1e-6);
    if length(ds) > length(xx)
        sk(ii) = 1;
    else
        warning('weeding out points')
    end
end

ss = ss(logical(sk));

[dmodel, dmc, dmd] = cokriging2(ss,real(lf(ss)), xx, real(hf(xx)),@regpoly0,@corrgauss,1e-6,1e2);
krig = dacefit(xx,real(hf(xx)),@regpoly0,@corrgauss,1e-6,1e-3,20);

ss = gridsamp([lb;ub],100);
cok = predict_cok2(ss, dmodel);
kg = predictor(ss,krig);

cok_plot = sortrows([ss, cok],1);
kg_plot = sortrows([ss, kg],1);
plot(cok_plot(:,1), cok_plot(:,2),'k')
plot(kg_plot(:,1), kg_plot(:,2),'k--')

% legend('High-Fidelity','Low-Fidelity','Samples','Co-Kriging','Reconstruction','Location','SouthEast')

y1 = hf(ss(:,1));
y2 = cok;
err(3) = sqrt(mean((y2 - y1).^2))/mean(y1);

axis off
grid off

set(gcf,'papersize',[4 3],'paperposition',[0 0 4 3])
saveas(gcf,'krigingExample_3.pdf')

return

%%
figure(6)
clf;
bar(err)
set(gca,'yscale','log')
ylabel('Normalized RMSE')
set(gca,'xTickLabel',{'Kriging','Co-Kriging','Hierarchical Kriging','Hierarchical-Co-Kriging'})
set(gcf,'papersize',2*[4 3],'paperposition',2*[0 0 4 3])
saveas(gcf,'krigingExample_6.pdf')

