% added a test comment for Git
% test continuous functions for permeabilities and microvascular exchanges

clear
close all

r = linspace(2-0.34,2+0.54+1,1000);
r12 = 2;
r23 = 2.54;

k1 = 1e-18; k2 = 1e-17; k3 = 1e-15;

k12 = @(r) k2 + (k1-k2)./(1+(r/r12).^200);
k = @(r) k3 + (k12(r)-k3)./(1+(r/r23).^200);

figure;hold on;set(gca,'YScale','log');
plot(r,k(r),'linew',1.5);