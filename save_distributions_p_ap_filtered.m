% get distributions of microvascular parameters after application of Lambda
% condition.

clear
close all

currdir = pwd;
k = strfind(currdir,'\');
homepath = currdir(1:k(end));

nBins = 30;
L_lim = [0.9 1.6];

% load healthy unfiltered results
load('output_healthy_fluid_unfiltered_parameters.mat');
id = find(and(L>L_lim(1),L<L_lim(2)));

% convert parameter structure of simulations
ct = 1;
for i = id
    k_i.Val(ct) = varpar(i).k_i;
    k_a.Val(ct) = varpar(i).k_a;
    l_pv.Val(ct) = varpar(i).l_pv;
    n_v.Val(ct) = varpar(i).n_v;
    n_l.Val(ct) = varpar(i).n_l;
    q_l.Val(ct) = varpar(i).q_l;
    ct = ct+1;
end

figure('position',[50 50 400 400],'color','w');
hold on;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
xlabel('$k_i ({\rm m}^2)$','Interpreter','latex');
h_k_i = histogram(1e6*k_i.Val,nBins,'Normalization','probability');
k_i.Edges = 1e-6*h_k_i.BinEdges;
k_i.CumCounts = [0 cumsum(h_k_i.Values)]; k_i.CumCounts(end) = 1;
k_i.base = makedist('PiecewiseLinear', 'x', k_i.Edges, 'Fx', k_i.CumCounts);

figure('position',[150 50 400 400],'color','w');
hold on;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
xlabel('$k_a ({\rm m}^2)$','Interpreter','latex');
h_k_a = histogram(k_a.Val,nBins,'Normalization','probability');
k_a.Edges = h_k_a.BinEdges;
k_a.CumCounts = [0 cumsum(h_k_a.Values)]; k_a.CumCounts(end) = 1;
k_a.base = makedist('PiecewiseLinear', 'x', k_a.Edges, 'Fx', k_a.CumCounts);

figure('position',[250 50 400 400],'color','w');
hold on;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
xlabel('$l_{p,v}$ (cm/s/mmHg)','Interpreter','latex');
h_l_pv = histogram(l_pv.Val,nBins,'Normalization','probability');
l_pv.Edges = h_l_pv.BinEdges;
l_pv.CumCounts = [0 cumsum(h_l_pv.Values)]; l_pv.CumCounts(end) = 1;
l_pv.base = makedist('PiecewiseLinear', 'x', l_pv.Edges, 'Fx', l_pv.CumCounts);

figure('position',[350 50 400 400],'color','w');
hold on;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
h_n_v = histogram(n_v.Val,nBins,'Normalization','probability');
xlabel('$N_v ({\rm mm}^{-2})$','Interpreter','latex');
n_v.Edges = h_n_v.BinEdges;
n_v.CumCounts = [0 cumsum(h_n_v.Values)]; n_v.CumCounts(end) = 1;
n_v.base = makedist('PiecewiseLinear', 'x', n_v.Edges, 'Fx', n_v.CumCounts);

figure('position',[450 50 400 400],'color','w');
hold on;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
h_n_l = histogram(n_l.Val,nBins,'Normalization','probability');
xlabel('$N_{\ell} ({\rm mm}^{-2})$','Interpreter','latex');
n_l.Edges = h_n_l.BinEdges;
n_l.CumCounts = [0 cumsum(h_n_l.Values)]; n_l.CumCounts(end) = 1;
n_l.base = makedist('PiecewiseLinear', 'x', n_l.Edges, 'Fx', n_l.CumCounts);

figure('position',[550 50 400 400],'color','w');
hold on;
set(gca,'fontsize',15,'TickLabelInterpreter','latex','XScale','log');
h_q_l = histogram(log10(q_l.Val),nBins,'Normalization','probability');
xlabel('$\tilde{q}_{\ell} ({\rm m}^2/{\rm s})$','Interpreter','latex');
q_l.Edges = h_q_l.BinEdges;
q_l.CumCounts = [0 cumsum(h_q_l.Values)]; q_l.CumCounts(end) = 1;
q_l.base = makedist('PiecewiseLinear', 'x', q_l.Edges, 'Fx', q_l.CumCounts);


save(['filtered_Lambda_' num2str(L_lim(1)) '_' num2str(L_lim(2)) '.mat'],'k_i','k_a','l_pv','n_v','n_l','q_l','nBins');
