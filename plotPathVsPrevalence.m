function plotPathVsPrevalence(matfile,approach,multiple_loci,for_paper)

addpath ./utils/
addpath ./utils/cbrewer

if nargin < 4
  for_paper = false;
end

switch multiple_loci
  case 0
    fld = 'singleloci';
  case 1
    fld = 'multiloci';
  case 2
    fld = 'ubloci';
end

appendix = [ '_' approach '_' fld];

load(matfile,'network');

net = network.(approach);
paths = net.(fld).unique_paths;
clusters = net.clusters.overlapping;
adj = net.(fld).adj_mat;

%% Get direction/severity of paths

U = length(paths);
m_i = zeros(U,1);
m_f = zeros(U,1);
ell = zeros(U,1);
for i = 1:length(paths)    
    m_i(i) = paths{i}(1);
    m_f(i) = paths{i}(end);
    ell(i) = length(paths{i});
end

x_i = 1./net.path_lengths(m_i);
x_f = 1./net.path_lengths(m_f);

y_i = clusters.mean_inc(m_i);
y_f = clusters.mean_inc(m_f);

ly_i = log10(y_i);

lr = x_f>x_i;
bt = y_f>y_i;

rl = ~lr;
tb = ~bt;

rl_bt = rl & bt;
lr_bt = lr & bt;
rl_tb = rl & tb;
lr_tb = lr & tb;

fprintf('Number of nodes in the directed network: %d; edges: %d; paths: %d\n',...
          sum(sum(adj,2)>0 | sum(adj,1)'>0), sum(sum(adj)), length(paths));
fprintf('Average direction:\n| %.2f | %.2f |\n| %.2f | %.2f |\n',...
          mean(lr_tb), mean(rl_tb), mean(lr_bt), mean(rl_bt) );
fprintf('Number in direction:\n| %d | %d |\n| %d | %d |\n',...
          sum(lr_tb), sum(rl_tb), sum(lr_bt), sum(rl_bt) );

ul = unique(ell);
mus = zeros(length(ul),1);
Ns = zeros(length(ul),1);
dirs = zeros(length(ul),4);
for i = 1:length(ul)
    id = ell == ul(i);
    mus(i) = mean(y_f(id)-y_i(id));
    Ns(i) = sum(id);
    dirs(i,:) = [mean(rl_bt(id)), mean(lr_bt(id)), mean(rl_tb(id)), mean(lr_tb(id))];
end

%% Set up grid and observations

if for_paper
  xcells = 300;
  ycells = 300;
else
  xcells = 100;
  ycells = 100;
end

% Inlclude starting points that don't have paths
x_i_all = 1./net.path_lengths;
ly_i_all = log10(clusters.mean_inc);

cmap = flipud(colormap('hot'));

[~,unconnected_ids] = setdiff([x_i_all,ly_i_all],[x_i,ly_i],'rows');

x_i_unconnected = x_i_all(unconnected_ids);
ly_i_unconnected = ly_i_all(unconnected_ids);

x_i_combined = [x_i;x_i_unconnected];
ly_i_combined = [ly_i;ly_i_unconnected];

[uniq_ss,~,uniq_idc] = unique([x_i_combined,ly_i_combined],'rows','stable');

gain = [y_f-y_i;zeros(length(unconnected_ids),1)];

gains = zeros(length(uniq_ss),1);
npaths = zeros(length(uniq_ss),1);

obs = [];
for i = 1:length(uniq_ss)
  gains(i) = mean(gain(uniq_idc == i));
  
  % Repeat observations of root node by number of paths
  npaths(i) = sum(uniq_idc == i);
  obs = [obs;repmat(uniq_ss(i,:),npaths(i),1)];
end

if for_paper
  [gridX,gridY] = meshgrid(linspace(0,0.037,xcells),...
                            linspace(0,log10(200),ycells));
else
  [gridX,gridY] = meshgrid(linspace(min(x_i_all),max(x_i_all)*1.05,xcells),...
                                  linspace(min(ly_i_all),max(ly_i_all)*1.05,ycells));
end
xi = [gridX(:),gridY(:)];

%% Plot the density 

[~,~,bw] = ksdensity(obs(:,1:2),xi);
% This gets a better match for the CDF
% bw = bw.*25;
bw = bw.*10;

density = ksdensity(obs(:,1:2),xi,'Bandwidth',bw);

grid_density = reshape(density,size(gridX));

figure( 'position', [383 506 849 414] );
imagesc([gridX(1),gridX(end)],...
               [gridY(1),gridY(end)],...
               grid_density);
hold on;

sz = npaths;
sz = (sz - min(sz)) ./ range(sz);
sz = sz.*50 + 1;
scatter(uniq_ss(:,1),uniq_ss(:,2),sz+1,'k',...
          'filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.8);
        
set(gca,'ydir','normal','fontsize',12,'box','on',...
        'ticklabelinterpreter','latex');
xlabel( 'Focal Node Centrality', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );
cmap = turbo(64); colormap(cmap); colorbar;

if for_paper
  set(gca,'ylim',[0 log10(200)],'xlim',[0 0.037],'ytick',[0 1 2]);
  yticks = get(gca,'YTickLabel');
  for i = 1:length(yticks)
    yticks{i} = sprintf('$10^{%s}$',yticks{i});
  end
  set(gca,'yticklabel',yticks);
  print( gcf, ['plots/density' appendix '.pdf'], '-dpdf', '-bestfit' );
end

%% Get expected values from kernel estimation

grid_exp_trunc = zeros(size(gridX));
grid_exp_pos = zeros(size(gridX));
grid_exp_neg = zeros(size(gridX));

gains_trunc = gains;
gains_trunc(gains_trunc<1) = 1;

pos_ids = gains>0;
neg_ids = gains<=0;

for i = 1:size(gridX,1)
  for j = 1:size(gridX,2)
    cX = ([gridX(i,j),gridY(i,j)] - uniq_ss) ./ bw;
    prob = mvnpdf(cX);
    grid_exp_trunc(i,j) = sum(prob .* gains_trunc);
    
    cX = ([gridX(i,j),gridY(i,j)] - uniq_ss) ./ (bw.*3.5);
    prob_pos = mvnpdf(cX);
    
    cX = ([gridX(i,j),gridY(i,j)] - uniq_ss) ./ (bw.*1.5);
    prob_neg = mvnpdf(cX);
    
    grid_exp_pos(i,j) = sum(prob_pos(pos_ids) .* gains(pos_ids));
    grid_exp_neg(i,j) = sum(prob_neg(neg_ids) .* gains(neg_ids));
  end
end

%% Plot the inferred expected value

figure( 'position', [383 506 849 414] );
imagesc([gridX(1),gridX(end)],...
               [gridY(1),gridY(end)],...
               grid_exp_trunc);
hold on;

sz = gains_trunc;
sz = (sz - min(sz)) ./ range(sz);
sz = sz.*50 + 1;
scatter(uniq_ss(:,1),uniq_ss(:,2),sz+1,'k',...
          'filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.8);

colormap(cmap);
colorbar;

if for_paper
  set(gca,'ytick',[0 1 2]);
  yticks = get(gca,'YTickLabel');
  for i = 1:length(yticks)
    yticks{i} = sprintf('$10^{%s}$',yticks{i});
  end
  set(gca,'YTickLabel',yticks);
end

set(gca,'ydir','normal','fontsize',12,'box','on',...
        'ticklabelinterpreter','latex');
xlabel( 'Focal Node Centrality', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );
if for_paper
  set(gca,'ylim',[0 log10(200)],'xlim',[0 0.037]);
  print( gcf, ['plots/heatmap' appendix '.pdf'], '-dpdf', '-bestfit' );
end

%% Plot the inferred conditional expectations (given virulent/benign)
divcmap = flipud(cbrewer('div','RdBu',127));

grid_exp_pos_norm = (grid_exp_pos - min(grid_exp_pos(:))) ./ range(grid_exp_pos(:));
grid_exp_neg_norm = (grid_exp_neg - min(grid_exp_neg(:))) ./ range(grid_exp_neg(:));

ind_exp_pos = ceil(grid_exp_pos_norm.*63)+1;
ind_exp_neg = ceil(grid_exp_neg_norm.*63)+1;

im_exp_pos = ind2rgb(ind_exp_pos,divcmap(64:end,:));
im_exp_neg = ind2rgb(ind_exp_neg,divcmap(1:64,:));

mv = max( abs(gains) );
          
figure( 'position', [383 506 849 414] );
impos = image([gridX(1),gridX(end)],...
               [gridY(1),gridY(end)],...
               im_exp_pos);

colormap(divcmap(64:end,:)); colorbar;
set(gca,'clim',[0 mv]);
set(gca,'ydir','normal','fontsize',12,...
          'box','on','ticklabelinterpreter','latex');
xlabel( 'Focal Node Centrality', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );
if for_paper
  set(gca,'ylim',[0 log10(200)],'xlim',[0 0.037],'ytick',[0 1 2]);
  yticks = get(gca,'YTickLabel');
  for i = 1:length(yticks)
    yticks{i} = sprintf('$10^{%s}$',yticks{i});
  end
  set(gca,'yticklabel',yticks);
end
        
hold on;

sz = gains(pos_ids);
sz = (sz - min(sz)) ./ mv;
sz = sz.*50 + 1;
shpos = scatter(uniq_ss(pos_ids,1),uniq_ss(pos_ids,2),round(sz+1),'r',...
                  'filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.8);

if for_paper
  print( gcf, ['plots/heatmap-pos' appendix '.pdf'], '-dpdf', '-bestfit' );
end
                
set(impos,'visible','off');
set(shpos,'visible','off');

imneg = image([gridX(1),gridX(end)],...
               [gridY(1),gridY(end)],...
               im_exp_neg);

sz = -gains(neg_ids);
sz = (sz - min(sz)) ./ mv;
sz = sz.*50 + 1;
scatter(uniq_ss(neg_ids,1),uniq_ss(neg_ids,2),sz+1,'b',...
          'filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.8);
        
colormap(divcmap(1:64,:));
set(gca,'clim',[-mv 0]);

if for_paper
  print( gcf, ['plots/heatmap-neg' appendix '.pdf'], '-dpdf', '-bestfit' );
end
        
colormap(divcmap); colorbar;
set(gca,'clim',[-mv mv]);

e10 = 1.2.^grid_exp_pos;
impos.AlphaData = (e10-min(e10(:))) ./ range(e10(:));
set(impos,'visible','on'); set(shpos,'visible','on');
e10 = 1.2.^(-grid_exp_neg);
imneg.AlphaData = (e10-min(e10(:))) ./ range(e10(:));

if for_paper
  print( gcf, ['plots/heatmap-posneg' appendix '.pdf'], '-dpdf', '-bestfit' );
end

figure;
h = cdfplot(gain);
set(h,'color','k','marker','.');
hold on;
h = cdfplot(grid_exp_trunc(:));
set(h,'color','r','linewidth',2);
lh = legend('Actual', 'Estimate','Location','Best');
xlabel('Change in prevalence','interpreter','latex');
ylabel('$F(x)$','interpreter','latex');
set(lh,'interpreter','latex','fontsize',12);
set(gca, 'fontsize', 12, 'ticklabelinterpreter', 'latex');
title('');
if for_paper
  print( gcf, ['plots/exp-trunc-cdf' appendix '.pdf'], '-dpdf', '-bestfit' );
end

figure;
h = cdfplot(gain(gain>0));
set(h,'color','k','marker','.');
hold on;
h = cdfplot(grid_exp_pos(:));
set(h,'color','r','linewidth',2);
lh = legend('Actual', 'Estimate','Location','Best');
xlabel('Change in prevalence','interpreter','latex');
ylabel('$F(x)$','interpreter','latex');
set(lh,'interpreter','latex','fontsize',12);
set(gca, 'fontsize', 12, 'ticklabelinterpreter', 'latex');
title('');
if for_paper
  print( gcf, ['plots/exp-pos-cdf' appendix '.pdf'], '-dpdf', '-bestfit' );
end

figure;
h = cdfplot(gain(gain<=0));
set(h,'color','k','marker','.');
hold on;
h = cdfplot(grid_exp_neg(:));
set(h,'color','r','linewidth',2);
lh = legend('Actual', 'Estimate','Location','Best');
xlabel('Change in prevalence','interpreter','latex');
ylabel('$F(x)$','interpreter','latex');
set(lh,'interpreter','latex','fontsize',12);
set(gca, 'fontsize', 12, 'ticklabelinterpreter', 'latex');
title('');
if for_paper
  print( gcf, ['plots/exp-neg-cdf' appendix '.pdf'], '-dpdf', '-bestfit' );
end

%% Path characteristics

figure( 'position',[1169 618 570 341] );
hold on;
linstyle = {'ko-', 'k^-','k.-','k--'};
for i = 1:4
  plot(ul,dirs(:,i),linstyle{i});
end
xlabel('Path length', 'interpreter', 'latex');
ylabel('Proportion of paths', 'interpreter', 'latex');
lh = legend( sprintf('RL-BT (%.2f)', mean(rl_bt)),...
              sprintf('LR-BT (%.2f)', mean(lr_bt)),...
              sprintf('RL-TB (%.2f)', mean(rl_tb)),...
              sprintf('LR-TB (%.2f)', mean(lr_tb)),...
              'location','NorthWest');

set(lh,'interpreter','latex','fontsize',12);
set(gca, 'fontsize', 12, 'ticklabelinterpreter', 'latex');
axis tight;

if for_paper
  print(gcf, ['plots/length-direction' appendix '.pdf'], '-dpdf');
end

figure( 'position',[1169 618 570 341] );
hold on;
linstyle = {'ko-', 'k^-','k.-','k--'};
for i = 1:4
  plot(ul,dirs(:,i).*Ns,linstyle{i});
end
xlabel('Path length', 'interpreter', 'latex');
ylabel('Number of paths', 'interpreter', 'latex');
lh = legend( sprintf('RL-BT (%d)', round(mean(rl_bt).*sum(Ns))),...
            sprintf('LR-BT (%d)', round(mean(lr_bt).*sum(Ns))),...
            sprintf('RL-TB (%d)', round(mean(rl_tb).*sum(Ns))),...
            sprintf('LR-TB (%d)', round(mean(lr_tb).*sum(Ns))),...
            'location','NorthWest');

set(lh,'interpreter','latex','fontsize',12);
set(gca, 'fontsize', 12, 'ticklabelinterpreter', 'latex');
axis tight;

if for_paper
  print(gcf, ['plots/length-direction-num' appendix '.pdf'], '-dpdf');
end

cgain = y_f-y_i;

[rho,pval] = corr(ell,cgain);
figure;
plot(ell,cgain,'k.');
hold on;
c = polyfit(ell,cgain,1);
y_est = polyval(c,ell);
plot(ell,y_est,'r-');
title('Path length to delta prevalence');
xlabel('Path length','interpreter','latex');
ylabel('Change in prevalence','interpreter','latex');
lh = legend('Data',sprintf('r = %.3g',rho),'Location','Best');

set(lh,'interpreter','latex');
set(gca,'fontsize',12,'box','on','ticklabelinterpreter','latex');
title('');

fprintf('Correlation of path length to prevalence: %.3g [M = %d, p = %.3g]\n',rho,length(ell),pval);

%% Path characteristics within basin

% Define basin as 75% of the max expectation
basin_cutoff = 0.75.*max(grid_exp_trunc(:));

within_basin = grid_exp_trunc >= basin_cutoff;

min_xi = min(gridX(within_basin));
max_xi = max(gridX(within_basin));

min_lyi = min(gridY(within_basin));
max_lyi = max(gridY(within_basin));

fprintf('Expected value cut-off for "basin": %.3f\n', basin_cutoff);
fprintf('Basin {centrality range, prevalence range (log)}: {%.3g--%.3g, %.3g--%.3g (%.3g--%.3g)}\n',...
            min_xi, max_xi, 10.^min_lyi, 10.^max_lyi, min_lyi, max_lyi);

paths_within_basin = x_i > min_xi & x_i < max_xi & ly_i > min_lyi & ly_i < max_lyi;

fprintf('Proportion of paths within basin: %.3f\n', mean(paths_within_basin))

x_ib = x_i(paths_within_basin);
x_fb = x_f(paths_within_basin);

y_ib = y_i(paths_within_basin);
y_fb = y_f(paths_within_basin);

unodes = unique([x_ib,y_ib],'rows');
unodes_orig = unique([x_i,y_i],'rows');

fprintf('Proportion of nodes within basin: %.3g\n', length(unodes) ./ length(uniq_ss));
fprintf('Proportion of nodes with >0 paths within basin: %.3g\n', length(unodes) ./ length(unodes_orig));

lrb = x_fb>x_ib;
btb = y_fb>y_ib;

rlb = ~lrb;
tbb = ~btb;

rl_btb = rlb & btb;
lr_btb = lrb & btb;
rl_tbb = rlb & tbb;
lr_tbb= lrb & tbb;

fprintf('[In basin] Average direction:\n| %.2f | %.2f |\n| %.2f | %.2f |\n',...
          mean(lr_tbb), mean(rl_tbb), mean(lr_btb), mean(rl_btb) );
fprintf('[In basin] Number in direction:\n| %d | %d |\n| %d | %d |\n',...
          sum(lr_tbb), sum(rl_tbb), sum(lr_btb), sum(rl_btb) );

ellb = ell(paths_within_basin);
        
ulb = unique(ell);
musb = zeros(length(ulb),1);
Nsb = zeros(length(ulb),1);
dirsb = zeros(length(ulb),4);
for i = 1:length(ulb)
    id = ellb == ulb(i);
    musb(i) = mean(y_fb(id)-y_i(id));
    Nsb(i) = sum(id);
    dirsb(i,:) = [mean(rl_btb(id)), mean(lr_btb(id)), mean(rl_tbb(id)), mean(lr_tbb(id))];
end

[rho,pval] = corr(ellb,y_fb-y_ib);

fprintf('[In basin] Correlation of path length to prevalence: %.3g [N = %d, p = %.3g]\n',rho,length(ellb),pval);
