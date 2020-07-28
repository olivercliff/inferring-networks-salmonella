function plotPathVsPrevalence(matfile,approach,multiple_loci)

addpath ./utils/

plot_geographic = false;
plot_case_studies = false;

if multiple_loci
  fld = 'multiloci';
else
  fld = 'singleloci';
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

lr = x_f>x_i;
bt = y_f>y_i;

rl = ~lr;
tb = ~bt;

rl_bt = rl & bt;
lr_bt = lr & bt;
rl_tb = rl & tb;
lr_tb= lr & tb;

fprintf('Number of nodes in the directed network: %d; number of paths: %d\n',...
          sum(sum(adj,2)>0), length(paths));
fprintf('Average direction:\n| %.2f | %.2f |\n| %.2f | %.2f |\n',...
          mean(lr_tb), mean(rl_tb), mean(lr_bt), mean(rl_bt) );

ul = unique(ell);
mus = zeros(length(ul),1);
Ns = zeros(length(ul),1);
dirs = zeros(length(ul),4);
i_coord = zeros(length(ul),2);
for i = 1:length(ul)
    id = ell == ul(i);
    mus(i) = mean(y_f(id)-y_i(id));
    Ns(i) = sum(id);
    dirs(i,:) = [mean(rl_bt(id)), mean(lr_bt(id)), mean(rl_tb(id)), mean(lr_tb(id))];
end

%% Contour state-space map (expected gain in prevalence)
xcells = 10; ycells = 10;

gain = y_f-y_i;

ly_i = log10(y_i);

[idX,edgeX] = discretize(x_i,xcells);
[idY,edgeY] = discretize(ly_i,ycells);

gridGain = zeros(length(edgeX), length(edgeY));
gridNum = zeros(length(edgeX), length(edgeY));
for i = 1:length(edgeY)
  for j = 1:length(edgeX)
    gridGain(i,j) = mean(gain(idX == j & idY == i));
    gridNum(i,j) = sum(idX == j & idY == i);
  end
end

% Unique starting points, do not repeat by prevalence
% heatmap_opt = 1;
% Unique starting points, repeat by prevalence
heatmap_opt = 2;

cmap = flipud(colormap('hot'));
% cmap = turbo(64);

[gridX,gridY] = meshgrid(edgeX,edgeY);

[uniq_ss,~,uniq_idc] = unique([x_i,ly_i],'rows');

if heatmap_opt == 2
    obs = [];
    inv_obs = [];
    reps = [];
    inv_reps = [];
end

gains = zeros(length(uniq_ss),1);
for i = 1:length(uniq_ss)
  gains(i) = max(gain(uniq_idc == i));
  
  if heatmap_opt == 2
    if gains(i) > 0
      reps = [reps; ceil(gains(i))];
      obs = [obs;repmat(uniq_ss(i,:),reps(end),1)];
    else
      inv_reps = [inv_reps; ceil(-gains(i))+1];
      inv_obs = [inv_obs;repmat(uniq_ss(i,:),inv_reps(end),1)];
    end   
  end
end

if heatmap_opt == 1
  obs = uniq_ss(gains > 0, :);
  inv_obs = uniq_ss(gains <=0, :);
end

[gridXfine,gridYfine] = meshgrid(linspace(edgeX(1),edgeX(end)*1.05,100),...
                                  linspace(edgeY(1),edgeY(end)*1.05,100));
xi = [gridXfine(:),gridYfine(:)];

f = ksdensity(obs(:,1:2),xi);
gridF = reshape(f,size(gridXfine));

figure;
im = imagesc([gridXfine(1),gridXfine(end)],...
               [gridYfine(1),gridYfine(end)],...
                gridF);
im.AlphaData = 0.8;

colormap(cmap);
hold on;
if heatmap_opt == 1
  plot(obs(:,1),obs(:,2),'k.');
else
  uniq_obs = uniq_ss(gains>0,:);
  scatter(uniq_obs(:,1),uniq_obs(:,2),reps.*2,'k',...
            'filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.8);
end
axis tight;
xlabel('Centrality');
ylabel('Prevalence (log)');
zlabel('Probability of increased prevalence');
set(gca,'ydir','normal');
title('Probability of virulent strain');

inv_f = ksdensity(inv_obs(:,1:2),xi);
gridinvF = reshape(inv_f,size(gridXfine));

figure;
im = imagesc([gridXfine(1),gridXfine(end)],...
              [gridYfine(1),gridYfine(end)],...
              gridinvF);
im.AlphaData = 0.8;

colormap(cmap);
hold on;
if heatmap_opt == 1
  plot(inv_obs(:,1),inv_obs(:,2),'k.');
else
  uniq_inv_obs = uniq_ss(gains<=0,:);
  scatter(uniq_inv_obs(:,1),uniq_inv_obs(:,2),inv_reps.*2,'k',...
            'filled','MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.8);
end
axis tight;
xlabel('Centrality');
ylabel('Prevalence (log)');
zlabel('Probability of increased prevalence');
set(gca,'ydir','normal')
title('Probability of benign strain');

figure;
h = pcolor(gridX, gridY, gridGain);
set(h,'edgecolor',[.7 .7 .7],'facealpha',.3);
hold on;

sz = abs(gains).*2;
scatter(uniq_ss(:,1),uniq_ss(:,2),sz+1,gains,'filled',...
          'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'MarkerEdgeColor','k');

set(gca,'clim',[-15 15],'ydir','normal');

cmap = redblue(64);
colormap(cmap);

xlabel('Centrality');
ylabel('Prevalence (log)');
        
% exploitative = gains > 0;
% Mdl = fitcsvm(uniq_ss, exploitative,'Weights',abs(gains)./max(abs(gains)),...
%                     'KernelFunction','rbf',...
%                     'Standardize',true,...
%                     'OptimizeHyperparameters','auto',...
%                     'KernelScale','auto');
% 
% exploitativehat = predict(Mdl,uniq_ss);
% 
% C = confusionmat(exploitative,exploitativehat);
% fprintf('TPR: %.3f\n', C(2,2) / (C(1,2) + C(2,2)));
% fprintf('FPR: %.3f\n', C(1,2) / (C(1,2) + C(2,2)));
% 
% figure;
% confusionchart(exploitative,exploitativehat);
% 
% % figure;
% sv = Mdl.SupportVectors;
% plot(sv(:,1),sv(:,2),'ko');

%% Path characteristics

figure( 'position',[1169 618 570 341] );
hold on;
linstyle = {'ko-', 'k^-','k.-','k--'};
for i = 1:4
  plot(ul,dirs(:,i),linstyle{i});
end
xlabel('Chain length', 'interpreter', 'latex');
ylabel('Proportion of chains', 'interpreter', 'latex');
lh = legend( sprintf('RL-BT (%.2f)', mean(rl_bt)),...
              sprintf('LR-BT (%.2f)', mean(lr_bt)),...
              sprintf('RL-TB (%.2f)', mean(rl_tb)),...
              sprintf('LR-TB (%.2f)', mean(lr_tb)),...
              'location','NorthWest');

set(lh,'interpreter','latex','fontsize',12);
set(gca, 'fontsize', 12, 'ticklabelinterpreter', 'latex');
axis tight;
print(gcf, ['plots/length-direction' appendix '.pdf'], '-dpdf');

figure( 'position',[1169 618 570 341] );
hold on;
linstyle = {'ko-', 'k^-','k.-','k--'};
for i = 1:4
  plot(ul,dirs(:,i).*Ns,linstyle{i});
end
xlabel('Chain length', 'interpreter', 'latex');
ylabel('Number of chains', 'interpreter', 'latex');
lh = legend( sprintf('RL-BT (%d)', round(mean(rl_bt).*sum(Ns))),...
            sprintf('LR-BT (%d)', round(mean(lr_bt).*sum(Ns))),...
            sprintf('RL-TB (%d)', round(mean(rl_tb).*sum(Ns))),...
            sprintf('LR-TB (%d)', round(mean(lr_tb).*sum(Ns))),...
            'location','NorthWest');

set(lh,'interpreter','latex','fontsize',12);
set(gca, 'fontsize', 12, 'ticklabelinterpreter', 'latex');
axis tight;
print(gcf, ['plots/length-direction-num' appendix '.pdf'], '-dpdf');

[rho,pval] = corr(ell,y_f-y_i);
figure;
plot(ell,y_f-y_i,'k.');
hold on;
plot(ul,mus,'k--');

fprintf('Correlation of path length to incidence: %.2f [p = %.3g]\n', rho,pval);
