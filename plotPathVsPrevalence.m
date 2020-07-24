function plotPathVsPrevalence(matfile,approach,multiple_loci)

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

figure;
plot(x_i,y_f - y_i,'k.');

figure;
plot(y_i,y_f - y_i,'k.');

figure;
gain = y_f-y_i;
exploit_ids = gain > 0;
% plot([x_i(exploit_ids) x_f(exploit_ids)], [y_i(exploit_ids) y_f(exploit_ids)],'color',[.8 .8 .8]);
hold on;
scatter(x_i(~exploit_ids),y_i(~exploit_ids),20,[.4 .4 .4],'filled');
scatter(x_i(exploit_ids),y_i(exploit_ids),20,gain(exploit_ids));
scatter(x_f(exploit_ids),y_f(exploit_ids),20,gain(exploit_ids),'filled');
set(gca,'yscale','log');
colormap('copper');

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
