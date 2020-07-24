function plotCDFs(matfile,approach,multiple_loci)

if multiple_loci
  fld = 'multiloci';
else
  fld = 'singleloci';
end

load( matfile, 'network' );
cdfs = network.(approach).(fld).cdfs;

D = size(cdfs.cdm_ptg,2);

pdf_lab = cell(D,1);
cdf_lab = cell(D,1);
for d = 1:D
  pdf_lab{d} = sprintf('G=%d',d);
  cdf_lab{d} = sprintf('G=%d',d);
end

% Plot time-lag CDF

maxD_to_plot = 4;

figure('position',[527   449   567   387]);
plot( cdfs.temporal_edges(2:end-1), cdfs.cdm_ptg(2:end,1:maxD_to_plot) );
lh = legend(cdf_lab{1:maxD_to_plot},'location','NorthWest');
set(lh,'interpreter','latex','box','off');
set(gca, 'yscale', 'log', 'xlim', [0 60], 'ylim', [1e-5 1e-1],...
    'fontsize',12,'ticklabelinterpreter', 'latex');
grid on;
xlabel('Days before first appearance','interpreter','latex');
ylabel('Probability','interpreter','latex');

print(gcf, [ 'plots/cdm_ptg_' approach '_' fld '.pdf'], '-dpdf' );

% Plot time-lead CDF

figure('position',[1102 416 342 415]);
plot( cdfs.temporal_edges(2:end-1), cdfs.cdm_ftg(2:end,1:maxD_to_plot) );
set(gca, 'yscale', 'log', 'xlim', [0 30], 'ylim', [1e-5 1e-1],...
  'fontsize',12,'ticklabelinterpreter', 'latex');
grid on;
xlabel('Days after first appearance','interpreter','latex');
ylabel('Probability','interpreter','latex');

print(gcf, [ 'plots/cdm_ftg_' approach '_' fld '.pdf'], '-dpdf' );