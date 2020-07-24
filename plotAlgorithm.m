close all
clear

max_distance = 4;

% Which processing method are we using?
% appendix1 = '-27-a-33'; % 5th loci has 27 and 33bp summed
% appendix1 = '-27-p-33'; % 5th loci has 27 and 33bp summed
appendix1 = '';

% appendix2 = '-UBL';
appendix2 = '-ML';
% appendix2 = '-SL';

appendix3 = sprintf('-D%i', max_distance);

appendix = [ appendix1 appendix2 appendix3 ];

plot_geographic = false;

%% Load and preprocess data
load( [ 'overlapping-clusters' appendix1 '.mat' ] );
load( [ 'windows' appendix1 appendix2 '.mat' ] );
load( [ 'network' appendix '.mat' ] );
load( [ 'distance-windows' appendix '.mat' ] );
load( [ 'centralities' appendix1 '.mat' ] );
load( [ 'superbug-chains' appendix1 appendix2 '.mat'] );

M = size(mlvas_adjusted,1);
N = 3287; % Number of days between 1/1/09 and 31/12/16

st_dn = datenum('01/01/2008','dd/mm/YYYY');

point_process = zeros(M,N);
first_appearance = zeros(M,1);

for i = 1:M
  cdates = mlva_dates{i};
  
  point_process(i,cdates-st_dn+1) = 1;
  first_appearance(i) = cdates(1);
end

[~,sort_ids] = sort(first_appearance,'descend');

sb_chain = sb_chains{4};

[ys, xs] = find(point_process(sort_ids,:));
figure('position',[302   651   418   366]);
plot(xs+st_dn-1, ys, 'k.','markersize',1)
hold on;

mlva_names = cell(length(sb_chain),1);

sb_col = {[0 0 1], [1 0 1], [1 0 0]};

sb_pp = point_process(sb_chain,:);
for s = 2:4
  csb_chain = sb_chains{s};
  for i = 1:length(csb_chain)
    csb_xs = find(point_process(csb_chain(i),:));
    csb_ys = find(sort_ids == csb_chain(i));
    sh = scatter(csb_xs+st_dn-1, csb_ys.*ones(size(csb_xs)), 5, sb_col{s-1},'filled');
    sh.MarkerFaceAlpha = 0.5;

    if s == 4
      mlva_names{i} = mat2str(mlvas_adjusted( csb_chain(i), : ) );
      mlva_names{i} = strrep( mlva_names{i}(2:end-1), ' ', '-' );
    end
  end
end
axis tight;
set(gca,'ytick',[],'TickLabelInterpreter','latex','fontsize',12,'box','off');
ylabel('Unique MLVAs','Interpreter','latex');
xlabel('Date','Interpreter','latex');
datetick('x','yyyy');

print( gcf, ['point-process' appendix '.pdf'], '-dpdf', '-bestfit' );

figure('position',[681   711   840   359]);
col_order = get(gca,'colororder');
[sb_ys, sb_xs] = find(sb_pp);

hold on;
for i = 2:length(sb_chain)
  [cys, cxs] = find(sb_pp(i,:));
  
  dist_to_prev = D(sb_chain(i),sb_chain(i-1));
  cpre = pre_windows(sb_chain(i),dist_to_prev);
  
  cpost = post_windows(sb_chain(i),dist_to_prev);
  
  fill_x = [cxs(1), cxs(1)-cpre, cxs(1)+cpost]+st_dn-1;
  fill_y = [i, i-1, i-1]-.1;
  fh = fill( fill_x, fill_y, col_order(dist_to_prev,:));
  set(fh,'facealpha',0.25,'edgecolor','none');
  
%   plot( fill_x(1:2), fill_y(1:2), '-',...
%         'color', col_order(dist_to_prev,:));
%   plot( fill_x([1 3]), fill_y([1 3]), '-',...
%         'color', col_order(dist_to_prev,:));
end

sh = scatter(sb_xs+st_dn-1, sb_ys, 50, 'k', 'filled');
sh.MarkerFaceAlpha = 0.3;
axis tight;

set(gca,'ylim',[0 6],'ydir','reverse','ytick',1:5,'YTickLabel', mlva_names,'TickLabelInterpreter','latex','fontsize',12);

datetick('x','yyyy');

print( gcf, ['algorithm' appendix '.pdf'], '-dpdf', '-bestfit' );