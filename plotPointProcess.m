function plotPointProcess(matfile,approach,multiple_loci,plot_case_studies)

  if nargin < 4
    plot_case_studies = false;
  end

  if multiple_loci
    fld = 'multiloci';
  else
    fld = 'singleloci';
  end

  appendix = [ '_' approach '_' fld];
  
  addpath('utils');

  load(matfile,'network','events');

  st_dn = min(events.dates);
  
  M = size(network.unadjusted.profile,1);
  N = days(max(events.dates) - st_dn);
  
  point_process = zeros(M,N);
  first_appearance = zeros(M,1);

  for m = 1:M
    cdates = events.dates(events.(approach).occurances(:,m))';

    ids = days(cdates-st_dn)+1;
    point_process(m,ids) = 1;
    first_appearance(m) = ids(1);
  end

  [~,sort_ids] = sort(first_appearance,'descend');


  [ys, xs] = find(point_process(sort_ids,:));
  figure('position',[302   651   418   366]);
  plot(xs+st_dn-1, ys, 'k.','markersize',1)

  if plot_case_studies
    
    if ~strcmp(approach,'unadjusted')
      fprintf('Cannot plot case studies unless using ''unadjusted'' method.');
      return;
    end
    
    net = network.(approach);
    paths = net.(fld).unique_paths;

    hold on;

    case_studies = [3,10,7,14,523;
                    3,12,9,10,550;
                    3,16,9,11,523;];

    cs_col = {[1 0 1], [0 0 1], [1 0 0]};
    
    [~,cs_ids] = intersect(network.unadjusted.profile,case_studies,'rows');
      
    % Find the longest paths that starts with these MLVAs
    cs_path_ids = findLongestPath(paths,cs_ids);
    
    for s = 1:size(case_studies,1)
      cs_path = paths{cs_path_ids(s)};
      
      for m = 1:length(cs_path)
        id = cs_path(m);
        csb_xs = find(point_process(id ,:));
        csb_ys = find(sort_ids == id );
        sh = scatter(csb_xs+st_dn-1, csb_ys.*ones(size(csb_xs)), 5, cs_col{s},'filled');
        sh.MarkerFaceAlpha = 0.5;
      end
    end
  end
  axis tight;
  set(gca,'ytick',[],'TickLabelInterpreter','latex','fontsize',12,'box','off');
  ylabel('Unique MLVA Profiles','Interpreter','latex');
  xlabel('Date','Interpreter','latex');
  datetick('x','yyyy');

  print( gcf, ['plots/point-process' appendix '.pdf'], '-dpdf', '-bestfit' );

%   if plot_case_studies
%   
%     figure('position',[681   711   840   359]);
%     col_order = get(gca,'colororder');
% 
%     [sb_ys, sb_xs] = find(cs_pp);
%     hold on;
%     for m = 2:length(cs)
%       [~, cxs] = find(cs_pp(m,:));
% 
%       dist_to_prev = D(cs(m),cs(m-1));
%       cpre = pre_windows(cs(m),dist_to_prev);
% 
%       cpost = post_windows(cs(m),dist_to_prev);
% 
%       fill_x = [cxs(1), cxs(1)-cpre, cxs(1)+cpost]+st_dn-1;
%       fill_y = [m, m-1, m-1]-.1;
%       fh = fill( fill_x, fill_y, col_order(dist_to_prev,:));
%       set(fh,'facealpha',0.25,'edgecolor','none');
%     end
% 
%     sh = scatter(sb_xs+st_dn-1, sb_ys, 50, 'k', 'filled');
%     sh.MarkerFaceAlpha = 0.3;
%   
%     axis tight;
% 
%     set(gca,'ylim',[0 6],'ydir','reverse','ytick',1:5,'YTickLabel', mlva_names,'TickLabelInterpreter','latex','fontsize',12);
% 
%     datetick('x','yyyy');
% 
%     print( gcf, ['algorithm' appendix '.pdf'], '-dpdf', '-bestfit' );
%   end
end
