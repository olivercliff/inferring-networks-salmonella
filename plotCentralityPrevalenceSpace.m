function plotCentralityPrevalenceSpace(matfile,approach,multiple_loci,marker_scaler,plot_case_studies)

if multiple_loci
  fld = 'multiloci';
else
  fld = 'singleloci';
end

if nargin < 5
  plot_case_studies = false;
end

appendix = [ '_' approach '_' fld];

load(matfile,'network');

net = network.(approach);
clusters = net.clusters.overlapping;
adj = net.(fld).adj_mat;

%% Plot chains on genetic graph
marker_alpha = .7;


% Get the `superbug' chains
in_graph = sum(adj,2)>0 | sum(adj,1)'>0;

% Get unique MLVA: find(all(ismember(mlvas,[3,12,9,10,550]),2))
% Get dates: datestr(dns(is(find(js==10,1,'first'))))
mp = max(clusters.mean_inc(in_graph));
mpn = find(clusters.mean_inc == mp);

fprintf('Most prevalent node: %s (prevalence: %.3g)\n',...
          mat2str(net.profile(mpn,:)), clusters.mean_inc(mpn));

% Plot which nodes are in the graph
try
    G = digraph(adj,net.names);
catch
    G = digraph(adj);
end
components = conncomp(G,'Type','weak')';
unique_comp = unique(components);
size_comp = zeros(size(components));
% sum_comp_inc = zeros(size(components));
for c = 1:length(unique_comp)
    ids = components == unique_comp(c);
    size_comp(ids) = sum(ids);
%     sum_comp_inc(ids) = sum(mean_inc(ids).*num_members(ids));
end

fprintf('Number of edges: %d\n', height(G.Edges));
fprintf('Number of nodes: %d\n', sum(in_graph));

load('utils/yellow-green_cmap.mat','cmap');

distance_mode = true;

if nargin < 4
    zoom_modes = 0:2; % 0 = no zoom, 1 == north ryde, 2 == cabra/turra/longest
else
    zoom_modes = 3;
end

node_sz_range = 15;
node_sz_min = 5;

for i = 1:length(zoom_modes)
  
  figure( 'position', [383 506 849 414] );
  
  zoom_mode = zoom_modes(i);

  if nargin < 4
      marker_scaler = 1;
      if zoom_mode == 1
        marker_scaler = 10;
      elseif zoom_mode == 2
        marker_scaler = 5;
      end
  end

  node_sz = (clusters.num_members - min(clusters.num_members)) ./ range(clusters.num_members);
  node_sz = node_sz .* node_sz_range + node_sz_min;
  
  node_sz(isnan(node_sz)) = 1;

  hold on;
  scatter( 1 ./ net.path_lengths(~in_graph), clusters.mean_inc(~in_graph),...
            node_sz(~in_graph) .* marker_scaler, [.7 .7 .7],...
            'filled', 'markerfacealpha', marker_alpha );
  if distance_mode
    scatter( 1 ./ net.path_lengths(in_graph), clusters.mean_inc(in_graph),...
              node_sz(in_graph) .* marker_scaler, net.pdist( in_graph, mpn ),...
              'filled', 'markerfacealpha', marker_alpha );
    colormap(cmap);
  else
    scatter( 1 ./ net.path_lengths(in_graph), clusters.mean_inc(in_graph),...
              node_sz(in_graph) .* marker_scaler, log10( size_comp( in_graph ) ),...
              'filled', 'markerfacealpha', marker_alpha );
    colormap(flipud(cmap));
    caxis([0.3010 2.4409]);
  end

  if plot_case_studies
    
    paths = net.(fld).unique_paths;
    
    case_studies = [3,10,7,14,523;
                    3,12,9,10,550;
                    3,16,9,11,523];

    cs_col = {[1 0 1], [0 0 1], [1 0 0], [0 0 0]};
    
    [~,cs_ids] = intersect(network.unadjusted.profile,case_studies,'rows');
      
    % Find the longest paths that starts with these MLVAs
    [~,cs_paths] = findLongestPath(paths,cs_ids);
    cs_ls = {'-','-','-',':'};
    
    [~,mpid] = max(cellfun(@length,paths))
    cs_paths = [cs_paths;paths(mpid)];

    for s = 1:length(cs_paths)
      cs_xs = 1./net.path_lengths(cs_paths{s});
      cs_ys = clusters.mean_inc(cs_paths{s});

      mk_sz = 3;
      plot( cs_xs, cs_ys, cs_ls{s}, 'color', cs_col{s}, 'linewidth', 1 );
      plot( cs_xs(1), cs_ys(1), 'o', 'color', cs_col{s}, 'markersize', mk_sz .* marker_scaler );
    end
  end

  set( gca, 'fontsize', 12, 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
  xlabel( 'Focal Node Centrality', 'interpreter', 'latex' );
  ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );

  % axis tight;
  apx2 = '';
  if distance_mode
    apx2 = '-dist';
  end

  if zoom_mode < 3
    set(gca,'DataAspectRatioMode', 'manual', 'DataAspectRatio',[1 8000 50]);
  end
  if zoom_mode == 2
    set(gca,'xlim',[0.0297 0.0365],'ylim',[15.0807 69.5206]);
  elseif zoom_mode == 1
    set(gca,'xlim',[0.01855 0.019],'ylim',[18.6089 22.5968]);
  elseif zoom_mode == 0
    set(gca,'ylim',[0 200],'xlim',[0 0.037]);
  end

  if zoom_mode == 0
    print( gcf, ['plots/network-trace' appendix apx2 '.pdf'], '-dpdf', '-bestfit' );
  elseif zoom_mode == 1
    print( gcf, ['plots/network-trace' appendix apx2 '_north-ryde.pdf'], '-dpdf', '-bestfit' );
  elseif zoom_mode == 2
    print( gcf, ['plots/network-trace' appendix apx2 '_cabra-turra.pdf'], '-dpdf', '-bestfit' );
  elseif zoom_mode == 3
    colorbar;
    print( gcf, ['plots/network-trace' appendix apx2 sprintf('-D%d',network.cluster_D) '.pdf'], '-dpdf', '-bestfit' );
  end
end

%% Plot correlations

figure('position',[257   660   570   357]);
plot( net.pdist( in_graph, mpn ), clusters.mean_inc(in_graph), 'k.' );

set( gca, 'fontsize', 12, 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
xlabel( 'Distance to Most Prevalent Node', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );

% curve = fit( D( in_graph, mp_node ), log10(mean_inc( in_graph )), 'poly1' );
% hold on;
% min_X = min(D( in_graph, mp_node ));
% max_X = max(D( in_graph, mp_node ));
% plot( [min_X, max_X], [curve(min_X), curve(round(max_X/2))], 'r-' );

[rho,pval] = corr(net.pdist(in_graph,mpn),log10(clusters.mean_inc(in_graph)));
fprintf('Correlation of distance-to-MPN to log-prevalence: %.3g [N = %d, p=%.3g]\n', rho, sum(in_graph), pval);


print( gcf, ['plots/dist_prev' appendix apx2 '.pdf'], '-dpdf', '-bestfit' );

figure('position',[832   659   570   358]);
plot( size_comp(in_graph), clusters.mean_inc(in_graph), 'k.' );

set( gca, 'fontsize', 12, 'xscale', 'log', 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
xlabel( 'Size of Component', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );

[rho,pval] = corr(log10(size_comp(in_graph)),log10(clusters.mean_inc(in_graph)));
fprintf('Correlation of log-component-size to log-prevalence: %.3g [N = %d, p=%.3g]\n', rho, sum(in_graph), pval);
