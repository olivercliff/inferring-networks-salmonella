function plotChains(matfile,approach,multiple_loci)

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

%% Load and preprocess data

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

%% Plot chains on genetic graph
marker_alpha = .7;


% Get the `superbug' chains
in_graph = sum(adj,2)>0 | sum(adj,1)'>0;

% Get unique MLVA: find(all(ismember(mlvas,[3,12,9,10,550]),2))
% Get dates: datestr(dns(is(find(js==10,1,'first'))))

[~,mp_node] = max(clusters.mean_inc);
mp = max(clusters.mean_inc(in_graph));
mp_node = find(clusters.mean_inc == mp);


% Plot which nodes are in the graph
G = digraph(adj);
components = conncomp(G,'Type','weak')';
unique_comp = unique(components);
size_comp = zeros(size(components));
% sum_comp_inc = zeros(size(components));
for c = 1:length(unique_comp)
    ids = components == unique_comp(c);
    size_comp(ids) = sum(ids);
%     sum_comp_inc(ids) = sum(mean_inc(ids).*num_members(ids));
end

load('utils/yellow-green_cmap.mat');

% figure( 'position', [460 369 1110 564] );
figure( 'position', [383 506 849 414] );

distance_mode = true;
zoom_mode = 1; % 0 = no zoom, 1 == north ryde, 2 == cabra/turra/longest

node_sz_range = 15;
node_sz_min = 5;

marker_scaler = 1;
if zoom_mode == 1
  marker_scaler = 10;
elseif zoom_mode == 2
  marker_scaler = 5;
end

node_sz = (clusters.num_members - min(clusters.num_members)) ./ range(clusters.num_members);
node_sz = node_sz .* node_sz_range + node_sz_min;

hold on;
scatter( 1 ./ net.path_lengths(~in_graph), clusters.mean_inc(~in_graph),...
          node_sz(~in_graph) .* marker_scaler, [.7 .7 .7],...
          'filled', 'markerfacealpha', marker_alpha );
if distance_mode
  scatter( 1 ./ net.path_lengths(in_graph), clusters.mean_inc(in_graph),...
            node_sz(in_graph) .* marker_scaler, net.pdist( in_graph, mp_node ),...
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
  north_ryde = 346; % Vitali north ryde, isolated chain (length 9)
  rockdale = 754; % Rockdale, 2014
  turramurra = 158; % Turramurra, 2013
  longest_chain = 354; % Longest chain
  cabra = 691; % Cabra, 2015

  % sb_nodes = longest_chain;
  % sb_nodes = turramurra;
  sb_nodes = [longest_chain north_ryde turramurra cabra];

  sb_chains = cell(length(sb_nodes),1);

  sb_col = {[0 0 0],[0 0 1], [1 0 1], [1 0 0]};
  sb_ls = {':','-','-','-'};

  % sb_col = {[1 0 1]};
  % sb_ls = {'-'};
  for s = 1:length(sb_nodes)
    in_paths = false(U,1);
    for i = 1:U
      if sum(ismember(paths{i},sb_nodes(s))) == 1
        in_paths(i) = true;
      end
    end

    superbug_chains = find(in_paths);
    [~,max_sb_chain] = max(ell(superbug_chains));

    sb_chain = paths{superbug_chains(max_sb_chain)};

    sb_xs = 1./net.path_lengths(sb_chain);
    sb_ys = clusters.mean_inc(sb_chain);

    mk_sz = 3;
    plot( sb_xs, sb_ys, sb_ls{s}, 'color', sb_col{s}, 'linewidth', 1 );
    plot( sb_xs(1), sb_ys(1), 'o', 'color', sb_col{s}, 'markersize', mk_sz .* marker_scaler );

    sb_chains{s} = sb_chain;
  end
end

set( gca, 'fontsize', 12, 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
% xlabel( 'Centrality', 'interpreter', 'latex' );
xlabel( 'Focal Node Centrality', 'interpreter', 'latex' );
% ylabel( 'Prevalence', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );

set( gca, 'fontsize', 12, 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
xlabel( 'Focal Node Centrality', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );

% axis tight;
apx2 = '';
if distance_mode
  apx2 = '-dist';
end

set(gca,'DataAspectRatioMode', 'manual', 'DataAspectRatio',[1 8000 50]);
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
end

%% Plot geographic map

% figure( 'position', [460 369 1110 564] );
% scatter( 1 ./ path_lengths(~in_graph), mean_inc(~in_graph), num_members(~in_graph) .* member_scaler, [.7 .7 .7], 'filled', 'markerfacealpha', marker_alpha );
% hold on;
% scatter( 1 ./ path_lengths(in_graph), mean_inc(in_graph), num_members(in_graph) .* member_scaler, 'b', 'filled', 'markerfacealpha', marker_alpha );
% 
% set( gca, 'fontsize', 12, 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
% % xlabel( 'Centrality', 'interpreter', 'latex' );
% xlabel( 'Focal Node Centrality', 'interpreter', 'latex' );
% % ylabel( 'Prevalence', 'interpreter', 'latex' );
% ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );

% Plot the chains going towards the `superbug'
% figure( 'position', [460 369 1110 564] );
% scatter( 1 ./ path_lengths, mean_inc, num_members .* member_scaler, [.7 .7 .7], 'filled', 'markerfacealpha', marker_alpha );
% hold on;

if plot_geographic
  nsw_lat = [-37.75 -28];
  nsw_long = [140.75 154];

  shape_coords = dlmread( [ database 'POA_2016_AUST_XY.csv' ], ',', 1, 0 );
  pc_to_latlongarea = containers.Map( num2cell( shape_coords( :, 1 ), 2 ),...
                                  num2cell( shape_coords( :, [5 4 3] ), 2 ) );

  % Ignore entries that are:
  % outside NSW; not in postcode map; weren't loaded correctly from the CSV (hacky - fix)
  valid_idxs = postcodes >= 2000 & postcodes <= 3000 ...
                  & isKey( pc_to_latlongarea, num2cell( postcodes ) );

  M_0 = sum( valid_idxs );

  % Scatter people around postcodes so we don't have overlap
  lat_long_areas = nan(M,3);
  lat_long_areas(valid_idxs,:) = cell2mat( values( pc_to_latlongarea, num2cell( postcodes(valid_idxs) ) ) );
  latlongs = lat_long_areas( :, 1:2 );

  figure( 'units', 'pixels', 'position', [732 1 1920 1090] );
  cmap = colormap('hot');
  cmap = cmap( 1:end-25, : );
  set( gca, 'color', [0.9 0.9 1.0] )
  postcode_shp_file = [ database '/shape_files/POA_2016_AUST.shp' ];
  ste_shp_file = [ database '/shape_files/STE06aAUST.shp' ];

  postcode_shps = shaperead( postcode_shp_file,...
                              'UseGeoCoords', true,...
                              'Selector', {@(v) (str2double(v) >= 2000) && (str2double(v) <= 3000), 'POA_CODE16' } );
  ext_shps = shaperead( postcode_shp_file,...
                              'UseGeoCoords', true,...
                              'Selector', {@(v) (str2double(v) < 2000) || (str2double(v) > 3000), 'POA_CODE16' } );

  pch = geoshow( ext_shps );
  pc_ph = get( pch, 'children' );
  set( pc_ph, 'FaceColor', [255 255 225] ./ 255, 'FaceAlpha', 1, 'EdgeColor', [.7 .7 .7] );

  hold on;

  pch = geoshow( postcode_shps );
  pc_ph = get( pch, 'children' );
  set( pc_ph, 'FaceColor', [255 255 225] ./ 255, 'FaceAlpha', 1 );
  xticklabels = get( gca, 'xticklabels' );
  for i = 1:length( xticklabels )
      xticklabels{i} = [ '$' xticklabels{i} '^{\circ}$' ];
  end
  yticklabels = get( gca, 'yticklabels' );
  for i = 1:length( yticklabels )
      yticklabels{i} = [ '$' yticklabels{i} '^{\circ}$' ];
  end
  set( gca, 'Layer', 'top', 'fontsize', 14, 'tickdir', 'out', 'ticklabelinterpreter', 'latex',...
              'xticklabels', xticklabels, 'yticklabels', yticklabels );
  box( gca, 'on' );

  % sh = scatter( latlongs( sb_chain, 2 ), latlongs( sb_chain, 1 ), 250, 'r', 'filled' );
  % alpha( sh, 0.7 );

  [~,max_fs] = sort(y_f(ell>10)-y_i(ell>10),'descend');
  max_m_fs = unique(m_f(max_fs));

  cols = [1 0 0;
          1 .5 .8;
          .2 .3 .8;
          .6 .3 1];
  for sb_id = 1:4
    focus_node = max_m_fs(sb_id);

    superbug_chains = find(m_f == focus_node);
    [~,max_sb_chain] = max(ell(superbug_chains));

%     sb_chain = unique_paths{superbug_chains(max_sb_chain)};
    sb_chain = unique_paths{6717};

    for m_i = 1:length(sb_chain)-1
      p1 = latlongs( sb_chain(m_i), [2 1] );
      p2 = latlongs( sb_chain(m_i+1), [2 1] );
      dp = p2-p1;
      quiver( p1(1), p1(2), dp(1), dp(2), 0, 'linewidth', 2, 'color', cols(sb_id,:), 'marker', '.', 'markersize', 20 );
    end
    break;
  end

  axis equal;
  set( gca, 'xlim', nsw_long, 'ylim', nsw_lat );
end

%% Plot correlations

figure('position',[257   660   570   357]);
plot( net.pdist( in_graph, mp_node ), clusters.mean_inc( in_graph ), 'k.' );

set( gca, 'fontsize', 12, 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
xlabel( 'Distance to Most Prevalent Node', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );

% curve = fit( D( in_graph, mp_node ), log10(mean_inc( in_graph )), 'poly1' );
% hold on;
% min_X = min(D( in_graph, mp_node ));
% max_X = max(D( in_graph, mp_node ));
% plot( [min_X, max_X], [curve(min_X), curve(round(max_X/2))], 'r-' );

[rho,pval] = corr( net.pdist( in_graph, mp_node ), log10(clusters.mean_inc(in_graph)) )


print( gcf, ['dist_prev' appendix apx2 '.pdf'], '-dpdf', '-bestfit' );

figure('position',[832   659   570   358]);
plot( size_comp(in_graph), clusters.mean_inc(in_graph), 'k.' );

set( gca, 'fontsize', 12, 'xscale', 'log', 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
xlabel( 'Size of Component', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );

[rho,pval] = corr( log10(size_comp(in_graph)), log10(clusters.mean_inc(in_graph)) )
