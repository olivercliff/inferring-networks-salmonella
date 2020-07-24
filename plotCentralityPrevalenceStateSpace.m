function plotCentralityPrevalenceStateSpace(matfile,approach,member_scaler)

load( matfile, 'network' );

cluster_D = network.cluster_D;
max_cluster = network.max_cluster;
mynet = network.(approach);

marker_col = [.4 .4 .4];
marker_alpha = .7;
if nargin < 3
  member_scaler = 0.8;
end

%% Overlapping Approach

marker_sz = 30;

clusters = mynet.clusters.overlapping;

[ ~, focus_node ] = max( clusters.mean_inc );

figure( 'position', [460 369 1110 564] );
scatter( 1 ./ mynet.path_lengths, clusters.mean_inc,...
          clusters.num_members .* member_scaler, log10( mynet.pdist( :, focus_node ) ),...
          'filled', 'markerfacealpha', marker_alpha );
        
set( gca, 'fontsize', 12, 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
xlabel( 'Focal Node Centrality', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );

edges = linspace( min( 1 ./ mynet.path_lengths ), max( 1 ./ mynet.path_lengths ), 50 );
bins = discretize( 1 ./ mynet.path_lengths, edges );
avgs = zeros( length( edges ), 1 );
stds = zeros( length( edges ), 1 );
for i = 1:length( edges )
    avgs(i) = mean( clusters.mean_inc( bins == i ) );
    stds(i) = std( clusters.mean_inc( bins == i ) );
end

hold on;
plot( edges, movmean( avgs, 5, 'omitnan' ), 'k.-' );
plot( edges, movmean( avgs+stds, 5, 'omitnan' ), 'k--' );
% plot( edges, movmean( avgs-0.5.*stds, 5, 'omitnan' ), 'k--' );

cbh = colorbar;
set( cbh, 'fontsize', 12, 'limits', [0 max(get(cbh, 'limits'))], 'ticklabelinterpreter', 'latex', 'ticks', [0 1], 'ticklabels', {'$10^0$', '$10^1$'} );
axis tight;
print( gcf, sprintf( './plots/trace-evolution%s-D%d.pdf', approach, cluster_D ), '-dpdf', '-bestfit' );
print( gcf, sprintf( './plots/trace-evolution%s-D%d.eps', approach, cluster_D ), '-depsc2', '-r0', '-painters' );

%% Dendogram Approach

clusters = mynet.clusters.partition;

[ ~, focus_cluster ] = max( clusters.mean_inc );

cdist = zeros( length( clusters.ids ), 1 );
for i = 1:length( clusters.ids )
    dist_mat = mynet.pdist( clusters.ids{focus_cluster}, clusters.ids{i} );
    cdist(i) = mean( dist_mat(:) );
end

% member_scaler2 = 20;

figure( 'position', [460 369 1110 564] );
scatter( clusters.mean_cen, clusters.mean_inc,...
          clusters.num_members .* member_scaler, log10( cdist ), 'filled' );
set( gca, 'fontsize', 12, 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
xlabel( 'Average Cluster Centrality', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );

edges = linspace( min( clusters.mean_cen ), max( clusters.mean_cen ), 10 );
bins = discretize( clusters.mean_cen, edges );
avgs = zeros( length( edges ), 1 );
stds = zeros( length( edges ), 1 );
for i = 1:length( edges )
    avgs(i) = mean( clusters.mean_inc( bins == i ) );
    stds(i) = std( clusters.mean_inc( bins == i ) );
end

hold on;
plot( edges, movmean( avgs, 5, 'omitnan' ), 'k.-' );
plot( edges, movmean( avgs+stds, 5, 'omitnan' ), 'k--' );

cbh = colorbar;
set( cbh, 'fontsize', 12, 'limits', [0 max(get(cbh, 'limits'))], 'ticklabelsmode', 'auto', 'ticklabelinterpreter', 'latex', 'ticks', [0 1], 'ticklabels', {'$10^0$', '$10^1$'} );
axis tight;

print( gcf, sprintf( './plots/dend-trace-evolution%s-M%d.pdf', approach, max_cluster ), '-dpdf', '-bestfit' );
print( gcf, sprintf( './plots/dend-trace-evolution%s-M%d.eps', approach, max_cluster ), '-depsc2', '-r0', '-painters' );