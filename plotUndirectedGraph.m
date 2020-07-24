function plotUndirectedGraph(matfile,approach)

load( matfile, 'network' );

mynet = network.(approach);

marker_sz = 30;
marker_col = [.4 .4 .4];
marker_alpha = .7;

%% Overlapping Approach

[ ~, focus_node ] = max( mynet.mean_inc );

member_scaler = 0.8;

figure( 'position', [460 369 1110 564] );
scatter( 1 ./ mynet.path_lengths, mynet.mean_inc,...
          mynet.num_members .* member_scaler, log10( mynet.pdist( :, focus_node ) ),...
          'filled', 'markerfacealpha', marker_alpha );
        
set( gca, 'fontsize', 12, 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
xlabel( 'Focal Node Centrality', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );

edges = linspace( min( 1 ./ mynet.path_lengths ), max( 1 ./ mynet.path_lengths ), 50 );
bins = discretize( 1 ./ mynet.path_lengths, edges );
avgs = zeros( length( edges ), 1 );
stds = zeros( length( edges ), 1 );
for i = 1:length( edges )
    avgs(i) = mean( mynet.mean_inc( bins == i ) );
    stds(i) = std( mynet.mean_inc( bins == i ) );
end

hold on;
plot( edges, movmean( avgs, 5, 'omitnan' ), 'k.-' );
plot( edges, movmean( avgs+stds, 5, 'omitnan' ), 'k--' );
% plot( edges, movmean( avgs-0.5.*stds, 5, 'omitnan' ), 'k--' );

cbh = colorbar;
set( cbh, 'fontsize', 12, 'limits', [0 max(get(cbh, 'limits'))], 'ticklabelinterpreter', 'latex', 'ticks', [0 1], 'ticklabels', {'$10^0$', '$10^1$'} );
axis tight;
print( gcf, sprintf( './plots/trace-evolution%s-D%d.pdf', appendix, cluster_D ), '-dpdf', '-bestfit' );
print( gcf, sprintf( './plots/trace-evolution%s-D%d.eps', appendix, cluster_D ), '-depsc2', '-r0', '-painters' );

if ~exist( 'mlvas_adjusted', 'var' )
    mlvas_adjusted = mlvas( :, 1:4 );
end

figure( 'position', [460 369 1110 564] );
for i = 1:size( mlvas_adjusted, 2 )
    cD = squareform( pdist( mlvas_adjusted( :, i ), 'cityblock' ) );
    
    scatter( 1 ./ path_lengths, mean_inc, num_members .* member_scaler, log10( cD( :, focus_node ) ), 'filled', 'markerfacealpha', marker_alpha );
    set( gca, 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
    xlabel( 'Focal Node Centrality', 'interpreter', 'latex' );
    ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );
    
    cbh = colorbar;
    set( cbh, 'ticklabelinterpreter', 'latex' );
    axis tight;
    print( gcf, sprintf( './plots/trace-evolution%s-L%d-D%d.pdf', appendix, i, cluster_D ), '-dpdf', '-bestfit' );
    print( gcf, sprintf( './plots/trace-evolution%s-L%d-D%d.eps', appendix, i, cluster_D ), '-depsc2', '-r0', '-painters' );
end

save( [ 'overlapping-clusters' appendix '.mat'], 'cluster_ids', 'num_members', 'mean_cen', 'mean_inc' ); 

%% Dendogram Approach

max_cluster = 21;
marker_sz = 150;

cluster_idxs = clusterdata( mlvas_adjusted, 'distance', 'cityblock', 'linkage', 'weighted', 'maxclust', max_cluster );
uniq_clusters = unique( cluster_idxs );

cluster_ids = cell( length( uniq_clusters ), 1 );
sum_centrality = zeros( length( uniq_clusters ), 1 );
sum_incidence = zeros( length( uniq_clusters ), 1 );

for i = 1:length( uniq_clusters )
    cluster_ids{i} = find( cluster_idxs == uniq_clusters(i) );
%     sum_centrality(i) = sum( pagerank( idxs{i} ) );
    sum_centrality(i) = sum( 1 ./ path_lengths( cluster_ids{i} ) );
    sum_incidence(i) = sum( n_incidences( cluster_ids{i} ) );
end

num_members = cellfun( @numel, cluster_ids );
mean_cen = sum_centrality ./ num_members;
mean_inc = sum_incidence ./ num_members;

% figure; hist( num_members );

[ ~, focus_cluster ] = max( mean_inc );

cdist = zeros( length( cluster_ids ), 1 );
for i = 1:length( cluster_ids )
    dist_mat = pdist2( mlvas_adjusted( cluster_ids{focus_cluster}, : ), mlvas_adjusted( cluster_ids{i}, : ), 'cityblock' );
    cdist(i) = mean( dist_mat(:) );
end

% member_scaler2 = 20;

figure( 'position', [460 369 1110 564] );
scatter( mean_cen, mean_inc, num_members .* member_scaler, log10( cdist ), 'filled' );
set( gca, 'fontsize', 12, 'yscale', 'log', 'box', 'on', 'ticklabelinterpreter', 'latex' );
xlabel( 'Average Cluster Centrality', 'interpreter', 'latex' );
ylabel( 'Average Cluster Prevalence', 'interpreter', 'latex' );

edges = linspace( min( mean_cen ), max( mean_cen ), 10 );
bins = discretize( mean_cen, edges );
avgs = zeros( length( edges ), 1 );
stds = zeros( length( edges ), 1 );
for i = 1:length( edges )
    avgs(i) = mean( mean_inc( bins == i ) );
    stds(i) = std( mean_inc( bins == i ) );
end

hold on;
plot( edges, movmean( avgs, 5, 'omitnan' ), 'k.-' );
plot( edges, movmean( avgs+stds, 5, 'omitnan' ), 'k--' );

cbh = colorbar;
set( cbh, 'fontsize', 12, 'limits', [0 max(get(cbh, 'limits'))], 'ticklabelsmode', 'auto', 'ticklabelinterpreter', 'latex', 'ticks', [0 1], 'ticklabels', {'$10^0$', '$10^1$'} );
axis tight;

print( gcf, sprintf( './plots/dend-trace-evolution%s-M%d.pdf', appendix, max_cluster ), '-dpdf', '-bestfit' );
print( gcf, sprintf( './plots/dend-trace-evolution%s-M%d.eps', appendix, max_cluster ), '-depsc2', '-r0', '-painters' );

save( [ 'hierarchical-clusters' appendix '.mat'], 'cluster_ids', 'num_members', 'mean_cen', 'mean_inc' ); 