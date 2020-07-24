function getClusters(matfile,cluster_D,max_cluster)

% Configure the cluster sizes
if nargin < 3
  max_cluster = 21; % Max number of non-overlapping (agglomerative) clusters
  if nargin < 2
    cluster_D = 5; % Max distance for overlapping clusters
  end
end

%% Load data
try
  load( matfile, 'network' );
catch
  warning('MAT-file %s does not exist yet. Have you run ''preprocess.m''?', matfile);
  return;
end

network.max_cluster = max_cluster;
network.cluster_D = cluster_D;

approaches = {'unadjusted','Larsson_plus','Larsson_joint'};

for i = 1:length(approaches)
  approach = approaches{i};

  mynet = network.(approach);
  N = length(mynet.profile);

  %% Overlapping Approach

  cluster_ids = cell( N, 1 );
  sum_centrality = zeros( N, 1 );
  sum_incidence = zeros( N, 1 );
  for j = 1:N
      cluster_ids{j} = find( mynet.pdist( j, : ) <= cluster_D );
      sum_centrality(j) = nansum( 1 ./ mynet.path_lengths( cluster_ids{j} ) );
      sum_incidence(j) = nansum( mynet.n_incidences( cluster_ids{j} ) );
  end

  mynet.clusters.overlapping.num_members = cellfun( @numel, cluster_ids );
  mynet.clusters.overlapping.ids = cluster_ids;

  mynet.clusters.overlapping.num_members = cellfun( @numel, cluster_ids );
  mynet.clusters.overlapping.mean_cen = sum_centrality ./ mynet.clusters.overlapping.num_members;
  mynet.clusters.overlapping.mean_inc = sum_incidence ./ mynet.clusters.overlapping.num_members;

  %% Agglomerative Approach

  cluster_idxs = clusterdata( mynet.profile,...
                              'distance','cityblock',...
                              'linkage','weighted',...
                              'maxclust',max_cluster );
  uniq_clusters = unique( cluster_idxs );

  cluster_ids = cell( length( uniq_clusters ), 1 );
  sum_centrality = zeros( length( uniq_clusters ), 1 );
  sum_incidence = zeros( length( uniq_clusters ), 1 );

  for j = 1:length( uniq_clusters )
      cluster_ids{j} = find( cluster_idxs == uniq_clusters(j) );
      sum_centrality(j) = sum( 1 ./ mynet.path_lengths( cluster_ids{j} ) );
      sum_incidence(j) = sum( 1 ./ mynet.n_incidences( cluster_ids{j} ) );
  end

  mynet.clusters.partition.num_members = cellfun( @numel, cluster_ids );
  mynet.clusters.partition.ids = cluster_ids;
  mynet.clusters.partition.mean_cen = sum_centrality ./ mynet.clusters.partition.num_members;
  mynet.clusters.partition.mean_inc = sum_incidence ./ mynet.clusters.partition.num_members;

  %% Save back to network file

  network.(approach) = mynet;
end

save( matfile, 'network', '-append' ); 