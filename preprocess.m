function preprocess(csvfile)

approaches = {'unadjusted','Larsson_plus','Larsson_joint'};

addpath ./utils

database = fileparts(csvfile);

dat = readtable(csvfile);
N = height(dat);

% Read in data, throwing out the incorrectly coded values
invalid = false(N,1);
mlva_profile = zeros(N,5);
name = cell(N,1);
for i = 1:N
  name{i} = dat.MLVAType{i};
  out = sscanf(name{i},'%d-%d-%d-%d-%d',[1,5]);
  try
    mlva_profile(i,:) = out;
  catch
    warning('Row %d incorrectly input.\n',i);
    invalid(i) = true;
  end
end

mlva_profile(invalid,:) = [];
name(invalid) = [];

% Initialise the events table with the genotype and date collected
mlva_type = dat.MLVAType(~invalid);
date = dat.DateCollected(~invalid);

events.dates = date;
events.profile_name = mlva_type;

% From Larsson's paper
fragment_size = [337 370 436 451 463 469 490 496 517 523 544 550 572 616];
repeats_27bp = [0 0 0 3 1 0 2 1 3 2 4 3 5 3];
repeats_33bp = [8 9 11 9 11 12 11 12 11 12 11 12 11 14];
unknown_sizes = setdiff( mlva_profile( :, 5 ), fragment_size )';

for i = 1:length(approaches)

  fprintf('Preprocessing the %s approach...\n', approaches{i});
  
  % There are a few different ways to define the MLVA profiles (using
  % Larsson's approach)
  if i == 1
    % Don't adjust, i.e., use the last loci
    fragment_size = [ fragment_size, unknown_sizes ];
    sz_cell = num2cell( fragment_size );
    sz_to_repeats = containers.Map( sz_cell, sz_cell );
  elseif i == 2
    % Adjust the last loci so that it's a summation of the 27&33 BPs
    repeats_summed = repeats_27bp + repeats_33bp;
    repeats_summed = [ repeats_summed, round( unknown_sizes ./ 37 ) ];
    sz_to_repeats = containers.Map( sz_cell, repeats_summed );
  elseif i == 3
    % Adjust include 2 extra locis so that we have both the 27&33 BPs
    repeats_joint = [repeats_27bp; repeats_33bp];
    repeats_joint = [ repeats_joint, [ zeros(1,length(unknown_sizes)); round(unknown_sizes./37)] ];
    joint_cell = num2cell(repeats_joint',2);
    sz_to_repeats = containers.Map( sz_cell, joint_cell );
  end

  vals = cell2mat( values( sz_to_repeats, num2cell( mlva_profile( :, 5 ) ) ) );

  capproach_profile = [ mlva_profile(:,1:4), vals ];
  [unique_profiles,~,id] = unique(capproach_profile,'rows');
  
  % Map the event...
  events.(['id_' approaches{i}]) = id;
  
  % ...to the unique profile
  network.(approaches{i}).profile = unique_profiles;
  
  network.(approaches{i}).names = cell(length(unique_profiles),1);
  for p = 1:size(unique_profiles,1)
    cname = [];
    for l = 1:size(unique_profiles,2)-1
      cname = [cname sprintf('%d-',network.(approaches{i}).profile(p,l))];
    end
    network.(approaches{i}).names{p} = [cname sprintf('%d',network.(approaches{i}).profile(p,end))];
  end
  
  % Number of times each profile was recorded
  network.(approaches{i}).n_incidences = histcounts(id,(1:size(unique_profiles,1)+1)-0.5);
  
  % Pairwise distances between all profiles using the L1-norm
  pd = squareform( pdist( unique_profiles, 'cityblock' ) );
%   pd(pd < 1 & pd > 0) = 1;
  network.(approaches{i}).pdist = pd;

  % Get edge weights
  ew = 1 ./ network.(approaches{i}).pdist;
  ew( isinf( ew ) ) = 1;
  ew( ew == 1 ) = 1 - 1e-16;

  % Make (undirected) graph to use matlab centrality functions
  G = graph( ew, 'OmitSelfLoops' );
  G_0 = graph( network.(approaches{i}).pdist, 'OmitSelfLoops' );

  % Compute a bunch of network stats
  network.(approaches{i}).deg_ranks = centrality( G, 'degree', 'Importance', G.Edges.Weight );
  network.(approaches{i}).closeness = centrality( G, 'closeness', 'Cost', 1 ./ G.Edges.Weight );
  network.(approaches{i}).betweenness = centrality( G, 'betweenness', 'Cost', 1 ./ G.Edges.Weight );
  network.(approaches{i}).pagerank = centrality( G, 'pagerank', 'Importance', G.Edges.Weight );
  [~, network.(approaches{i}).clustering_coefficients] = weightedClusteringCoefficient( network.(approaches{i}).pdist );
  network.(approaches{i}).path_lengths = mean( distances( G_0 ) )';
  network.(approaches{i}).small_world_propensity = network.(approaches{i}).clustering_coefficients ./ network.(approaches{i}).path_lengths;
  
  fprintf('Done.\n');
end

save( [database '/network.mat'], 'network', 'events' ); 