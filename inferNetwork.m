function inferNetwork(matfile,approach,multiple_loci,alpha)

if nargin < 4
  alpha = 0.05;
end

switch multiple_loci
  case 0
    fld = 'singleloci';
  case 1
    fld = 'multiloci';
  case 2
    fld = 'ubloci';
end

% Load and preprocess data
load( matfile, 'network', 'events' );

mynet = network.(approach);

cdfs = mynet.(fld).cdfs;

%% Build or load adjacency matrix

% Build adjacency matrix
fprintf('Building adjacency matrix...\n');

event_profile_ids = events.(['id_' approach]);

adj = buildAdjacency(mynet,events.dates,event_profile_ids,events.(approach).(fld).genetic_distances,cdfs,4,alpha);
fprintf('Done.\n');

fprintf('Nodes in the directed network: %d\n', sum(sum(adj,2)>0 | sum(adj,1)'>0));

%% Find all paths

% Create node names
M = size(mynet.profile,1);
node_name = cell(M,1);

for i = 1:M
  node_name{i} = sprintf('%d-%d-%d-%d-%d',mynet.profile(i,:));
end

G = digraph(adj,node_name,'OmitSelfLoops'); % Create graph object

% Having two for loops like this should speed up the parallelisation a bit
tic

% Find all potential (evolutionary) paths through the adjacency
paths = cell(M,1);
parfor m_i = 1:M
  paths{m_i} = findAllPaths(G, m_i);
  fprintf('Computed all paths for source node %d/%d\n', m_i, M);
end

all_paths = {};
for m_i = 1:M
  all_paths = [all_paths; paths{m_i}];
end
unique_paths = findUniquePaths(all_paths);

toc

mynet.(fld).adj_mat = adj;
mynet.(fld).unique_paths = unique_paths;

% Save adjacency and paths back to network file
network.(approach) = mynet;
save( matfile, 'network', '-append' );