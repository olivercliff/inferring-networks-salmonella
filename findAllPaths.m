function [paths,visited] = findAllPaths(G,s)
% Find all paths in digraph G. Known NP-hard problem so it might not finish
% Inputs
%   - matlab G digraph with N nodes
%   - s source node from which to search (1 <= s <= N)
% Outputs
%   - paths is cell array of all paths found outgoing from node s
%   - visited is an array of size Nx1 that gives which nodes exist in any
%   array in paths

N = G.numnodes;
visited = false(N,1);

if visited(s)
    error('Already visited node %d\n', s);
end

c_i = successors(G,s);

% Unfortunately, no neat way to pre-allocate memory for `paths'
% (maybe change to allocating blocks at a time)
paths = cell(length(c_i),1);
visited(s) = true;

if isempty(paths)
    return;
end

for i = 1:length(c_i)
    paths{i} = [s c_i(i)];
end

% Same size as `paths', just determines whether we can explore more
at_root = false(length(c_i),1);

while ~all(at_root)
    for i = 1:length(paths)
        if ~at_root(i)
            cpath = paths{i};
            p_i = cpath(end);
            tts = successors(G,p_i);
            
            new_node_ids = ~ismember(tts,cpath);
            if isempty(tts) || ~any(new_node_ids)
                at_root(i) = true;
            else
                valid_nodes = tts(new_node_ids);
                for j = 1:length(valid_nodes)
                    % If there's only one, just add it, otherwise make new
                    % branch
                    if j == 1
                        paths{i} = [cpath valid_nodes(j)];
                    else
                        paths{end+1,1} = [cpath valid_nodes(j)];
                        at_root(end+1,1) = false;
                    end
                end
            end
            visited(p_i) = true;
        end
    end
end