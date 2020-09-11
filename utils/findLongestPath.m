function [path_ids,longest_paths] = findLongestPath(paths,ids)

% Find the longest paths that starts with these MLVAs
path_ids = nan(length(ids),1);
for i = 1:length(paths)
  [~,~,match] = intersect(paths{i}, ids);
  if ~isempty(match)
    if isnan(path_ids(match)) || length(paths{i}) > length(paths{path_ids(match)})
      path_ids(match) = i;
    end
  end
end

longest_paths = paths(path_ids);