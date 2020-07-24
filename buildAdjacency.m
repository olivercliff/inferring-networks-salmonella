function A = buildAdjacency(net,event_dates,event_profile_ids,distances,cdfs,max_distance,alpha,bf_correct_per_distance)

if nargin < 8
  bf_correct_per_distance = false;
  if nargin < 7
    alpha = 0.05;
  end
end

Ds = 1:max_distance;

A = zeros(size(net.profile,1));

for m = 1:size(net.profile,1)
  
  first_date = event_dates(find(event_profile_ids == m,1,'first'))';
  
  if ~isempty(first_date)
    
    candidate_time_ids = find(distances(:,m) > 0 & distances(:,m) <= max_distance);

    if isempty(candidate_time_ids)
      continue;
    end
    
    candidate_distances = distances(candidate_time_ids,m);
    [candidate_profile_ids,uniq_ids] = unique(event_profile_ids(candidate_time_ids));
    candidate_distances = candidate_distances(uniq_ids);
    
    % Bonferroni correct for number of tests (per distance or total?)
    if bf_correct_per_distance
      n_tests = sum(Ds == candidate_distances,1);
      bf_alphas = alpha ./ n_tests;
    else
      n_tests = length(candidate_distances);
      bf_alpha = alpha ./ n_tests;
    end
    
    for c_id = 1:length(candidate_profile_ids)  
      % Obtain all dates of this candidate profile...
      all_dates = event_dates(event_profile_ids == candidate_profile_ids(c_id));
      
      % ...and get the difference to the first date of the focus profile
      day_diffs = days(all_dates - first_date);
      cdist = candidate_distances(c_id);
      
      % Bonferroni correct the 
      if bf_correct_per_distance
        bf_alpha = bf_alphas(cdist);
      end
      
      pre_window = find(cdfs.cdm_ptg(:,cdist)>=bf_alpha,1,'first');
      post_window = find(cdfs.cdm_ftg(:,cdist)>=bf_alpha,1,'first');

      % Which nodes occured immediately before (parent nodes)
      parent_node = any( day_diffs >= -pre_window & day_diffs <= 0 );

      % Which nodes co-evolved (sibling nodes)
      sibling_node = any( day_diffs <= post_window & day_diffs > 0 ) && day_diffs(1) >= -pre_window;
      if parent_node || sibling_node
        A(candidate_profile_ids(c_id),m) = 1;
      end
    end
  end

end