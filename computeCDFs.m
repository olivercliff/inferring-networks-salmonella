function computeCDFs(matfile,approach,multiple_loci,up_to_ml)
% Build the cumua

if nargin < 5
  up_to_ml = false;
end

if multiple_loci
  fld = 'multiloci';
else
  fld = 'singleloci';
end

load( matfile, 'network', 'events' );
mynet = network.(approach);

mynet.(fld).cdfs.up_to_ml = up_to_ml;

%% Build frequency distributions

M = size(mynet.profile,1); % Number of unique profiles
L = size(mynet.profile,2); % Number of locis
T = length(events.dates); % Length of dataset

ds = [0.5 1:L inf];

mynet.(fld).cdfs.distance_edges = ds;

edges = 0:1:T;
mynet.(fld).cdfs.temporal_edges = edges;

dp_hist = zeros(length(edges)-1,L,M);
df_hist = zeros(length(edges)-1,L,M);

indp_hist = zeros(length(edges)-1,L,M);
indf_hist = zeros(length(edges)-1,L,M);

fprintf('Building frequency distributions...\n');

events.(approach).occurances = false(T,M);
events.(approach).genetic_distances = inf(T,M);

event_profile_ids = events.(['id_' approach]);
event_profiles = mynet.profile(event_profile_ids,:);

for m = 1:M
  
  cmlva_adjusted = mynet.profile(m,:);
  events.(approach).occurances(:,m) = event_profile_ids == m;
  cdates = events.dates(events.(approach).occurances(:,m))';
  
  if length(cdates)>1
    % Take the difference between MLVA profiles
    loci_changes = abs( bsxfun(@minus, event_profiles, cmlva_adjusted) );

    n_loci_changed = L - sum(loci_changes == 0, 2 );
    l1_norm = sum(loci_changes, 2);

    if multiple_loci
      if up_to_ml
        l1_norm(n_loci_changed > l1_norm) = inf;
      else
        l1_norm(n_loci_changed ~= l1_norm) = inf;
      end
    else
      l1_norm(n_loci_changed > 1) = inf;
    end
    events.(approach).genetic_distances(:,m) = l1_norm;

    for l = 1:L

      % Find MLVAs that are within genetic distance L
      valid_dates = events.dates(l1_norm>ds(l) & l1_norm<=ds(l+1));

      if isempty(valid_dates)
        continue;
      end

      diffs_p = nan(length(cdates),1);
      diffs_f = nan(length(cdates),1);
      for i = 1:length(cdates)
        pi = find(valid_dates<cdates(i),1,'first');
        fi = find(valid_dates>cdates(i),1,'first');

        if ~isempty(pi)
          % Occurred before this instance
          diffs_p(i) = days(cdates(i) - valid_dates(pi));
        end
        if ~isempty(fi)
          % Occurred after this instance
          diffs_f(i) = days(valid_dates(fi) - cdates(i));
        end
      end

      indp_hist(:,l,m) = histcounts(diffs_p(2:end),edges);
      indf_hist(:,l,m) = histcounts(diffs_f(2:end),edges);

      % Include after FIRST instance (potential mutation) and...
      dp_hist(:,l,m) = histcounts(diffs_p,edges);

      % ...include before FIRST instance (potential mutation).
      df_hist(:,l,m) = histcounts(diffs_f,edges);
    end
  end

end
fprintf('Done.\n');

sindp_hist = nansum(indp_hist,3);
sindf_hist = nansum(indf_hist,3);

sdp_hist = nansum(dp_hist,3);
sdf_hist = nansum(df_hist,3);

% p(past T|G,~M), probability of temporal differences given NOT a mutation
% and genetic distance G
ppt_gm = sindp_hist ./ sum(sindp_hist);
ppt_gm = ppt_gm ./ sum(ppt_gm);

% p(future T|G,~M), probability of temporal differences given NOT a
% mutation and genetic distance G
pft_gm = sindf_hist ./ sum(sindf_hist);
pft_gm = pft_gm ./ sum(pft_gm);

% p(past T|G), probability of past temporal differences given genetic
% distance G
ppt_g = sdp_hist ./ sum(sdp_hist);
ppt_g = ppt_g ./ sum(ppt_g);

% p(future T|G), probability of future temporal differences given genetic
% distance (i.e., G>0)
pft_g = sdf_hist ./ sum(sdf_hist);
pft_g = pft_g ./ sum(pft_g);

% f(~M | past T, G), frequency distribution of NOT a mutation, given
% past temporal differences and genetic distance
fm_ptg = ppt_gm./ppt_g;

% f(~M | future T, G), frequency distribution of NOT a mutation, given
% future temporal differences and genetic distance
fm_ftg = pft_gm./pft_g;

% p(~M | G ), probability of NOT a mutation given genetic distance)
pm_g_pt = zeros(1,L);
pm_g_ft = zeros(1,L);
for l = 1:L
  pm_g_pt(l) = 1./nansum( fm_ptg(~isinf(fm_ptg(:,l)),l) );
  pm_g_ft(l) = 1./nansum( fm_ftg(~isinf(fm_ftg(:,l)),l) );
end

pm_ptg = pm_g_pt.*fm_ptg;
pm_ptg(isinf(pm_ptg)) = nan;

pm_ftg = pm_g_ft.*fm_ftg;
pm_ftg(isinf(pm_ftg)) = nan;

% Cumulate along the genetic-distance dimensions
cdm_ptg = cumsum(pm_ptg,'omitnan');
cdm_ftg = cumsum(pm_ftg,'omitnan');

% Cumulate along the time-delay dimension
% (Past)
cdm_ptg = cumsum(cdm_ptg,2);
mynet.(fld).cdfs.cdm_ptg = cdm_ptg./cdm_ptg(end);

% (Future)
cdm_ftg = cumsum(cdm_ftg,2);
mynet.(fld).cdfs.cdm_ftg = cdm_ftg./cdm_ftg(end);

network.(approach) = mynet;

save(matfile,'network','events','-append');