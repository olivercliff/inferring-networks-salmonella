fn1 = nan;
fn2 = nan;
fn3 = 346;

in_paths = false(U,1);
for i = 1:U
  if sum(ismember(unique_paths{i},[fn1 fn2 fn3])) == 1
    in_paths(i) = true;
  end
end

% in_paths = false(U,1);
% for i = 1:U
%   if any(ismember(unique_paths{i},fn3))
%     in_paths(i) = true;
%   end
% end