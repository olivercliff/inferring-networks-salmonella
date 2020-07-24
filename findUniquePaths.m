function unique_paths = findUniquePaths(all_paths)

fprintf('Converting paths to strings for matching...\n');

% Easier to match substrings (no inbuilt command for checking ordered
% subsets? I'm also just being lazy by doing this conversion)
P = length(all_paths);
str_paths = cell(P,1);
parfor i = 1:P
    str_paths{i} = mat2str(all_paths{i});
    
    if mod(i,1000) == 0
        fprintf('%d/%d\n', i, P);
    end
end
fprintf('Done.\n');

fprintf('Pruning non-overlapping paths...\n')

uniq_ids = true(P,1);
parfor i = 1:P
    for j = 1:P
        if i == j
            continue;
        end
        
        if isequal(all_paths{i},all_paths{j})
            if i > j
                uniq_ids(i) = false;
            end
            continue;
        end
        
        if contains(str_paths{j}(2:end-1),str_paths{i}(2:end-1))
            uniq_ids(i) = false;
            break;
        end
    end
    if mod(i,1000) == 0
        fprintf('%d/%d.\n', i, P);
    end
end

unique_paths = all_paths(uniq_ids);
fprintf('Done.\n');