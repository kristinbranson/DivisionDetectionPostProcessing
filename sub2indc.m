function i = sub2indc(sz,c)

nd = numel(c);
ns = cellfun(@numel,c);
if any(ns == 0),
  i = [];
  return;
elseif nnz(ns > 1) == 1,
  n = ns(ns>1);
  for j = find(ns(:)'<=1),
    c{j} = repmat(c{j},[n,1]);
  end
else
  [c{ns>1}] = ndgrid(c{ns>1});
  n = numel(c{find(ns>1,1)});
  for j = 1:nd,
    if ns(j) == 1,
      c{j} = repmat(c{j},[n,1]);
    else
      c{j} = c{j}(:);
    end
  end  
end

i = sub2indv(sz,cat(2,c{:}));
