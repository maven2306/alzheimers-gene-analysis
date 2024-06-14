function indx_all = make_fair(list_patients, list_set)
% function indx_all = make_fair(list_patients, list_set)
% Make control/patient ratio equal.

indx_pat = find(list_set & list_patients);
indx_con = find(list_set & ~list_patients);


% Countnumber of subjects in largest group.
npat = length(indx_pat);
ncon = length(indx_con);
nmax = max(npat, ncon);

% split controls and patients into two variables.

% repeat the smallest group x times, such that has more datapoints than
% the other (previously larger) group. And cut, it such that both have
% equal length.
indx_con = repmat(indx_con, ...
    [ceil(npat ./ ncon) 1]);
indx_con = indx_con(1:nmax);

indx_pat = repmat(indx_pat, ...
    [ceil(ncon ./ npat) 1]);
indx_pat = indx_pat(1:nmax);

% Put controls and patients indices back together.
indx_all = [indx_con; ...
    indx_pat];

end