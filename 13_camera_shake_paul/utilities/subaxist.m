function h=subaxist(rows,cols,celly,cellx,settings)
    % wrapper for subaxis function. transposes the indexing
    % and allows you to pass the cell settings array in w/o doing {:}.
    if nargin < 5, settings = {}; end
    h = subaxis(rows,cols,cellx,celly,settings{:});
end