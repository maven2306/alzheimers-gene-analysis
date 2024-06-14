function plotSurfaceBoth(regions, values, cm, varargin)
% PLOTSURFACEBOTH Create SVG of human hemisphere, with colored regions.
%
% PLOTSURFACEBOTH(regions, values, cm) - regions describes a list (Nx1) of
% regions to color. Values is an vector with the values associated with
% each region (Nx1). Variable cm is the colormap (Mx3) to which values are
% mapped
%
% PLOTSURFACEBOTH(regions, values, cm, ...) plots surfaces with
% optional arguments:
% 'limits',    array with [xmin xmax]
%              - default: [min(values) max(values)]
% 'viewer',    scalar with 1/true (show) or 0/false (no show)
%              - default: 1 (show)
% 'savePath'   path to save file 'savePath'_DK114_combined.svg
%              - default: temporary dir, deleted afterwards.
% 'scaling'    set scaling of image. Original scaling is very large, but
%              smaller scalings show small white lines
%              - default: '0.1' (10%)
% 'atlas'      chose which atlas to use ('DK114', 'MR250', 'BB38human', 'BB38chimp' or 'WBB47')
%              - default: 'DK114'

% PARSE INPUT
% parse obligatory input regions, values, cm
assert(iscellstr(regions), ...
    'regions must be array of strings with region names');
assert(isnumeric(values), 'values must be numeric');
regions = regions(:); % make sure it is a row-vector
assert(isvector(regions) && isvector(values) && ...
    numel(regions) == numel(values), ...
    'regions and values must be same size');
assert(isnumeric(cm) && (size(cm,2) == 3), ...
    'colormap must be numeric Nx3 matrix');

%-- HACK FOR BIGDATA COURSE: map DK68 regionvals to DK114
RegionLabels114=load('RegionNamesDK114.mat');
values2=zeros(1,114);
for ii=1:size(regions,1)
   tmp68=regions{ii};
   tmpDK114_indx=find(cellfun(@isempty, strfind(RegionLabels114.regionDescriptions,tmp68))==0); 
   values2(tmpDK114_indx) = values(ii);
end
regions=RegionLabels114.regionDescriptions;
values=values2;
%-- 

% parse optional arguments (limits and viewer, save path)
while ~isempty(varargin)
    if numel(varargin) == 1
        error('lscatter:missing_option', ...
            'Optional arguments must come in pairs.');
    end
    
    switch lower(varargin{1})
        case 'limits'
            assert(isnumeric(varargin{2}), (numel(varargin{2}) == 2))
            values_min = varargin{2}(1);
            values_max = varargin{2}(2);
        case 'viewer'
            assert(numel(varargin{2}) == 1, ...
                'viewer must be 0 (no show) or 1 (show)');
            viewer = varargin{2}>0;
        case 'savepath'
            assert(ischar(varargin{2}), ...
                'save_file must be string with file path');
            save_file = varargin{2};
        case 'scaling'
            assert(isnumeric(varargin{2}), ...
                'scaling must be a scalar > 0 and <= 1');
            assert((varargin{2} > 0) & (varargin{2} <= 1), ...
                'scaling must be a scalar > 0 and <= 1');
            scaling = varargin{2};     
        case 'atlas'
            assert(ischar(varargin{2}), ...
                'atlas must be string with file path');
            atlas = varargin{2};            
    end
    
    % remove the two entries we have just dealt with
    varargin(1:2) = [];
    
end

% set defaults
if ~exist('save_file', 'var')
    tmpDir = tempdir; % path to Documents
    tmpDir = [tmpDir 'matlab_plot_surface_' num2str(randi(10000))];
    mkdir(tmpDir);
    save_file = tmpDir;
end

if ~exist('values_min', 'var')
    values_min = nanmin(values);
    values_max = nanmax(values);
end

if ~exist('viewer', 'var')
    viewer = true;
end

if ~exist('scaling', 'var')
    scaling = 0.1;
end

if ~exist('atlas', 'var')
    atlas = 'DK114';
end



% INITIALIZE
cb_path = [save_file '_cb.png'];
combined_svg_path = fullfile(save_file,['plot_',atlas '_combined.svg']);
original_svg_path = [strrep(which(mfilename),[mfilename '.m'],''), ...
    atlas '_combined.svg'];

assert(exist(original_svg_path, 'file') > 0, ...
    'Atlas not found (%s)', original_svg_path);

switch atlas
    case 'WBB47'
        regions = strcat(regions, '_');
    case 'BB38human'
        regions = strcat(regions, '-');        
    case 'BB38chimp'
        regions = strcat(regions, '_');        
end

% COLORING
% from value2Color
% cut-off values at values_min and values_max
values(values<values_min) = values_min;
values(values>values_max) = values_max;

n_color = size(cm,1);
if length(unique(values)) > 1
    %normalize values to between 0 and 1
    values_norm = (values-values_min)./(values_max - values_min);
    
    % project on coloring
    values_color = round(values_norm*(n_color-1))+1;
    coloring = cm(values_color,:);
    
else
    % if only one value is specified, pick one color
    coloring = cm(repmat(n_color, [length(values) 1]), :);
end

% from cm2rgb
% write cm as hex rgb code
coloring_rgb = cell(size(coloring,1),1);
for i = 1:size(coloring,1)
    red = coloring(i,1);
    green = coloring(i,2);
    blue = coloring(i,3);
    htmlcolor = strcat('#', ...
        dec2hex(round(255*red), 2), ...
        dec2hex(round(255*green), 2), ...
        dec2hex(round(255*blue), 2));
    coloring_rgb(i,1) = cellstr(htmlcolor);
end

% SURFACES
makeSurface(original_svg_path, combined_svg_path, regions, coloring_rgb);

% COLORBAR
cm2(:,1,:) = cm;
cm2 = repmat(cm2, [ 1 200 1]);
cm2 = flipud(cm2);
imwrite(cm2,cb_path);

% WRITE
% replace values and path to colorbar
fi = fopen(combined_svg_path,'r');
fo = fopen([combined_svg_path '.tmp'],'w');

% go read through original SVG
while ~feof(fi)
    tline = fgets(fi);
    % real colorbar path
    tline = strrep(tline, 'colorbarpath', cb_path);
    
    % min and max values
    tline = strrep(tline, 'minvalue', sprintf('%4.2f', values_min));
    tline = strrep(tline, 'maxvalue', sprintf('%4.2f', values_max));
    
    % scaling
    tline = strrep(tline, 'height="1756mm"', ...
        sprintf('height="%fmm"', scaling*1756));
    tline = strrep(tline, 'width="4224.5mm"', ...
        sprintf('width="%fmm"', scaling*4224.5));    
    
    fwrite(fo,tline);
    
end

fclose(fo);
fclose(fi);
movefile([combined_svg_path '.tmp'], combined_svg_path);

% VIEWER
if viewer
    web(combined_svg_path, '-new');
    
    if exist('tmpDir', 'var')
        pause(2);
        rmdir(tmpDir, 's');
    end
end

end


function makeSurface(inname, outname, ID_list, coloring)

% initialize
fi = fopen(inname,'r');
fo = fopen(outname,'w');
tline = fgets(fi);

% go read through original SVG
while ischar(tline)
    
    % find out which element is on the line we are reading
    ID_selected = ~cellfun(@(x) isempty(strfind(tline, x)), ...
        ID_list(:,1));
    %         strcat('"', ID_list(:,1), '"'));
    
    if any(ID_selected)
        if ~isempty(coloring(ID_selected))
            
            % If object is a node:
            position = strfind(tline, 'fill=');
            tline = [tline(1:position+5) coloring{ID_selected}, ...
                tline(position+13:end)];
            
        end
        
        fwrite(fo, tline);
        
    else
        
        % if no specified object is on the line, then write normal line.
        fwrite(fo,tline);
        
    end
    
    tline = fgets(fi);
    
end

% close
fclose(fo);
fclose(fi);

end
