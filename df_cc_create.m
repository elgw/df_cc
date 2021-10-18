function varargout = df_cc_create(varargin)
% Create correction coefficients for chromatic aberrations
%
% Usage:
% df_cc_create('keyword', value, 'keyword', value, ...)
%
% Keywords:
%  'dots' 
%    A 2D struct array where D{a,b} are dots in channel a vs channel b.
%    each entry D{a,b} should be a Nx6 matrix where the first three columns
%    are (x,y,z) from channel a, and the last three columns are (x,y,z) from
%    channel b.
%  'channels'
%    A cell array with the names of the channels, e.g., {'a594', 'tmr'}.
%  'settings'
%    A struct, use df_cc_create('getDefaults') to get the default values.
%    polyorder: Polygon order, 2 or 3
%    pixelsize: pixel size in nm for x, y, z
%  'getDefaults'
%
%  see df_cc_demo1 etc for example usage

s.polyorder = 2;

for kk = 1:numel(varargin)
    if strcmpi(varargin{kk}, 'dots')
        D = varargin{kk+1};
    end
    if strcmpi(varargin{kk}, 'channels')
        channels = varargin{kk+1};
    end
    if strcmpi(varargin{kk}, 'settings')
        s = varargin{kk+1};
    end
    if strcmpi(varargin{kk}, 'getDefaults')
        % Return default settings
        s.filename = [];
        s.pixelsize = [100, 100, 200];
        s.plot = 0;
        varargout{1} = s;
        return
    end
end

if(s.pixelsize(1) ~= s.pixelsize(2))
    error('Pixel size in x and y has to be the same')
end
assert(size(D,1) == size(D,2));
assert(numel(channels) == size(D,1));
assert(s.polyorder > 0);
assert(s.polyorder < 4);

Cx = cell(size(D,1));
Cy = cell(size(D,1));
dz = cell(size(D,1));
E0 = zeros(size(D,1)); % MSE errors before
E = zeros(size(D,1)); % MSE errors after fitting
E3 = zeros(size(D,1)); % MSE errors after fitting (3D)
N = zeros(size(D,1)); % Number of dots


%% channel aa vs channel bb
for aa = 1:size(D,1)
    % Create matrix for the polynomial model from D{aa}
                
    for bb=1:size(D,1)        
        if aa ~= bb
                        
        Daa = D{aa,bb}(:,1:3);
        Dbb = D{aa,bb}(:,4:6);                
        
        N(aa,bb) = size(D{aa,bb}, 1);
        
        fprintf('%d dots in %s, %d dots in %s\n',...
            size(Daa,1), channels{aa}, ...
            size(Dbb,1), channels{bb});
        
        if size(D{aa,bb},1) < 10
            warning('Using identity transformation. Error will be set to nan.');
            Daa  = 1024*rand(100, 3);
            Dbb = Daa;
        end
            
        % Possibly do a sub-selection on the dots
        % i.e., skip row kk if  numel(sum([Daa(kk,1:3), Dbb(kk,1:3)] ==
        % nan)> 0
        
        MXY1 = df_cc_poly2mat(Daa(:,1:2), s.polyorder);
        
        % Polynomial coefficients        
        Cx{aa,bb} = MXY1\Dbb(:,1);
        Cy{aa,bb} = MXY1\Dbb(:,2);
        dz{aa,bb} = mean(Dbb(:,3)) - mean(Daa(:,3));
                
        Ft = zeros(size(MXY1,1), 2);
        Ft(:,1)=MXY1*Cx{aa,bb};
        Ft(:,2)=MXY1*Cy{aa,bb};
        
        aniso = s.pixelsize(3)/s.pixelsize(1);
        Ft(:,3) = (Daa(:,3)+dz{aa,bb})*aniso;
        
        D2 = df_cc_eudist(Ft(:,1:2), Dbb(:,1:2));
        D3 = df_cc_eudist(Ft(:,1:3), [Dbb(:,1:2), Dbb(:,3)*aniso]);
        D20 = df_cc_eudist(Daa(:,1:2), Dbb(:,1:2));
        D30 = df_cc_eudist([Daa(:, 1:2), aniso*Daa(:,3)], [Dbb(:,1:2), Dbb(:,3)*aniso]);
        
        if s.plot
           error_before_2 = D20;
           error_after_xy = D2;
           error_before_xyz = D30;
           error_after_xyz = D3;
           
            figure,             
            subplot(2,2,1), histogram(s.pixelsize(1)*error_before_2)
            title('Before. XY')
            xlabel('distance [nm]')
            ylabel('# bead pairs')
            legend({sprintf('mean: %1.f nm', mean(s.pixelsize(1)*error_before_2) )});
            
            subplot(2,2,2), histogram(s.pixelsize(1)*error_after_xy)
            title('After. XY')
            xlabel('distance [nm]')
            ylabel('# bead pairs')
            legend({sprintf('mean: %1.f nm', mean(s.pixelsize(1)*error_after_xy) )});
            
            subplot(2,2,3), 
            histogram(s.pixelsize(1)*error_before_xyz)
            title('Before. XYZ')
            xlabel('distance [nm]')
            ylabel('# bead pairs')
            legend({sprintf('mean: %1.f nm', mean(error_before_xyz*s.pixelsize(1)) )});
            
            subplot(2,2,4), 
            histogram(s.pixelsize(1)*error_after_xyz)
            title('After. XYZ')
            xlabel('distance [nm]')
            ylabel('# bead pairs')
            legend({sprintf('mean: %1.f nm', mean(error_after_xyz*s.pixelsize(1)) )});
            dprintpdf(sprintf('%s_vs_%s', channels{aa}, channels{bb}), '--publish')
        end
        
        if size(D{aa,bb},1) < 10
            E(aa,bb) = nan;
            E0(aa,bb) = nan;
            E3(aa,bb) = nan;
        else
            E(aa,bb) = mean(D2.^2);                
            E3(aa,bb) = mean(D3.^2);                
            E0(aa,bb) = mean(D20.^2);
        end
            
        end
    end
end

M.creationDate = datestr(now(), 'YYmmDD');
if exist('df_version')
    M.dotterVersion = df_version();
else
    M.dotterVersion = 'standalone';
end

if isfield(s, 'filename')
    if numel(s.filename) > 0
        fprintf('Writing to %s\n', s.filename)
        M.pixelsize = s.pixelsize;         
        save(s.filename, 'Cx', 'Cy', 'dz', 'E', 'channels', 'M', 'E0', 'N', 'E3', '-v7.3');
        varargout{1} = s.filename;
    end
else
    varargout{1} = '';
end
    
end