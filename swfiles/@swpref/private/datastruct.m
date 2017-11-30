function d = datastruct()
% Function called to create a default preference object.
%
%  {{warning Internal function for the Spin preferences.}}
%
% ### Syntax
%
% 'prefs = datastruct()'
%
% ### Description
%
% 'prefs = datastruct()' creates a structure with the following fields:
%
% 'Name' a cell array giving the name of each dynamic property.
%
% 'Validation' a cell array with functions to evaluate when a
% property is set.
%
% 'Value' a cell array giving the default value of a dynamic property
%
% 'Label' a cell array giving a short description on the dynamic
% property.
%

d.Name = {
    'fid',...
    'expert',...
    'tag',...
    'nmesh',...
    'maxmesh',...
    'npatch',...
    'fontsize',...
    'tid',...
    'colormap',...
    'usemex',...
    'docurl'
    };

Size = {
    [1, 1],...
    [1, 1],...
    [    ],...
    [1, 1],...
    [1, 1],...
    [1, 1],...
    [1, 1],...
    [1, 1],...
    [1, 1],...
    [1, 1],...
    [    ]
    };

d.Validation = {
    {@isnumeric, @mustBeInteger, @mustBeNonnegative, @(x) check_size(x,Size{1}), @(x) mustBeLessThan(x,256)},...
    {@islogical, @(x) check_size(x,Size{2})},...
    {@ischar},...
    {@isnumeric, @mustBeInteger, @mustBeNonnegative, @(x) check_size(x,Size{4}), @(x) mustBeGreaterThanOrEqual(x,1), @(x) mustBeLessThan(x,256)},...
    {@isnumeric, @mustBeInteger, @mustBeNonnegative, @(x) check_size(x,Size{5}), @(x) mustBeGreaterThan(x,1), @(x) mustBeLessThan(x,256)},...
    {@isnumeric, @mustBeInteger, @mustBeNonnegative, @(x) check_size(x,Size{6}), @(x) mustBeGreaterThan(x,1), @(x) mustBeLessThan(x,256)},...
    {@isnumeric, @mustBeInteger, @mustBeNonnegative, @(x) check_size(x,Size{7}), @(x) mustBeGreaterThan(x,4), @(x) mustBeLessThan(x,256)},...
    {@isnumeric, @mustBeInteger, @mustBeNonnegative, @(x) check_size(x,Size{8}), @(x) mustBeLessThan(x,256)},...
    {@(x) check_size(x,Size{9}), @check_mex, @(x) isa(x,'function_handle')},...
    {@(x) check_size(x,Size{10}), @islogical},...
    {@ischar, @(x) strfind(x,'http')}
    };

d.Value = {
    1,...
    false,...
    'swplot',...
    1,...
    6,...
    20,...
    12,...
    1,...
    @cm_inferno,...
    false,...
    'https://tsdev.github.io/spinwdoc'
    };

d.Label =  {
    'file identifier for text output, default value is 1 (Command Window)'...
    'expert mode (1) gives less warnings (not recommended), default value is 0'...
    'defines the tag property of the crystal structure plot figures'...
    'default number of subdivision of the icosahedron for plotting'...
    'maximum number of subdivision of the icosahedron for plotting'...
    'number of edges for patch object'...
    'fontsize for plotting'...
    'identifier how the timer is shown, default value is 1 (Command Window), value 2 gives graphical output'...
    'default colormap'...
    'if true, mex files are used in the spin wave calculation'...
    'url to the documentation server'...
    };

    function out = check_size(obj,S)
        % checks to see if an object is the wrong size.
        %
        %  {{warning Internal function for the Spin preferences.}}
        %
        % ### Syntax
        % 
        % 'logical = check_size(toBeChecked,size)'
        %
        % ### Description
        %
        % 'logical = check_size(toBeChecked,size)' checks to see if an 
        % object 'obj 'is the expected size given by 'size'. An error is
        % thrown if there is a difference.
        %
        
        sz = size(obj);
        if ~all(sz == S)
            error('spref:WrongSize','Value to be asigned is the wrong size [%i, %i] not [%i, %i]',sz(1), sz(2), S(1), S(2))
        else
            out = 1;
        end
    end

    function out = check_mex(~)
        % checks to see if mex files are available.
        %
        %  {{warning Internal function for the Spin preferences.}}
        %
        % ### Syntax
        % 
        % 'logical = check_mex(obj)'
        %
        % ### Description
        %
        % 'logical = check_mex(obj)' checks to see if files 'chol_omp' and 
        % 'eig_omp' are present in the MATLAB path.An error is thrown if 
        % they do not exist.
        %
        
        if ~(exist('chol_omp','file')==3 && exist('eig_omp','file')==3)
            error('spref:MissingMex','Necessary mex files are missing, compile them!')
        else
            out = 1;
        end
    end
end