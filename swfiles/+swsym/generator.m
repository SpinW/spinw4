function [symOp, symInfo] = generator(sym, fid0)
% returns symmetry operators of a given space group
%
% [symOp, symInfo] = SWSYM.GENERATOR(sym, {fid})
%
% It gives the symmetry elements based on the space group number or given
% list of symmetry operators. Without arguments, returns the name of all
% space groups stored in symmetry.dat file.
%
% Input:
%
% sym           Either the label of the space group or the index from
%               the International Tables of Crystallography or string
%               containing the space group operators in the same format as
%               used in the symmetry.dat file.
% fid           For printing the symmetry operators:
%                   0   no printed output (Default)
%                   1   standard output (Command Line)
%                   fid text file opened before with the fid = fopen(path)
%
% Output:
%
% symOp         Matrices defining the symmetry operators, dimensions are 
%               [3 4 nOp].
% symInfo       Structure containing additional information about the space 
%               group with the following fields:
%   name            Name of the space group in string. If function called
%                   with no input, name stores the name of all spase groups
%                   from symmetry.dat in a cell.
%   str             The string of the symmetry operations.
%   num             The index of the symmetry in the symmetry.dat file.
%
% See also SWSYM.ADD, SPINW, SPINW.GENCOUPLING, SWSYM.POSITION.
%

symInfo = struct('name','','str','','num',-1);

if nargin < 2
    fid0 = 0;
end

if nargin > 0
    symNum = -1;
    
    if iscell(sym)
        % handle input where the operators are already given
        symOp = [sym{1} permute(sym{2},[1 3 2])];        
        return
    elseif swsym.isop(sym)
        if ~isempty(sym)
            symOp = sym;
        else
            symOp = [eye(3) zeros(3,1)];
        end        
        return
        
    elseif ~ischar(sym) && numel(sym)~=1
        error('generator:WrongInput','Wrong symmetry definition!');
    end
end

% Open the symmetry definition file.
symPath = [sw_rootdir 'dat_files' filesep 'symmetry.dat'];
fid = fopen(symPath);

if fid == -1    
    error('spinw:sw_gensym:FileNotFound',['Symmetry definition file not found: '...
        regexprep(symPath,'\' , '\\\') '!']);
end

% Just returns the name of all space groups.
if nargin == 0
    ii = 1;
    symName = cell(1,230);
    while ~feof(fid)
        textLine    = fgetl(fid);
        symName{ii} = [textLine(7:17) sprintf(' (%3i)',ii)];
        ii          = ii+1;
    end
    fclose(fid);
    symOp  = zeros(3,4,0);
    symInfo.str  = '';
    symInfo.name = symName;
    symInfo.num  = 1:(ii-1);
    return
    
else
    
    if ischar(sym)
        if isempty(sym)
            symStr  = 'x,y,z';
            symName = 'P 1';
        elseif any(sym=='x') && any(sym=='y') && any(sym=='z') && any(sym==',')
            symStr = sym;
            symName = '';
        else
            % find symmetry label
            symName = sym;
            symName(end+1:11) = 32;
            symIdx = 0;
            ii     = 1;
            while (symIdx == 0) && ~feof(fid)
                textLine    = fgetl(fid);
                if strfind(symName,textLine(7:17))
                    symIdx = ii;
                end
                ii = ii+1;
            end
            fclose(fid);
            if symIdx == 0
                error('generator:WrongInput','Symmetry name does not exists (case insensitive)!');
            end
            symNum = symIdx;
            symStr = textLine(20:end);
        end
        
    else
        symNum = sym;
        
        if symNum<0
            fclose(fid);
            error('generator:WrongInput','Symmetry number has to be positive integer!');
        elseif symNum == 0
            fclose(fid);
            symOp        = [eye(3) zeros(3,1)];
            symInfo.name = 'No sym';
            symInfo.str  = 'x,y,z';
            symInfo.num  = 0;
            return
        end
        ii = 1;
        while (ii<=symNum) && ~feof(fid)
            textLine = fgetl(fid);
            ii = ii+1;
        end
        fclose(fid);
        if ii <= symNum
            error('generator:WrongInput','Symmetry number not found!')
        end
        symStr  = textLine(20:end);
        symName = textLine(7:17);
    end
end

symOp = zeros(3,4,30);
vNew   = zeros(3,1);

nNew  = 1;
nOp   = 1;
nSign = 1;

ii=1;
while(ii<=length(symStr))
    if symStr(ii)==','
        symOp(nNew,1:3,nOp) = vNew;
        vNew  = vNew*0;
        nSign = 1;
        nNew  = mod(nNew,3)+1;
    elseif symStr(ii)==';'
        symOp(nNew,1:3,nOp) = vNew;
        vNew = vNew*0;
        nSign = 1;
        nNew  = 1;
        nOp   = nOp+1;
    elseif symStr(ii)=='x'
        vNew(1) = nSign;
    elseif symStr(ii)=='y'
        vNew(2) = nSign;
    elseif symStr(ii)=='z'
        vNew(3) = nSign;
    elseif symStr(ii)=='-'
        nSign = -1;
    elseif symStr(ii)=='+'
        nSign = 1;
    elseif (symStr(ii)=='1')||(symStr(ii)=='2')||(symStr(ii)=='3')
        symOp(nNew,4,nOp) = (symStr(ii)-'0')/(symStr(ii+2)-'0');
        ii = ii+2;
    end
    ii = ii+1;
    
end

symOp(nNew,1:3,nOp) = vNew;
symOp = symOp(:,:,1:nOp);

% cut trailing spaces from symName
if isnan(symName)
    symName = '';
else
    symName = strtrim(symName);
end

% print symmetry elements
if fid0~=0
    fprintf(fid0, 'Generators of space group: %s\n',symName);
    
    sStr = strtrim(strsplit(swsym.str(symOp),';')');
    for ii = 1:numel(sStr)
        fprintf(fid0,'(%02d) %s\n',ii,sStr{ii});
    end
end

% generate output
symInfo.name = symName;
symInfo.str  = strtrim(symStr);
symInfo.num  = symNum;

end