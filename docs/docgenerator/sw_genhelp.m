function doctree = sw_genhelp(varargin)
% generates markdown files from function help
%
% SW_GENHELP('option1', value1, ...)
%

inpForm.fname  = {'sourcepath' 'outpath' 'sidebar'    'fun'      'verstr' 'recalc' 'done'     'docpath'};
inpForm.defval = {{}           ''        'sw_sidebar' zeros(1,0) struct() true     zeros(1,0) ''       };
inpForm.size   = {[1 -1]       [1 -5]    [1 -2]       [1 -3]     [1 1]    [1 1]    [1 -4]     [1 -6]   };

param = sw_readparam(inpForm, varargin{:});

if ~iscell(param.sourcepath)
    path0 = {param.sourcepath};
else
    path0 = param.sourcepath;
end

nPath = numel(path0);

%docroot  = [sw_rootdir 'docs' filesep];
docroot  = [param.outpath filesep];
doctree  = struct('name',cell(1,nPath),'folder',[],'content',[]);

% loop over all path to generate help files
for ii = 1:nPath
    % name of the parent folder
    [~,pp1,pp2] = fileparts(path0{ii});
    doctree(ii).fullname = [pp1 pp2];
    doctree(ii).folder   = path0{ii};
    
    % type of folder
    doctree(ii).isPackage = false;
    doctree(ii).isClass   = false;
    
    switch doctree(ii).fullname(1)
        case '+'
            % package
            doctree(ii).isPackage = true;
            doctree(ii).name      = doctree(ii).fullname(2:end);
        case '@'
            % class
            doctree(ii).isClass = true;
            doctree(ii).name = doctree(ii).fullname(2:end);
        otherwise
            doctree(ii).name = doctree(ii).fullname;
    end
    
    if doctree(ii).isClass
        % class
        name        = doctree(ii).name;
        
        propNames   = properties(name);
        methodNames = methods(name);
        
        funList = [name; cellfun(@(C)[name '.' C],[propNames;methodNames],'UniformOutput',false)];
        doctree(ii).content = struct('fun',cell(1,numel(funList)),'isProp',[],'file',[]);
        [doctree(ii).content(:).fun] = funList{:};
        isProp = num2cell([false true(1,numel(propNames)) false(1,numel(methodNames))]);
        [doctree(ii).content(:).isProp] = isProp{:};
    elseif doctree(ii).isPackage
        % package
        name        = doctree(ii).name;
        
        % find all *.m files in the folder
        fList = dir([path0{ii} filesep '*.m']);
        fList = {fList(:).name};
        % remove .m
        nList = cellfun(@(C)[name '.' C(1:end-2)],fList,'UniformOutput',false);
        doctree(ii).content = struct('file',cell(1,numel(fList)),'fun',[],'isProp',[]);
        [doctree(ii).content(:).file] = fList{:};
        [doctree(ii).content(:).fun]  = nList{:};
    else
        % find all *.m files in the folder
        fList = dir([path0{ii} filesep '*.m']);
        fList = {fList(:).name};
        % remove .m
        nList = cellfun(@(C)C(1:end-2),fList,'UniformOutput',false);
        doctree(ii).content = struct('file',cell(1,numel(fList)),'fun',[],'isProp',[]);
        [doctree(ii).content(:).file] = fList{:};
        [doctree(ii).content(:).fun]  = nList{:};
    end
    
end

if ~isempty(param.fun)
    for ii = 1:nPath
        % keep only the selected function
        doctree(ii).content = doctree(ii).content(ismember({doctree(ii).content.fun},param.fun));
    end
    doctree = doctree(cellfun(@(C)~isempty(C),{doctree.content}));
    nPath   = numel(doctree);
end

for ii = 1:nPath
    % load the help text for each file
    for jj = 1:numel(doctree(ii).content)
        if isempty(doctree(ii).content(ii).file)
            doctree(ii).content(jj).text = strsplit(help(doctree(ii).content(jj).fun),newline);
        else
            doctree(ii).content(jj).text = strsplit(help([doctree(ii).folder filesep doctree(ii).content(jj).file]),newline);
        end
        % remove common leading spaces
        doctree(ii).content(jj).text = sw_rmspace(doctree(ii).content(jj).text);
    end
    
    % find Contents.m file and put it to the first place
    if doctree(ii).isClass
        isContents = ismember({doctree(ii).content(:).fun},doctree(ii).name);
    else
        isContents = ismember({doctree(ii).content(:).file},'Contents.m');
    end
    
    isContents = num2cell(isContents);
    [doctree(ii).content(:).isContents] = isContents{:};
    
    % put Contents.m file to the first place
    [~,idx] = sort(cell2mat(isContents),'descend');
    doctree(ii).content = doctree(ii).content(idx);
end

% convert Contents.m files and mine out titles
for ii = 1:nPath
    cIdx = find([doctree(ii).content(:).isContents]);
    for jj = 1:numel(cIdx)
        [doctree(ii).content(cIdx(jj)), tStruct] = helpcontentfun(doctree(ii).content(cIdx(jj)),doctree(ii));
    end
    % assign title from tList
    fList  = {doctree(ii).content(:).fun};
    
    if all(cellfun(@(C)isempty(C),{tStruct.title}))
        % no titles are defined
        if doctree(ii).isClass
            title0 = 'Methods';
        else
            title0 = 'Files';
        end
        title0 = repmat({title0},1,numel(tStruct));
        [tStruct(:).title] = title0{:};
    end
    
    doctree(ii).utitle = ['Properties' unique({tStruct.title},'stable')];
    % add properties
    if doctree(ii).isClass && any(~[doctree(ii).content.isProp])
        pList = {doctree(ii).content([doctree(ii).content.isProp]).fun};
        [tStruct(end+(1:numel(pList))).fun] = pList{:};
        title0 = repmat({'Properties'},1,numel(pList));
        [tStruct(end-(numel(pList):-1:1)+1).title] = title0{:};
    end
    
    tfList = {tStruct(:).fun};
    
    for jj = 1:numel(fList)
        idx = find(ismember(tfList,fList{jj}));
        if numel(idx) == 1
            doctree(ii).content(jj).title = tStruct(idx).title;
        else
            doctree(ii).content(jj).title = 'Miscellaneous';
            if ~doctree(ii).content(jj).isContents && ~strcmp(doctree(ii).content(jj).fun,doctree(ii).name)
                warning([doctree(ii).content(jj).fun ' has no corresponding title defined, please append Contents.m or class description!'])
            end
        end
    end
    
end

for ii = 1:nPath
    doctree(ii).content(1).frontmatter = struct;
    
    for jj = 1:numel(doctree(ii).content)
        content = doctree(ii).content(jj);
        
        % title
        if content.isContents
            if doctree(ii).isPackage
                content.frontmatter.title = [doctree(ii).name ' package'];
            elseif doctree(ii).isClass
                content.frontmatter.title = [doctree(ii).name ' class'];
            else
                content.frontmatter.title = ['Functions in ' doctree(ii).name];
            end
            content.frontmatter.link = content.frontmatter.title;
        else
            if ~doctree(ii).isPackage && ~doctree(ii).isClass
                %content.frontmatter.title = [content.fun '( )'];
                content.frontmatter.title = content.fun;
            elseif doctree(ii).isClass && content.isProp
                content.frontmatter.title = [content.fun ' property'];
            elseif doctree(ii).isClass
                content.frontmatter.title = [content.fun ' method'];
            else
                content.frontmatter.title = content.fun;
            end
            content.frontmatter.link = content.fun;
        end
        % summary
        content.frontmatter.summary = strtrim(content.text{1});
        % keywords
        content.frontmatter.keywords  = 'sample';
        % sidebar
        content.frontmatter.sidebar   = param.sidebar;
        % permalink
        if content.isContents
            %content.frontmatter.permalink = [doctree(ii).name '.html'];
            content.frontmatter.permalink = doctree(ii).name;
        else
            %content.frontmatter.permalink = [strrep(content.fun,'.','_'),'.html'];
            content.frontmatter.permalink = strrep(content.fun,'.','_');
        end
        % folder
        content.frontmatter.folder    = doctree(ii).name;
        % mathjax
        isMath    = ~isempty(regexp(content.text,'\$\$','once'));
        falsetrue = {'false' 'true'};
        content.frontmatter.mathjax   = falsetrue{isMath+1};
        doctree(ii).content(jj)       = content;
    end
end

% remove H1 line from all help files and last 2 lines from class files
for ii = 1:nPath
    %temp = cellfun(@(C)C(2:end),{doctree(ii).content(~[doctree(ii).content.isContents]).text},'UniformOutput',false);
    %[doctree(ii).content(~[doctree(ii).content.isContents]).text] = temp{:};
    temp = cellfun(@(C)C(2:end),{doctree(ii).content.text},'UniformOutput',false);
    [doctree(ii).content.text] = temp{:};
    
    temp = cellfun(@(C)C(1:end-2),{doctree(ii).content([doctree(ii).content.isContents] & doctree(ii).isClass).text},'UniformOutput',false);
    [doctree(ii).content([doctree(ii).content.isContents] & doctree(ii).isClass).text] = temp{:};
    
end

% YAML java
javaaddpath(YAML.jarfile);
% Load yaml into java obj
snakeyaml = org.yaml.snakeyaml.Yaml;


% add doc to the beginning of doctree
doctree = doctree([1 1:end]);
doctree(1).name      = 'documentation';
doctree(1).folder    = [param.docpath filesep 'source'];
doctree(1).fullname  = 'docs';
doctree(1).isPackage = false;
doctree(1).isClass   = false;
doctree(1).utitle    = {};
% remove existing content
doctree(1).content   = doctree(1).content([]);

% load all the documentation stored in .md files
docFiles = dir([param.docpath filesep 'source' filesep '*.md']);
% remove index.md from the list
docFiles = docFiles(~ismember({docFiles(:).name},'index.md'));

for ii = 1:numel(docFiles)
    docText = strsplit(fileread([param.docpath filesep 'source' filesep docFiles(ii).name]),'---');
    % parse the frontmatter
    if numel(docText) == 3
        docFrontMatter = docText{2};
        docText        = docText{3};
        
        docFrontMatter = YAML.load(docFrontMatter);
    elseif numel(docText) == 1
        docFrontMatter = struct;
        docText        = docText{1};
        docFrontMatter.title = 'NO TITLE DEFINED';
    else
        error('sw_genhelp:ParseDocFrontMatter','Cannot parse the documentation frontmatter!')
    end
    % change permalink, remove .html
    docText = regexprep(docText,[newline newline] ,[newline ' ' newline]);
    doctree(1).content(ii).text   = strsplit(docText,newline);
    [~,docFunName] = fileparts(docFiles(ii).name);
    % remove leading numbers: "DD_...."
    docFunName = regexprep(docFunName,'^\d\d\_','');
    docFrontMatter.permalink = docFunName;
    doctree(1).content(ii).fun    = docFunName;
    doctree(1).content(ii).isProp = false;
    doctree(1).content(ii).file   = docFiles(ii).name;
    doctree(1).content(ii).frontmatter = docFrontMatter;
end

isDoc = num2cell([true false(1,nPath)]);
[doctree(:).isDoc] = deal(isDoc{:});

nPath = nPath+1;

%doctree(:).is
% remove old folders and make new ones
for ii = 1:nPath
    try
        rmdir([docroot 'pages' filesep doctree(ii).name],'s')
    catch
    end
    % add new empty folder
    mkdir([docroot 'pages' filesep doctree(ii).name]);
end

% copy index.md straight into the doc folder
if ~isempty(dir([param.docpath filesep 'source' filesep 'index.md']))
    copyfile([param.docpath filesep 'source' filesep 'index.md'],[docroot 'index.md']);
end

% copy resources
copyfile([param.docpath filesep 'resources' filesep '*'],[docroot 'images' filesep]);

% run the examples
imgPath = [docroot 'images' filesep 'generated'];
if param.recalc
    if ~isempty(dir(imgPath))
        rmdir(imgPath,'s')
    end
    mkdir(imgPath)
end
% close all figures open before the examples
close('all');

for ii = 1:numel(doctree)
    content = doctree(ii).content;
    for jj = 1:numel(content)
        % find ``` code blocks
        text0    = content(jj).text;
        blockIdx = find(cellfun(@(C)~isempty(C),regexp(text0,'```')));
        if mod(numel(blockIdx),2)
            error('Wrong number of block comments in %s',content(jj).fun)
        end
        
        blockIdx = reshape(blockIdx',2,[]);
        for kk = size(blockIdx,2):-1:1
            % code
            if ~isempty(regexp(text0{blockIdx(1,kk)},'none','ONCE'))
                text0{blockIdx(1,kk)} = regexprep(text0{blockIdx(1,kk)},'none','');
            else
                text0{blockIdx(1,kk)} = [text0{blockIdx(1,kk)} 'matlab'];
            end
            selBlock = text0((blockIdx(1,kk)+1):(blockIdx(2,kk)-1));
            if any(cellfun(@(C)~isempty(C),regexp(text0,'^>>')))
                %
                imgName = ['generated/' content(jj).frontmatter.permalink(1:(end-5))];
                [newText,~] = sw_example(selBlock,[docroot 'images' filesep],imgName,param.recalc,content(jj).fun);
                if ~isempty(regexp(newText{end},'```','once'))
                    text0 = [text0(1:blockIdx(1,kk)) newText(1:(end-1))' text0((blockIdx(2,kk)+1):end)];
                else
                    text0 = [text0(1:blockIdx(1,kk)) newText' text0(blockIdx(2,kk):end)];
                end
            end
        end
        content(jj).text = text0;
        %if any(cellfun(@(C)~isempty(C),regexp(text0,'>>')))
        %    imgName = ['generated/' content(jj).frontmatter.permalink(1:(end-5))];
        %    [content(jj).text,~] = sw_example(text0,[docroot 'images' filesep],imgName);
        %end
    end
    doctree(ii).content = content;
end

% get all the help text
content = [doctree.content];
allhelp = {content.text};

for ii = 1:numel(allhelp)
    allhelp{ii} = cellfun(@(C)[C newline],allhelp{ii},'UniformOutput',false);
    allhelp{ii} = [allhelp{ii}{:}];
end

% % temporary convert
% no need any more
% for ii = 1:numel(allhelp)
%     if isempty(regexp(allhelp{ii},'###','once'))
%         allhelp{ii} = sw_convhelp(allhelp{ii},false);
%     end
% end

% generate all links
pLink = cellfun(@(C)C.permalink,{content.frontmatter},'UniformOutput',false);
fun   = {content.fun};
% replace function names with links
for ii = 1:numel(pLink)
    allhelp = regexprep(allhelp,['\[(' fun{ii} ')\]'],['[$1](' pLink{ii} ')']);
    %allhelp(~isCont) = regexprep(allhelp(~isCont),['(\W+?)(@?' fun{ii} '(\(\))?)(\W+)'],['$1[$2](' pLink{ii} ')$3']);
    %allhelp = regexprep(allhelp,['\n([^\<^\<].*?\W+?)(@?' fun{ii} '(\(\))?)(\W+)'],['$1[$2](' pLink{ii} ')$3']);
end
% exchange $ --> $$ for math, make it easier to write MarkDown
allhelp = regexprep(allhelp,'\$','$$');
% exchange text into symbols, e.g. \\Angstrom --> A
%sText = {'Angstrom' 'hbar' 'alpha' 'beta' 'gamma' 'degree' 'sigma' 'deg'};
%cText = {'ang' 'hbar' 'alpha' 'beta' 'gamma' 'deg' 'sigma' 'deg'};
%for ii = 1:numel(sText)
%    allhelp = regexprep(allhelp,['\\\\' sText{ii}],symbol(cText{ii}));
%end

% exchange any symbol referenced by \\label into the output of the symbol
% command symbol(label)
allhelp = regexprep(allhelp,'\\\\(\w+)','${symbol($1,true)}');
% substitue [matlab.funname] links into
% "[funname](https://www.mathworks.com/help/matlab/ref/funname.html)"
allhelp = regexprep(allhelp,'\[matlab\.(\w+?)\]','[$1](https://www.mathworks.com/help/matlab/ref/$1.html)');

allhelp = regexprep(allhelp,'{{(\w+?) ','{% include $1.html content=" ');
allhelp = regexprep(allhelp,'}}','" %}');

% restore standard KramDown behaviour by converting standard KramDown into
% {include } like shit
% image with caption and url
regexp0 = '\[\!\[([^\!\[\]]+?)\]\((\S+?)\)\{(.+?)\}\]\((\S+?)\)';
allhelp = regexprep(allhelp,regexp0,'{% include image.html file="$2" alt="$1" caption="$3" url="$4" %}');

% image with caption only
regexp0 = '\!\[([^\!\[\]]+?)\]\((\S+?)\)\{(.+?)\}';
allhelp = regexprep(allhelp,regexp0,'{% include image.html file="$2" alt="$1" caption="$3"%}');

% image without caption
regexp0 = '\!\[([^\!\[\]]+?)\]\((\S+?)\)';
allhelp = regexprep(allhelp,regexp0,'{% include image.html file="$2" alt="$1" %}');

% convert section names
allhelp = regexprep(allhelp,'\[section\.(\S+?)\]','\[$1\]($1)');

idx = 1;
% generate the .md files
for ii = 1:nPath
    for jj = 1:numel(doctree(ii).content)
        
        content     = doctree(ii).content(jj);
        %frontmatter = YAML.dump(content.frontmatter)
        
        % Convert to matlab object
        frontmatter = char(snakeyaml.dump(YAML.dump_data(content.frontmatter)));
        % save the text as .md file
        if isfield(content.frontmatter,'folder')
            newFile = [docroot 'pages' filesep content.frontmatter.folder filesep content.fun '.md'];
        else
            newFile = [docroot 'pages' filesep content.fun '.md'];
        end
        fid = fopen(newFile,'w');
        % add newline
        %helpText = cellfun(@(C)[C newline],content.text,'UniformOutput',false);
        %        help1 = allhelp{idx};
        
        %helpText = cellfun(@(C)[C newline],help1,'UniformOutput',false);
        
        %fprintf(fid,'%s',['---' newline frontmatter newline '---' newline helpText{:}]);
        fprintf(fid,'%s',['---' newline frontmatter newline '---' newline allhelp{idx} newline '{% include links.html %}' newline]);
        fclose(fid);
        idx = idx + 1;
    end
end

% generate sidebar YAML file
if  all(isfield(param.verstr,{'Version' 'Revision'}))
    if isempty(param.verstr.Version)
        verStr = ['R' param.verstr.Revision];
    else
        verStr = param.verstr.Version;
    end
else
    verStr = '';
end

sidebar = struct;
sidebar.entries.title   = 'Sidebar';
sidebar.entries.product = 'SpinW';
sidebar.entries.version = verStr;
% documentation
%sidebar.entries.folders(1).title  = 'Documentation';
%sidebar.entries.folders(1).output = 'web, pdf';

if isempty(param.done)
    done = false(1,numel(doctree));
else
    done = param.done;
end

okSymbol = symbol('ok');

for ii = 1:nPath
    isDoc = doctree(ii).isDoc;
    
    if done(ii)
        sidebar.entries.folders(ii).title  = [doctree(ii).name ' ' okSymbol];
    else
        sidebar.entries.folders(ii).title  = doctree(ii).name;
    end
    
    if isDoc
        % first letter upp
        tempStr = sidebar.entries.folders(ii).title;
        sidebar.entries.folders(ii).title  = [upper(tempStr(1)) tempStr(2:end)];
    end
    
    sidebar.entries.folders(ii).output = 'web, pdf';
    
    if ~isDoc
        sidebar.entries.folders(ii).folderitems(1).title  = 'Description';
        %sidebar.entries.folders(ii).folderitems(1).url    = ['/' doctree(ii).name '.html'];
        sidebar.entries.folders(ii).folderitems(1).url    = ['/' doctree(ii).name];
        sidebar.entries.folders(ii).folderitems(1).output = 'web, pdf';
        
        % find unique titles
        tCell   = {doctree(ii).content.title};
        %tUnique = unique(tCell,'stable');
        tUnique = doctree(ii).utitle;
        mIdx = find(ismember(tUnique,'Miscellaneous'));
        if ~isempty(mIdx)
            tUnique = tUnique([1:(mIdx-1) (mIdx+1):end mIdx]);
        end
        
        for jj = 1:numel(tUnique)
            idx = find(ismember(tCell,tUnique{jj}) & ~[doctree(ii).content.isContents]);
            if ~isempty(idx)
                % title
                sidebar.entries.folders(ii).folderitems(1).subfolders(jj).title = tUnique{jj};
                sidebar.entries.folders(ii).folderitems(1).subfolders(jj).output = 'web, pdf';
                
                for kk = 1:numel(idx)
                    content = doctree(ii).content(idx(kk));
                    sidebar.entries.folders(ii).folderitems(1).subfolders(jj).subfolderitems(kk).title  = content.frontmatter.link;
                    sidebar.entries.folders(ii).folderitems(1).subfolders(jj).subfolderitems(kk).url    = ['/' content.frontmatter.permalink];
                    sidebar.entries.folders(ii).folderitems(1).subfolders(jj).subfolderitems(kk).output = 'web, pdf';
                end
            end
        end
    else
        for kk = 1:numel(doctree(ii).content)
            content = doctree(ii).content(kk);
            sidebar.entries.folders(ii).folderitems(kk).title  = content.frontmatter.title;
            sidebar.entries.folders(ii).folderitems(kk).url    = ['/' content.frontmatter.permalink];
            sidebar.entries.folders(ii).folderitems(kk).output = 'web, pdf';
        end

        
    end
    
end

%yamlStr = YAML.dump(sidebar);
yamlStr  = char(snakeyaml.dump(YAML.dump_data(sidebar)));
yamlStr1 = strsplit(yamlStr,newline);

% add extra newline where {} unit is inline with some previous text
bLine = regexp(yamlStr1,': {');
bLine2 = find(cellfun(@(C)~isempty(C),bLine));
bLine  = bLine(bLine2);
for ii = 1:numel(bLine2)
    newLine = yamlStr1{bLine2(ii)}(bLine{ii}+2:end);
    yamlStr1{bLine2(ii)} = yamlStr1{bLine2(ii)}(1:bLine{ii});
    yamlStr1 = yamlStr1([1:bLine2(ii) bLine2(ii) (bLine2(ii)+1):end]);
    nSpace = sum(yamlStr1{bLine2(ii)}==' ');
    yamlStr1{bLine2(ii)+1} = [repmat(' ',1,nSpace) '- ' newLine];
    bLine2(ii+1:end) = bLine2(ii+1:end)+1;
end

% add extra '-' that is required by Jekyll
nSpace   = zeros(1,numel(yamlStr1));
for ii = 1:numel(yamlStr1)
    temp = find(diff(yamlStr1{ii}~=' '),1,'first');
    if isempty(temp)
        nSpace(ii) = 0;
    else
        nSpace(ii) = temp;
    end
end

iLine  = find(diff(nSpace))+1;
nSpace = nSpace(iLine);

iLine(nSpace==0)  = [];
nSpace(nSpace==0) = [];

for ii = 1:numel(iLine)
    if yamlStr1{iLine(ii)}(nSpace(ii)+1) ~= '-' && yamlStr1{iLine(ii)-1}(nSpace(ii)-1) ~= '-'
        yamlStr1{iLine(ii)}(nSpace(ii)-1) = '-';
    end
end
yamlStr1 = cellfun(@(C)[C newline],yamlStr1,'uniformoutput',false);
yamlStr1 = [yamlStr1{:}];

fid = fopen([docroot filesep '_data' filesep 'sidebars' filesep 'sw_sidebar.yml'],'w');
fprintf(fid,yamlStr1);
fclose(fid);

end

function [doccontent, tList] = helpcontentfun(doccontent,doctree)
% convert Contents.m / class help text
fprintf('%s\n',doccontent.fun);
text = doccontent.text;

header = {'### Files' '### Methods'};
% find line "Files" / "Methods"
lIdx  = ismember(doccontent.text,header);
lIdx  = find(lIdx);
tList = struct('fun',{},'title',{});

% generate H1 lines
H1 = cellfun(@(C)C(1),{doctree.content.text});
% generate all functions in the folder
fList = {doctree.content.fun};

if numel(lIdx) == 1
    endoftoc = false;
    
    title = '';
    idx = lIdx+1;
    
    while idx<=numel(text) && ~endoftoc
        if isempty(strtrim(text{idx}))
        elseif text{idx}(1)=='#' && text{idx}(4)~='#'
            % line starting with ###
            endoftoc = true;
        elseif text{idx}(1)=='#' && strcmp(text{idx}(1:4),'####')
            % sub header line
            title = strtrim(text{idx}(5:end));
        else
            % function reference
            fun = strtrim(text{idx});
            % add H1 line
            fIdx = find(ismember(fList,fun),1);
            if ~isempty(fIdx)
                text{idx} = ['* [' fun '] ' H1{fIdx}];
            end
            tList(end+1).fun = fun; %#ok<AGROW>
            tList(end).title = title;
        end
        idx = idx+1;
    end
    
    doccontent.text = text;
else
    error('sw_genhelp:MissingContent','Contents.m file is missing for %s!',name);
end

end