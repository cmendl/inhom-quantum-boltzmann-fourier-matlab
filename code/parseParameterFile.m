function params = parseParameterFile(filename)
%PARSEPARAMETERFILE - Parse simulation parameter file
%
%    params = PARSEPARAMETERFILE(filename)

fid = fopen(filename);
C = textscan(fid,'%s %s %s');
fclose(fid);

if length(C)~=3
	error('Expecting file entries of the form "name = value"');
end

% partition into parameter name, equal sign and parameter value
name  = C{1};
eq    = C{2};
value = C{3};

for j=1:length(name)
	if ~strcmp(eq{j},'=')
		error('Expecting equal sign in the middle');
	end
	switch name{j}
		case {'J','M','L','R','h','numVol','dt','numsteps','lambda'}
			params.(name{j}) = str2double(value{j});
		case {'boundaryType','filenameWinit','filenameWevolv','filenameBext','filenameWMaxwL','filenameWMaxwR'}
			params.(name{j}) = value{j};
		otherwise
			warning('Unrecognized parameter name "%s"',name{j});
	end
end
