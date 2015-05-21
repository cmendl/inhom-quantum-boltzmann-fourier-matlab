function data = loadData(fn,type,dim)
%LOADDATA - Load binary data in 'double' format from disk
%
%    data = LOADDATA(fn,type,dim)
%
%    Input:    fn		filename
%              type     'real' or 'complex'
%              dim      data in the file is interpreted as array of dimension 'dim'
%
%    Output:   data     data in file

fid = fopen(fn,'rb');
data = fread(fid,'double');
fclose(fid);

% interleave real and imaginary entries
if strcmpi(type,'complex')
	data = data(1:2:end-1) + 1i*data(2:2:end);
end

data = reshape(data,dim);
