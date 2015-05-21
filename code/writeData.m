function writeData(fn,data)
%WRITEDATA - Write binary data in 'double' format to disk
%
%    WRITEDATA(fn,data)

% row vector
data = data(:).';

% sort interleaved real and imaginary entries into a long array
if ~isreal(data)
	data = reshape([real(data);imag(data)],1,[]);
end

fid = fopen(fn,'wb');
fwrite(fid,data,'double'); 
fclose(fid);
