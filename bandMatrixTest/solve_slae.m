input_file = getenv('INPUT_MATRIX_FILE');

if isempty(input_file)
    error('Error: Environment variable INPUT_MATRIX_FILE is not set.');
end

fid = fopen(input_file, 'r');
if fid == -1
    error('Error: Cannot open input file %s', input_file);
end

n = fscanf(fid, '%lf', 1); 
fscanf(fid, '%lf', 1); 
A = fscanf(fid, '%lf', [n, n])'; 
b = fscanf(fid, '%lf', [1, n])'; 
fclose(fid); 

x = A \ b; 

fid = fopen('msolution.txt', 'w');
fprintf(fid, '%.15f ', x); 
fclose(fid);