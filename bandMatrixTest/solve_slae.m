input_file = getenv('INPUT_MATRIX_FILE');

if isempty(input_file)
    error('Error: Environment variable INPUT_MATRIX_FILE is not set.');
end

fid = fopen(input_file, 'r');
if fid == -1
    error('Error: Cannot open input file %s', input_file);
end

n = fscanf(fid, '%lf', 1); 
b_size = fscanf(fid, '%lf', 1); 
A = fscanf(fid, '%lf', [n, n])'; 
b = fscanf(fid, '%lf', [1, n])'; 
fclose(fid); 

x = A \ b;

output_filename = sprintf('matlabSolutions/msolution%dx%d.txt', n, n);

fid = fopen(output_filename, 'w');
if fid == -1
    error('Error: Cannot open output file %s', output_filename);
end
fprintf(fid, '%.15f ', x); 
fclose(fid);

fprintf('Solution saved to %s\n', output_filename);