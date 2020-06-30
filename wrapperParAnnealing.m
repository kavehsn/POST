% this is a wrapper for a .m file to use with SLURM on Hamilton
% read the index of the data set to be processed from standard input and store in the variable data_set_index

data_set_index = input('');

% print it as a reminder (can be commented with %)

fprintf('\nIndex of data set is:  %d \n', data_set_index);

% compose the name of the inputfile(s) from a base name ,e.g. 'input' 
% and the index

input_file = strcat('SimAnneal_', int2str(data_set_index), '.m');

% and print it

fprintf( 'Input file for data set is: %s \n' , input_file );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here the code or call of a function should be included to import the data set from the file input_file 

run(input_file);