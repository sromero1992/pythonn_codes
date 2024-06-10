% Set the Python environment
pe = pyenv('Version', 'F:\Anaconda\envs\scanpy_env\python.exe');

% Check if the environment is loaded
if pe.Status == "NotLoaded"
    % Load the environment by executing a simple Python command
    py.exec('import sys');
end

% Display the environment details
disp(pyenv);



% Example adjacency matrix (replace with your data)
adj_matrix = [0 1 1; 1 0 1; 1 1 0];


% Save the adjacency matrix to a text file
save('adj_matrix.txt', 'adj_matrix', '-ascii');

% Path to the Python executable and the script
python_executable = 'F:\Anaconda\envs\scanpy_env\python.exe';  % Use 'python3' if necessary
python_script = 'run_leiden.py';

% Call the Python script with the adjacency matrix file as argument
system_command = sprintf('%s %s adj_matrix.txt', python_executable, python_script);
system(system_command);

% Load the clustering results
clusters = jsondecode(fileread('clusters.json'));

% Display the clustering results
disp('Leiden clustering results:');
disp(clusters);


