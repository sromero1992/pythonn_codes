function  sce = leiden_annotation(sce, method, species)
    % leiden_annotation computes leiden clustering interfaced from
    % python3.11 with Mutual Nearest Neighbors (MNN) or K-Nearest
    % Neighbors (KNN).
    % INPUT:
    % sce -------> SCE object 
    % method ----> Method  to find neighbors (mnn or knn)
    % OUTPUT:
    % sce ------> sce object containing Leiden clusters and corresponding
    %             annotation
    species = lower(species);
    method = lower(method);

    fprintf("WARNING: sce object should be prior QC if desired... \n")
    fprintf("Leiden annotation with method: %s \n", method);

    % Set the Python environment (Python 3.11)
    env_bin = 'C:\Users\ssromerogon\.conda\envs\scanpy_env_311\python.exe';
    pe = pyenv('Version', env_bin);
    
    % Check if the environment is loaded
    if pe.Status == "NotLoaded"
        % Load the environment by executing a simple Python command
        py.exec('import sys');
    end
    
    % Display the environment details
    disp(pyenv);
    
    switch method
        case 'mnn'
            n_neighbors = 50; % 20, 30, 50 (the more neighbors, the less clusters)
            adjX = adj_mat_construct(sce, 'mnn', n_neighbors);
    
        case 'knn'
            n_neighbors = 15; % 10, 15 % the less the neighbors produces more clusters
            adjX = adj_mat_construct(sce, 'knn', n_neighbors);

        otherwise
            error('Unknown method: %s. Method should be either ''mnn'' or ''knn''.', method);
    end

    % Save the adjacency matrix to a text file
    adj_file = 'adjX.txt';
    save(adj_file, 'adjX', '-ascii');
    
    % Path to the Python executable and the script
    python_executable = env_bin;  
    python_script = 'run_leiden.py';
    
    % Call the Python script with the adjacency matrix file as argument
    system_command = sprintf('%s %s adjX.txt', python_executable, python_script);
    [status, cmdout] = system(system_command);

    % Clean up the temporary adjacency matrix file
    if exist(adj_file, 'file')
        delete(adj_file);
    end    

    % Check for errors
    if status ~= 0
        disp('Error running the Python script:');
        disp(cmdout);
    else
        % Load the clustering results
        clusters = jsondecode(fileread('clusters.json'));
    
        delete('clusters.json')
    
        disp("Parsing Leiden:")
        disp(cmdout);
    
        % Display the clustering results
        nclus = length( unique(clusters) );
        fprintf('Number Leiden clusters %d \n', nclus);
    end

    % Assign clustering results to the SCE object
    sce.c_cluster_id = clusters;

    % Embed cells and assign cell types
    tic;
    sce = sce.embedcells('umap3d', true, false, 3);
    %sce = sce.embedcells('umap2d', true, false, 2);
    fprintf("Annotating species %s \n\n", species)
    sce= sce.assigncelltype(species, false);
    time_assign = toc;
    fprintf("Time for cell annotation and embedding: %f \n", time_assign);
end
