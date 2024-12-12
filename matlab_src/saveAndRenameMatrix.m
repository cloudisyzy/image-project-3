function renamedMatrix = saveAndRenameMatrix(matrix, newName, fileName)
    save(fileName, 'matrix');
    loadedData = load(fileName);
    originalName = 'matrix';
    if isfield(loadedData, originalName)
        tempMatrix = loadedData.(originalName);
    else
        error('error');
    end
    eval(sprintf('%s = tempMatrix;', newName));
    save(fileName, newName, '-append');
    renamedMatrix = eval(newName);
end
