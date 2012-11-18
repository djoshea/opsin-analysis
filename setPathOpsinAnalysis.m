function setPathOpsinAnalysis()
    % requires matlab-utils to be on the path already
    pathRoot = pathToThisFile(); 
    fprintf('Path: Adding opsin-analysis at %s\n', pathRoot);
    addPathRecursive(pathRoot);
end
