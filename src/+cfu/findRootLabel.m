function root = findRootLabel(labels,id)
    root = id;
    while(labels(root)~=root)
        root = labels(root);
    end
end