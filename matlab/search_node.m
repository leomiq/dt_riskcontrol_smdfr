% search_node.m
function node = search_node(tree,value)
    for i = 1:length(tree)
        current_node = tree{i};
        if current_node.id == value
            node = current_node;
        end
    end
end