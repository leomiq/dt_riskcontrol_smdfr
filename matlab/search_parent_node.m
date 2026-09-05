% search_parent_node.m
function parent_node = search_parent_node(tree,value)
    for i = 1:length(tree)
        current_node = tree{i};
        if current_node.id == value
            parent_node = current_node.left;
        end
    end
end