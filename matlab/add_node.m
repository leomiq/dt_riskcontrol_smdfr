% add_node.m
function tree = add_node(tree, parent_node, value, compl_value)
    new_node = struct('left', parent_node, 'right', [], 'id', value,'compl_branch_event',compl_value);
    % updating parent node
    for i = 1:length(tree)
        current_node = tree{i};
        if current_node.id == parent_node
            child_nodes_parent = current_node.right;
            child_nodes_parent(end+1) = value;
            current_node.right = child_nodes_parent;
            tree{i} = current_node;
        end
        % including new node
    end
    tree{end+1} = new_node;
end