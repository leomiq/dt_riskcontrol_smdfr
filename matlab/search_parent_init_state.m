% search_parent_init_state.m
function parent_init_state = search_parent_init_state(tree_data,parent_node)
    for i = 1:length(tree_data)
        branch_id = tree_data{i}.id;
        if branch_id == parent_node
            parent_init_state = tree_data{i}.init_state;
            break
        end
    end
end