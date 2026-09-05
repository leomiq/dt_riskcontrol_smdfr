% search_parent_finalOP.m
function parent_finalOP = search_parent_finalOP(tree_data,parent_node)
    for i = 1:length(tree_data)
        branch_id = tree_data{i}.id;
        if branch_id == parent_node
            parent_finalOP = tree_data{i}.final_OP;
            break
        end
    end
end