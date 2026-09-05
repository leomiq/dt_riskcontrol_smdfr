% search_final_state.m
function final_state = search_final_state(tree_data,value,Ncomps)
    node = 0;
    for i = 1:length(tree_data)
        current_node = tree_data{i};
        if current_node.id == value
            node = current_node;
        end
    end
    current_sys_state = node.data.current_sys_state.Data(end,:);
    if all(current_sys_state == ones(Ncomps))
        availability_state = 1;
    else
        availability_state = 0;
    end
    protection_systems = node.data.protection_systems.Data(end);
    final_state = [protection_systems, availability_state];
end