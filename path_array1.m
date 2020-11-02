function [path_arr] = path_array1(path,EE,n)
    path_arr = zeros(EE,1);Ap = zeros(n,n);
    for i=1:n
        for j=2:length(path)
            curr_node = path(j);
            prev_node = path(j-1);
            Ap(curr_node,prev_node)=1;
        end
    end
    path_arr=reshape(Ap',[EE,1]);
end