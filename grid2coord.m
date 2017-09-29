function coord_path = grid2coord(path, grid_center_coord)
    [path_size, num_dim] = size(path);
    coord_path = zeros(path_size,num_dim);
    for i = 1:1:path_size
        if num_dim == 2
            coord_path(i,:) = grid_center_coord{path(i,1), path(i,2)}; 
        else
            coord_path(i,:) = grid_center_coord{path(i,1), path(i,2), path(i,3)}; 
        end
    end
end