clc; clear all; close all;
%% parameters
cell_size_row = 0.16;
cell_size_col = 0.16;
grid_size_row = 14;
grid_size_col = 14;
%visualize occupancy grid
grid_left_top     = [0.0 grid_size_row * cell_size_row];
grid_center_coord = cell(grid_size_row, grid_size_col);
POLYORDER = 5;

temp_coord = zeros(2,1);
for i = 1:1:grid_size_row
   for j = 1:1:grid_size_col
       temp_coord(1) = (i-1) * cell_size_row + 0.5 * cell_size_row;
       temp_coord(2) = (j-1) * cell_size_col + 0.5 * cell_size_col;
       grid_center_coord{i,j} = temp_coord;
   end
end

kBasisMat =[ 1   26  66  26    1 0;
            -5  -50  0   50  5 0;
            10  20 -60  20  10 0;
           -10  20  0  -20  10 0;
             5 -20  30 -20  5 0;
            -1  5  -10  10 -5 1;];
kBasisMat = kBasisMat/factorial(5);         
         
%% check all the patterns
start_index = [7 7]; %locate initial point at center so that all the patterns are on the grid
direction_map = [ 
%  0  0;
  1  0;
  1 -1;
  0 -1;
 -1 -1;
 -1  0;
 -1  1;
  0  1;
  1  1;  
];

num_connect = size(direction_map,1);
num_gen = POLYORDER;

combs = permn(1:num_connect, num_gen);
num_total_comb = size(combs,1);
pattern_coord = cell(num_total_comb,1);
pattern_index = cell(num_total_comb,1);
for i = 1 : num_total_comb
   %for each comb/pattern
   current_pattern = zeros(num_gen + 1, 2);
   current_pattern(1,:) = start_index;
   current_comb = combs(i,:);
   for step = 1 : num_gen
        current_pattern(step+1,:) = current_pattern(step,: ) + direction_map(current_comb(step),:);
   end
   %transform current_pattern to coord
   current_coord_pattern = grid2coord(current_pattern, grid_center_coord);
   %store all the patterns
   pattern_coord{i}  = current_coord_pattern;
   pattern_index{i}  = current_pattern;
end

%for B-spline evaluation
ut_step = 0.02;
utVec = 0: ut_step :1.0;

delta_vec = 0.0: 0.01: 0.5 * cell_size_row;
for  delta = delta_vec
    is_exceed = 0; %for each delta, by default exceed flag is false;
    for i = 1: num_total_comb
        pattern = pattern_coord{i};
        current_index = pattern_index{i};
        for ut = utVec % for each sampling point of the trajectory of the pattern
            is_incell_found = 0; %by default, no in grid is found
            p_ut = [1 ut ut^2 ut^3 ut^4 ut^5] * kBasisMat * pattern; %evaluate position   
            for k = 1:6 %for all known collision free cells.
                temp_coord = grid_center_coord{current_index(k,1),current_index(k,2)}; %get coordinate
                if( abs(p_ut(1) -temp_coord(1) ) <= (0.5 * cell_size_row + delta) &&   abs(p_ut(2) -temp_coord(2) ) <= (0.5 * cell_size_col+delta)   ) %find a cell that contains current point with margin delta
                    is_incell_found = 1; 
                end
            end

            if(is_incell_found == 0) %if no cell find, then the trajectory deviates from collision-free space
               % disp('deviates from grid..');
                is_exceed = 1;
                break;
            end
        end
    end

    if(is_exceed) 
        delta
        disp('Collision-free condition is violated');
    else
        delta
        disp('All patterns are collision free!');
    end
end
