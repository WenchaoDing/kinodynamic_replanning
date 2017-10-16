clc; clear all; close all;
%% parameters
cell_size_row    = 0.16;
cell_size_col    = 0.16;
cell_size_height = 0.16;
grid_size_row    = 14;
grid_size_col    = 14;
grid_size_height = 14;
%visualize occupancy grid
POLYORDER = 5;

grid_center_coord  = cell(grid_size_row, grid_size_col, grid_size_height);
temp_coord = zeros(3,1);
for i = 1:1:grid_size_row
   for j = 1:1:grid_size_col
       for k = 1:1:grid_size_height
            %compute pixel index
            temp_coord(1) = (j - 1) * cell_size_col    + 0.5 * cell_size_col;
            temp_coord(2) = (i - 1) * cell_size_row    + 0.5 * cell_size_row;
            temp_coord(3) = (k - 1) * cell_size_height + 0.5 * cell_size_height;
            grid_center_coord{i,j,k} = temp_coord;
       end
   end
end

kBasisMat =[ 1   26  66   26    1  0;
            -5  -50   0   50    5  0;
            10   20 -60   20   10  0;
           -10   20   0  -20   10  0;
             5  -20  30  -20    5  0;
            -1    5 -10   10   -5  1;];
kBasisMat = kBasisMat/factorial(5);         
         
%% check all the patterns

start_index = [7 7 7]; %locate initial point at center so that all the patterns are on the grid
direction_map = [ 
  1  1   -1;  
  1  1    0;
  1  1    1;      
  1  0  -1;
  1  0   0;  
  1  0   1;  
  1 -1  -1;
  1 -1   0;
  1 -1   1;  
  0  0  -1;
  0  0   1;    
  0 -1  -1;
  0 -1   0;
  0 -1   1;  
  0  1   -1;
  0  1    0;
  0  1   1;    
 -1 -1  -1;
 -1 -1   0;
 -1 -1   1; 
 -1  0  -1;
 -1  0   0;
 -1  0   1; 
 -1  1  -1;
 -1  1   0;
 -1  1   1; 
  ];

num_connect = size(direction_map,1);
num_gen = POLYORDER;

combs = permn(1:num_connect, num_gen);
num_total_comb = size(combs,1);
pattern_coord = cell(num_total_comb,1);
pattern_index = cell(num_total_comb,1);

%for B-spline evaluation
ut_step = 0.02;
utVec = 0: ut_step :1.0;
%delta_vec = 0.0: 0.01: 0.5 * cell_size_row;
delta_vec = 0.2 * cell_size_row;

run = 0; %turn on this switch to run this part, may need a long time to output results.
%or you can find the results generated below.
if run
    disp('Start to generate pattern info..');
    for i = 1 : num_total_comb
       %for each comb/pattern
       current_pattern = zeros(num_gen + 1, 3);
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

    disp('Pattern info recorded..');


    for  delta = delta_vec

        max_pattern_deviation = -inf; %compute the maximum deviation for all possible patterns
        max_pattern_index     = 0;

        is_exceed = 0; %for each delta, by default exceed flag is false;
        for i = 1: num_total_comb
            disp(['current comb: ' num2str(i) ]);
            disp(['Total comb: '   num2str(num_total_comb) ]);
            pattern_deviation = -inf;  %for each pattern, compute the deviation
            pattern = pattern_coord{i};
            current_index = pattern_index{i};
            for ut = utVec % for each sampling point of the trajectory of the pattern
                sample_deviation = inf; %for each sampling point, compute the deviation
                is_incell_found = 0; %by default, no in cell is found
                p_ut = [1 ut ut^2 ut^3 ut^4 ut^5] * kBasisMat * pattern; %evaluate position   

                for k = 1:6 %for all known collision free cells.
                    temp_coord = pattern(k,:); %get coordinate
                    if( abs(p_ut(1) -temp_coord(1) ) <= (0.5 * cell_size_row    + delta) && ...
                        abs(p_ut(2) -temp_coord(2) ) <= (0.5 * cell_size_col    + delta) && ...
                        abs(p_ut(3) -temp_coord(3) ) <= (0.5 * cell_size_height + delta)) %find a cell that contains current point with margin delta
                        is_incell_found = 1; 
                    end
                    cur_cell_deviaton =  max( max(abs(p_ut(1) -temp_coord(1)), abs(p_ut(2) -temp_coord(2))), abs(p_ut(3) -temp_coord(3)) ); %deviation metric
                    if(cur_cell_deviaton < sample_deviation)
                        sample_deviation       =  cur_cell_deviaton;
                    end

                end

                if(sample_deviation > pattern_deviation)
                    pattern_deviation = sample_deviation;
                end

                if(is_incell_found == 0) %if no cell find, then the sampling point deviates from collision-free space
                   % disp('deviates from grid..');
                    is_exceed = 1;
                end
            end
            if(pattern_deviation > max_pattern_deviation)
                max_pattern_deviation = pattern_deviation;
                max_pattern_index     = i;
            end

        end

        delta
        max_pattern_deviation
        max_pattern_index
        max_pattern_show = pattern_index{max_pattern_index}

        if(is_exceed) 
            disp('Collision-free condition is violated');
        else
            disp('All patterns are collision free!');
        end
    end

end


% delta =
% 
%     0.0320

% max_pattern_deviation =
% 
%     0.0999

% max_pattern_index =
% 
%      2847367

% max_pattern_show =
% 
%      7     7     7
%      8     6     6
%      9     5     5
%     10     6     4
%     11     7     5
%     12     8     6

%% visualization

max_pattern_deviation = 0.0999;

max_pattern_show = [  
     7     7     7
     8     6     6
     9     5     5
    10     6     4
    11     7     5
    12     8     6];

disp(['Current inflation level:' ...
       num2str(delta_vec(1)) ]);

disp(['Minimum inflation required:' ...
       num2str(max_pattern_deviation - 0.5 * min(min(cell_size_col, cell_size_row), cell_size_height)) ]);

max_pattern_coord = grid2coord(max_pattern_show, grid_center_coord);
figure;


plot3(max_pattern_coord(:,1), max_pattern_coord(:,2),max_pattern_coord(:,3),'.','MarkerSize',10 ); hold on;
plot3(max_pattern_coord(:,1), max_pattern_coord(:,2),max_pattern_coord(:,3),'-gs','LineWidth',2);

all_points = [];
for ut = utVec % for each sampling point of the trajectory of the pattern
    p_ut = [1 ut ut^2 ut^3 ut^4 ut^5] * kBasisMat * max_pattern_coord; %evaluate position   
    all_points = [all_points; p_ut];
end
plot3(all_points(:,1), all_points(:,2),all_points(:,3),'r','LineWidth',1);

set(gca,'Xtick',0:cell_size_col:cell_size_col*grid_size_col);
set(gca,'Ytick',0:cell_size_row:cell_size_row*grid_size_row);
set(gca,'Ztick',0:cell_size_height:cell_size_col*grid_size_height); 
set(gca,'XTickLabel','');
set(gca,'YTickLabel','');
set(gca,'ZTickLabel','');
xlim([0 grid_size_col*cell_size_col]);
ylim([0 grid_size_row*cell_size_row]);
zlim([0 grid_size_height*cell_size_height]);
grid on; box on;