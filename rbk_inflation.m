clc; clear all; close all;
% This program is to check the inflation needed to guarantee that the RBK
% trajectory is collision-free. 
% Author: Wenchao DING, HKUST, wdingae@connect.ust.hk
%% paremeters
cell_size_row = 0.16;
cell_size_col = 0.16;
grid_size_row = 14;
grid_size_col = 14;
POLYORDER = 5;

% input the inflation you what to check
%delta_vec = 0.0: 0.01: 0.03; %check a set of inflation.
delta_vec = 0.2 * cell_size_row; % check a given inflation.

raw_image = imread('grid_14x14.png'); %load background image 
resize_coeff = 0.5;
resized_image = imresize(raw_image,resize_coeff);
gray_image = rgb2gray(resized_image);
[pixel_size_row, pixel_size_col] = size(gray_image);
%for visualization, we use image background
pixel2grid_scale_col = pixel_size_col/grid_size_col;
pixel2grid_scale_row = pixel_size_row/grid_size_row; 

%record the coord of centers.
pixel_center_coord = cell(grid_size_row, grid_size_col);
grid_center_coord  = cell(grid_size_row, grid_size_col);
temp_coord = zeros(2,1);
for i = 1:1:grid_size_row
   for j = 1:1:grid_size_col
       %compute pixel index
        temp_coord(1) =  (j - 1) * cell_size_col + 0.5 * cell_size_col;
        temp_coord(2) =  (i - 1) * cell_size_row + 0.5 * cell_size_row;
        grid_center_coord{i,j} = temp_coord; % 2d world coordinate
        
        temp_coord(1) =  (j - 1) * pixel2grid_scale_col + 0.5 * pixel2grid_scale_col;
        temp_coord(2) =  (i - 1) * pixel2grid_scale_row + 0.5 * pixel2grid_scale_row;
        pixel_center_coord{i,j} = temp_coord; %image coordinate for visualization and plot       
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
idx = 1;
for  delta = delta_vec
    idx
    h = figure(idx);
    set(0,'CurrentFigure',h)
    rbk_layer = imshow(resized_image); hold on;
    idx = idx + 1;
    for i = 1:1:grid_size_row
       for j = 1:1:grid_size_col
           %compute pixel index
           temp_coord = pixel_center_coord{i,j};        

           left_top   = [temp_coord(1) - (pixel2grid_scale_col - delta * pixel2grid_scale_col/cell_size_col)*0.5,...
                         temp_coord(2) - (pixel2grid_scale_row - delta * pixel2grid_scale_col/cell_size_col)*0.5  ];
           right_top  = [temp_coord(1) + (pixel2grid_scale_col - delta * pixel2grid_scale_col/cell_size_col)*0.5,...
                         temp_coord(2) - (pixel2grid_scale_row - delta * pixel2grid_scale_col/cell_size_col)*0.5  ];
           right_down = [temp_coord(1) + (pixel2grid_scale_col - delta * pixel2grid_scale_col/cell_size_col)*0.5,...
                         temp_coord(2) + (pixel2grid_scale_row - delta * pixel2grid_scale_col/cell_size_col)*0.5  ];
           left_down  = [temp_coord(1) - (pixel2grid_scale_col - delta * pixel2grid_scale_col/cell_size_col)*0.5,...
                         temp_coord(2) + (pixel2grid_scale_row - delta * pixel2grid_scale_col/cell_size_col)*0.5  ];
           X = [left_top(1), right_top(1), right_down(1), left_down(1)];
           Y = [left_top(2), right_top(2), right_down(2), left_down(2)];       
           fill(X,Y,'b','LineStyle','none');
       end
    end    
    
    max_pattern_deviation = -inf; %compute the maximum deviation for all possible patterns
    max_pattern_index     = 0;
    
    is_exceed = 0; %for each delta, by default exceed flag is false;
    for i = 1: num_total_comb
        pattern_deviation = -inf;  %for each pattern, compute the deviation
        
        pattern = pattern_coord{i};
        current_index = pattern_index{i};
        for ut = utVec % for each sampling point of the trajectory of the pattern
            sample_deviation = inf; %for each sampling point, compute the deviation
            is_incell_found = 0; %by default, no in cell is found
            p_ut = [1 ut ut^2 ut^3 ut^4 ut^5] * kBasisMat * pattern; %evaluate position   
            
            for k = 1:6 %for all known collision free cells.
                temp_coord = pattern(k,:); %get coordinate
                if( abs(p_ut(1) -temp_coord(1) ) <= (0.5 * cell_size_row + delta) &&   abs(p_ut(2) -temp_coord(2) ) <= (0.5 * cell_size_col+delta)   ) %find a cell that contains current point with margin delta
                    is_incell_found = 1; 
                end
                cur_cell_deviaton =  max( abs(p_ut(1) -temp_coord(1) ),  abs(p_ut(2) -temp_coord(2) )  ); %deviation metric (l-infinite norn)
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
    disp(['Current inflation level:' ...
           num2str(delta) ]);
    
    disp(['Minimum inflation required:' ...
           num2str(max_pattern_deviation - 0.5 * min(cell_size_col, cell_size_row) ) ]);
    
    disp('The pattern that has the maximum deviation:');
    max_pattern_show = pattern_index{max_pattern_index}
    
    max_pattern_grid_coord  =  zeros(size(max_pattern_show) );
    max_pattern_pixel_coord =  zeros(size(max_pattern_show) );
    
    for i = 1 : size(max_pattern_grid_coord,1)
       max_pattern_grid_coord(i,:)  = grid_center_coord{max_pattern_show(i,1), max_pattern_show(i,2)};
       max_pattern_pixel_coord(i,:) = pixel_center_coord{max_pattern_show(i,1), max_pattern_show(i,2)};
       
       left_top   = [max_pattern_pixel_coord(i,1) - (pixel2grid_scale_col + delta * pixel2grid_scale_col/cell_size_col)*0.5,...
                     max_pattern_pixel_coord(i,2) - (pixel2grid_scale_row + delta * pixel2grid_scale_row/cell_size_row)*0.5  ];
       right_top  = [max_pattern_pixel_coord(i,1) + (pixel2grid_scale_col + delta * pixel2grid_scale_col/cell_size_col)*0.5,...
                     max_pattern_pixel_coord(i,2) - (pixel2grid_scale_row + delta * pixel2grid_scale_row/cell_size_row)*0.5  ];
       right_down = [max_pattern_pixel_coord(i,1) + (pixel2grid_scale_col + delta * pixel2grid_scale_col/cell_size_col)*0.5,...
                     max_pattern_pixel_coord(i,2) + (pixel2grid_scale_row + delta * pixel2grid_scale_row/cell_size_row)*0.5  ];
       left_down  = [max_pattern_pixel_coord(i,1) - (pixel2grid_scale_col + delta * pixel2grid_scale_col/cell_size_col)*0.5,...
                     max_pattern_pixel_coord(i,2) + (pixel2grid_scale_row + delta * pixel2grid_scale_row/cell_size_row)*0.5  ];                 
       X = [left_top(1), right_top(1), right_down(1), left_down(1)];
       Y = [left_top(2), right_top(2), right_down(2), left_down(2)];       
       fill(X,Y,'g','LineStyle','none');
    end
    
    plot(max_pattern_pixel_coord(:,1), max_pattern_pixel_coord(:,2),'o','MarkerSize',5,'MarkerFaceColor',[1.0 0.25 1.0]);
    plot(max_pattern_pixel_coord(:,1), max_pattern_pixel_coord(:,2),'k');
    
    all_points = [];
    for ut = utVec % for each sampling point of the trajectory of the pattern
        p_ut = [1 ut ut^2 ut^3 ut^4 ut^5] * kBasisMat * max_pattern_grid_coord; %evaluate position   
        temp_pixel_coord = [pixel2grid_scale_col/cell_size_col * p_ut(1), pixel2grid_scale_row/cell_size_row * p_ut(2)];
        all_points = [all_points; temp_pixel_coord];
    end    
    plot(all_points(:,1), all_points(:,2),'r','LineWidth',1.5);
    title(['current inflation: ' num2str(delta)]);
    hold off;
    if(is_exceed) 
        disp('Collision-free condition is violated');
    else
        disp('All patterns are collision free!');
    end
    
    pause(0.2)
end
