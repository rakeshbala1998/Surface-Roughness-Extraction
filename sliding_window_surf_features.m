clc;
dirpath="D:\Research\CRACK_LOC_FEATURES\g6a\files";
% Get a list of all the .txt files in the directory
file_list = dir(fullfile(dirpath, '*.csv'));
% Create an empty table with column names
results_table = table('Size',[0,8],'VariableTypes',{'string','double','double','double','double','double','double','double'},...
'VariableNames',{'Sample','angle','min_y','max_y','Ra','Ry','Rz','Rho'});
% Loop through the list and read each file
for m = 1:length(file_list)
    % Get the filename
    filename = file_list(m).name;
    fprintf(filename)
    
    % Read the file
    surface = readmatrix(fullfile(dirpath, filename));
    numStr = regexp(filename, '\d+', 'match');
    angle = str2double(numStr{1,2});
        
    sample='amga6';
    %find the maximum and minimum point 
    min_y=min(surface(:,2));
    max_y=max(surface(:,2));
    
    for j=min_y:5:max_y
    
        %Taking only the right side of the right
        x_indices = find(surface(:,1)>=3);
        %Extract the values
        x_filtered_matrix=surface(x_indices,:);
        y_indices=find(x_filtered_matrix(:,2)>j & x_filtered_matrix(:,2)<j+5);
        y_filtered_matrix=x_filtered_matrix(y_indices,:);
        %Defining x and y
        x=y_filtered_matrix(:,2);
        y=y_filtered_matrix(:,1)-3.5; %shift the roughness to 0 
        %check if it is empty 
        if isempty(x)
            continue;
        end
        if isempty(y)
            continue;
        end
        %Find the average value
        actual_avg=mean(y(:,1));
        % Find the indices of the peaks and valleys
        [pks, locs] = findpeaks(y);
        % Find the indices of the valleys
        [vls, vlocs] = findpeaks(-y);
        vls = -vls;
        % Sort the peaks and valleys by height
        [~, pk_sort_idx] = sort(pks, 'descend');
        [~, vl_sort_idx] = sort(vls, 'ascend');
        % Get the top 5 peaks and valleys
        top_5_pks = pks(pk_sort_idx(1:5));
        top_5_pk_locs = locs(pk_sort_idx(1:5));
        top_5_vls = vls(vl_sort_idx(1:5));
        top_5_vl_locs = vlocs(vl_sort_idx(1:5));
        all_vls = vls(vl_sort_idx(1:length(vl_sort_idx)));
        all_vl_locs = vlocs(vl_sort_idx(1:length(vl_sort_idx)));
        % Plot the surface roughness data and highlight the top 5 peaks and valleys
%         figure;
%         plot(x, y);
%         hold on
%         scatter(x(top_5_pk_locs), top_5_pks, 'ro', 'filled');
%         scatter(x(top_5_vl_locs), top_5_vls, 'go', 'filled');
        xlabel('X');
        ylabel('Y');
        legend('Surface Roughness', 'Top 5 Peaks', 'Top 5 Valleys');
        
        %Average Roughness
        max_x=max(y_filtered_matrix(:,1));
        min_x=min(y_filtered_matrix(:,1));
        Ry= max_x-min_x;
        fprintf("The Average Roughness Ra: %f mm\n",actual_avg)
        fprintf("Peak to Valley Height Ry: %f mm\n",Ry)
        b=sort(unique(y),'descend');
        max5=b(1:5);
        min5=b(end-4:end);
        Rz=1/5*(sum(max5(:))+sum(min5(:)));
        fprintf("10-point Roughness Rz: %f mm\n",Rz)
        root_radius=zeros(1,5);
              
        for k=1:5
            actual_index=top_5_vl_locs(k);
            
        
            
            % Search forward from the valley
            forward_indices = [];
            ignore_positives = true;
            ignore_zeros = true;
            last_slope = 0;
            for i = actual_index:length(y)
                slope = (y(i) - y(i-1))/(x(i)-x(i-1));
                %fprintf("slope:%d",slope)
                if ignore_positives && slope > 0
                    continue;  % Ignore initial positive slopes
                else
                    ignore_positives = false;
                end
                if ignore_zeros && slope == 0
                    continue;  % Ignore initial zeros
                else
                    ignore_zeros = false;
                end
                if slope ==Inf
                    continue;
                end
                if slope >= 0 && slope ~= +Inf || slope == -Inf  
                    break;  % Stop if slope is positive 
                else
                    forward_indices(end+1) = i;
                end
                
            end
        
            %plot(x(forward_indices(end)),y(forward_indices(end)),'go')
        
            % Search backward from the valley
            backward_indices = [];
            ignore_negatives = true;
            ignore_zeros = true;
            for i = actual_index:-1:1
                %fprintf("i = %d &  i+1 = %d",i, i+1)
                slope = (y(i+1) - y(i))/(x(i+1) - x(i));
                %fprintf("slope backward :%f",slope)
                if ignore_negatives && slope < 0
                    continue;  % Ignore initial negative slopes
                else
                    ignore_negatives = false;
                end
                if ignore_zeros && slope == 0
                    continue;  % Ignore initial zeros
                else
                    ignore_zeros = false;
                end
                if slope ~= +Inf && slope <= 0 || slope == -Inf % Stop if slope is Negative
                    break;
                else
                    backward_indices(end+1)=i;
               end
            end
        
            if isempty(backward_indices)
                continue;
            end
        
            if isempty(forward_indices)
                continue;
            end
        
            end_indices=forward_indices(end);
            start_indices=backward_indices(end);
            
            y_end_indices=y(end_indices);
            y_start_indices= y(start_indices);
            
            if y_start_indices < y_end_indices
                threshold = 0.001; % Set a threshold value
                new_end_indices = end_indices; % Use old end indices as default
                for u = 1:length(forward_indices)
                    if abs(y(forward_indices(u)) - y_start_indices) < threshold
                        new_end_indices = forward_indices(u);
                        break;
                    end
                end
            else
                new_end_indices = end_indices;
            end
            
            if y_start_indices > y_end_indices
                threshold = 0.001; % Set a threshold value
                new_start_indices = start_indices; % Use old start indices as default
                for u = 1:length(backward_indices)
                    if abs(y(backward_indices(u)) - y_end_indices) < threshold
                        new_start_indices = backward_indices(u);
                        break;
                    end
                end
            else
                new_start_indices = start_indices;
            end
        
        
            % Now we have the adaptive valley with the minimum point index
            x_valley = x(start_indices:end_indices);
            y_valley = y(start_indices:end_indices);
        
            x_valley_new = x(new_start_indices:new_end_indices);
            y_valley_new = y(new_start_indices:new_end_indices);
            
            
            
            
            % Fit a polynomial to the data
            p = polyfit(x_valley_new, y_valley_new, 2);
            
            % Find the first derivative of the polynomial using polyder
            p_prime = polyder(p);
            
            % Find the second derivative of the polynomial using polyder
            p_double_prime = polyder(p_prime);
            
            % Find the index of the lowest point in the valley
            [y_min, min_index] = min(y_valley_new);
            xmin = fminbnd(@(x_valley_new) polyval(p, x_valley_new), min(x_valley_new), max(x_valley_new));
        
        
            % Evaluate the second derivative of the polynomial for a given x1 value
            %d2y_dx2 = polyval(p_double_prime, x_valley(min_index));
            %dy_dx=polyval(p_prime, x_valley(min_index));
            curvature= abs(polyval(p_double_prime,xmin)) / (1 + polyval(p_prime, xmin)^2)^(3/2);
            
            %root radius
            root_radius(1,k)=(1/curvature);
            fprintf('Rho for valley %d: %f mm \n',k,root_radius(1,k));
            
            %plot the fiting 
            % Define the x-axis
            x_axis = linspace(min(x_valley_new), max(x_valley_new), 100);
            
            % Evaluate the polynomial at each x-value
            y_fit = polyval(p, x_axis);
            
            %centre of the circle
            radius=root_radius(1,k);
            center=[x_valley_new(min_index),y_valley_new(min_index)+radius];

            if radius>0.2
                root_radius(1,k)=0;
            end
            % Plot the circle and the fitted polynomial curve
%             figure;
%             plot(x(actual_index-5:actual_index+5),y(actual_index-5:actual_index+5),'r')
%             hold on
%             plot(x(actual_index),y(actual_index),'bo')
%             plot(x_axis, y_fit, 'g');
%             plot(x_valley, y_valley, 'b--')
%             plot(x_valley_new,y_valley_new,'black')
%             viscircles(center, radius);
%             legend('manual selection','lowest point', 'adaptive selection','polyfit','location','northwest')
%             xlabel('X');
%             ylabel('Y');
            
        end

            
           % Append the results for this interval to the table
        
        fprintf("average rho of dominant valleys : %f mm\n",mean(root_radius));
        results_table= [results_table; {sample, angle,j,j+5, actual_avg, Ry, Rz, root_radius}];
     end
        
        
    
    
end
%%
% Save the table with index
path= "C:/Users/rakes/Dropbox (ASU)/Rakesh AM FEA/Data extraction";
temp=[path+filesep+sample+'_surf_new_r'+'.csv'];  
writetable(results_table, temp, 'WriteVariableNames', true);