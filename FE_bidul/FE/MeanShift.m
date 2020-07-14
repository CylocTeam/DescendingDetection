function [pts_shifted,unconverged_pts] = MeanShift(pts, cluster_bw ,cluster_threshold, neighbour_dist, convergence_threshold)
cluster_bw = 2;
neighbour_dist = 1.7;
convergence_threshold = 0.1;
MAX_ITERS = 100;
% pts = pts(:,1:2);


if isvector(pts) && ~iscolumn(pts)
    pts = pts';
end
kernel_fcn = @(r) exp(-r.^2);

pts_shifted = pts;
NumOfPoints = size(pts, 1);
NumOfDims = size(pts, 2);

num_iters = 1;
unconverged_idx = 1:NumOfPoints;
while  ~isempty(unconverged_idx) && num_iters < MAX_ITERS
    
    unconverged_pts = pts_shifted(unconverged_idx,:);
    curr_dists = pdist2(pts, unconverged_pts);
    valid_neighbour_points_mat = (curr_dists <= neighbour_dist);
    curr_pts_weights_mat = kernel_fcn(curr_dists/cluster_bw);
    curr_pts_weights_mat(~valid_neighbour_points_mat) = 0;
    curr_pts_weights_mat = curr_pts_weights_mat ./ sum(curr_pts_weights_mat,1);
    newPointss = curr_pts_weights_mat' * pts ;
    shiftsteps = abs(newPointss - unconverged_pts);
    pts_shifted(unconverged_idx,:) = newPointss;
    
    has_converged = vecnorm(shiftsteps,2,2)/NumOfDims < convergence_threshold;
    unconverged_idx = unconverged_idx(~has_converged);
    num_iters = num_iters + 1;
    
%     curr_pt1 = pts_shifted(pts_idx, :);
%     curr_dist = pdist2(pts, curr_pt1);
%     valid_neighbour_points = (curr_dist <= neighbour_dist);
%     curr_pts_weights = kernel_fcn(curr_dist(valid_neighbour_points)/cluster_bw);
%     curr_pts_weights = curr_pts_weights / sum(curr_pts_weights);
%     newPoint = sum( repmat(curr_pts_weights, 1, NumOfDims).* pts(valid_neighbour_points,:), 1 );
%     shiftstep = abs(newPoint - pts_shifted(pts_idx, :));
%     pts_shifted(pts_idx, :) = newPoint;
end

unconverged_pts = pts_shifted(unconverged_idx,:);
    
% figure;
% scatter(pts(:,1),pts(:,2),80,'filled'); hold all
% scatter(pts_shifted(:,1),pts_shifted(:,2),50,'filled')
% grid minor
%% reference alg
% pts_shifted2 = pts;
% for pts_idx = 1:NumOfPoints
%     shiftstep = inf;
%     while  norm(shiftstep)/NumOfDims > convergence_threshold % for shift_iter = 1:mean_shift_iters %
%         curr_pt1 = pts_shifted2(pts_idx, :);
%         curr_dist = pdist2(pts, curr_pt1);
%         valid_neighbour_points = (curr_dist <= neighbour_dist);
%         curr_pts_weights = kernel_fcn(curr_dist(valid_neighbour_points)/cluster_bw);
%         curr_pts_weights = curr_pts_weights / sum(curr_pts_weights);
%         newPoint = sum( repmat(curr_pts_weights, 1, NumOfDims).* pts(valid_neighbour_points,:), 1 );
%         shiftstep = abs(newPoint - pts_shifted2(pts_idx, :));
%         pts_shifted2(pts_idx, :) = newPoint;
%     end
% end
% disp(norm(pts_shifted2 - pts_shifted))

