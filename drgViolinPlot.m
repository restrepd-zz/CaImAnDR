function drgViolinPlot(x_pos,all_points,edges,rand_offset,which_color_mean,which_color_point,point_size)
%This does a violin plot for several points

for ii=1:size(all_points,1)
    [mean_out, CIout]=drgViolinPoint(all_points(ii,:),edges,x_pos(ii),rand_offset,which_color_mean,which_color_point,point_size);
end

