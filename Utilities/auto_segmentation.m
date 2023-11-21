

function Simu_Labels = auto_segmentation(          Fetal_Brain_zp, ...
                                                      Fetal_Brain, ...
                                                           SimRes, ...
                                                  sampling_factor, ...
                                                   SliceThickness, ...
                                                       SubunitRes, ...
                                                          FOVRead, ...
                                                         FOVPhase, ...
                                                         NbSlices, ...
                                          interleavedSlices_index, ...
                                                         Sl_to_Sl, ...
                                          motion_corrupted_slices, ...
                                         translation_displacement, ...
                                                   rotation_angle, ...
                                                    rotation_axis)

% Input check
if nargin < 15
    error('Missing input(s).');
elseif nargin > 15
    error('Too many inputs.');
end

% Memory pre-allocation
Sl_Volume = zeros(size(Fetal_Brain_zp,1), size(Fetal_Brain_zp,2), NbSlices);

% Initialization of tranformation arrays
motion_index = 0;

% Loop through slices to apply motion tranform
for iSlice=1:length(interleavedSlices_index)
    disp(['Slice ', num2str(interleavedSlices_index(iSlice)), ' of ', num2str(length(interleavedSlices_index))])
    index = interleavedSlices_index(iSlice)*length(Sl_to_Sl)-(length(Sl_to_Sl)-1);
    % Random inter-slice motion: be careful, labels are interpolated!
    if any(interleavedSlices_index(iSlice)==motion_corrupted_slices)
        motion_index = motion_index + 1;
        Fetal_Brain_rotated = moving_3Dbrain(                             Fetal_Brain, ...
                                                                               SimRes, ...
                                             translation_displacement(motion_index,:), ...
                                                                            'nearest', ...
                                                         rotation_angle(motion_index), ...
                                                        rotation_axis(motion_index,:), ...
                                                                            'nearest', ...
                                                                               'crop');
        Fetal_Brain_rotated_upsampled = sampling_OoP(Fetal_Brain_rotated, ...
                                                         sampling_factor, ...
                                                               'nearest');
        Fetal_Brain_moved_zp = Resize_Volume(Fetal_Brain_rotated_upsampled, size(Fetal_Brain_zp));
        Fetal_Brain_zp(:,:,index:index+round(SliceThickness/SubunitRes)-1) = Fetal_Brain_moved_zp(:,:,index:index+round(SliceThickness/SubunitRes)-1);
        % Update the T2decay_zp variable as it is unlikely that the
        % fetus comes back to its initial position at iSlice+1 and
        % following slices after it has moved
        if iSlice < length(interleavedSlices_index)
            for nextSlice=iSlice+1:length(interleavedSlices_index)
                nextSlice_index = interleavedSlices_index(nextSlice)*length(Sl_to_Sl)-(length(Sl_to_Sl)-1);
                Fetal_Brain_zp(:,:,nextSlice_index:nextSlice_index+round(SliceThickness/SubunitRes)-1) = Fetal_Brain_moved_zp(:,:,nextSlice_index:nextSlice_index+round(SliceThickness/SubunitRes)-1);
            end
        end    
        clear Fetal_Brain_rotated
        clear Fetal_Brain_rotated_upsampled
    end
end

%

freq_labels = 0;
freq_labels_other = 0;

% Now that motion occured and that the labels are right, we keep only one
% label value within the slice thickness
for iSlice=1:length(interleavedSlices_index)
	disp(['Slice ', num2str(interleavedSlices_index(iSlice)), ' of ', num2str(length(interleavedSlices_index))])
    index = interleavedSlices_index(iSlice)*length(Sl_to_Sl)-(length(Sl_to_Sl)-1);
    % Interpolate the label for all voxels at the same (x,y) location
    % across the slice thickness direction
    for row=1:size(Fetal_Brain_zp,1)
        for col=1:size(Fetal_Brain_zp,2)
            slice_labels = Fetal_Brain_zp(row,col,index:index+round(SliceThickness/SubunitRes)-1);
            if unique(slice_labels)==0
                max_value = 0;
            else
                % Count the occurrences of each non-zero label within the
                % slice thickness
                [count, value] = groupcounts(squeeze(slice_labels(slice_labels>0)));
                max_index = find(count==max(count));
                if length(max_index)==1
                    max_value = value(max_index);
                %Handle cases where more than 1 label occurs at the same
                %frequency in the slice thickness
                else
                    freq_labels = freq_labels + 1;
                    if any(value(max_index)==Fetal_Brain_zp(row,col,index-1))
                        freq_labels_other = freq_labels_other + 1;
                        max_value = Fetal_Brain_zp(row,col,index-1);
                    elseif any(value(max_index)==Fetal_Brain_zp(row,col,index+round(SliceThickness/SubunitRes)))
                        freq_labels_other = freq_labels_other + 1;
                        max_value = Fetal_Brain_zp(row,col,index+round(SliceThickness/SubunitRes));
    %                 elseif any(value(max_index)==Fetal_Brain_zp(row-1,col,round(index+SliceThickness/SubunitRes-1)/2))
    %                     freq_labels_other = freq_labels_other + 1;
    %                     max_value = Fetal_Brain_zp(row-1,col,round(index+SliceThickness/SubunitRes-1)/2);
    %                 elseif any(value(max_index)==Fetal_Brain_zp(row+1,col,round(index+SliceThickness/SubunitRes-1)/2))
    %                     freq_labels_other = freq_labels_other + 1;
    %                     max_value = Fetal_Brain_zp(row+1,col,round(index+SliceThickness/SubunitRes-1)/2);
    %                 elseif any(value(max_index)==Fetal_Brain_zp(row,col-1,round(index+SliceThickness/SubunitRes-1)/2))
    %                     freq_labels_other = freq_labels_other + 1;
    %                     max_value = Fetal_Brain_zp(row,col-1,round(index+SliceThickness/SubunitRes-1)/2);
    %                 elseif any(value(max_index)==Fetal_Brain_zp(row,col+1,round(index+SliceThickness/SubunitRes-1)/2))
    %                     freq_labels_other = freq_labels_other + 1;
    %                     max_value = Fetal_Brain_zp(row,col+1,round(index+SliceThickness/SubunitRes-1)/2);
                    end
                end
            end
            Fetal_Brain_zp(row,col,index:index+round(SliceThickness/SubunitRes)-1) = max_value;
        end
    end
%     [count, value] = groupcounts(slice_labels);
%     max_index = find(count==max(count));
%     max_value = value(max_index);
%     Fetal_Brain_zp(100,100,index:index+SliceThickness/SubunitRes-1) = max_value;
    Sl_Volume(:,:,interleavedSlices_index(iSlice)) = Fetal_Brain_zp(:,:,index);
end

% Recover the motion-corrupted label map at the right simulated resolution
% of acquisition
Simu_Labels = Resize_Volume(Sl_Volume, [round(FOVRead)/SimRes(1), round(FOVPhase)/SimRes(2), NbSlices]);

end