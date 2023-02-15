function Write_Data(fileName,data,formatSpec)

fileID = fopen(fileName,'w');
for ii = 1:length(data)
    data_now = data(ii,:);
    fprintf(fileID,formatSpec,data_now);
end
fclose(fileID);

%     save(sprintf('Project1_SC%d_Measurements_Raw.txt',k),'measurement_data_raw','-ascii')
%     save(sprintf('Project1_SC%d_Measurements_Noisy.txt',k),'measurement_data_noisy','-ascii')
%     save(sprintf('Project1_SC%d_Trajectory_Data.txt',k),'trajectory_data','-ascii')