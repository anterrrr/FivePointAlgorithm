function [normalCoordinate] = Normalization(intrinsicParameters,coordinate)
%% 函数解释：已知内参数和畸变矫正过的像面坐标求解归一化的像面坐标
k = [-intrinsicParameters(3),0,0;
    0,-intrinsicParameters(3),0;
    0,0,1];

normalCoordinate = [];
for i=1:length(coordinate)
    temp = inv(k)*[coordinate(i,1);coordinate(i,2);1];
    normalCoordinate = [normalCoordinate;temp'];
end
end