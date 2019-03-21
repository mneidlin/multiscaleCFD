function A = readpv()
%Reads two column txt pv.txt and plots the shit!
clear all;
%fileID = fopen('pv.txt','r');
fileID = fopen('results.txt','r');
%formatSpec = '%f %f';
formatSpec = '%f %f %f';
%sizeA = [2 Inf];
sizeA = [3 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A=A';
%plot(A(:,1),A(:,2));
plot(A(:,1),A(:,2),A(:,1),A(:,3));
title('LV p-v Loop');
end

