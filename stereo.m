function [ ] = stereo( Left, Right, Scale, Map )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

occlusion = Scale;
scale = 1;
%L=Left;
%Left=Right;
%Right=L;
% ===== Initializing ====== %

% Read the input Image Files
RightImage = imread(Right);
LeftImage = imread(Left);

RI = double(RightImage);LI = double(LeftImage);
for i=1:size(RI,1),
    for j=1:size(RI,2),
        RI(i,j,:) = RI(i,j,:)/norm([RI(i,j,1),RI(i,j,2),RI(i,j,3)]);
        LI(i,j,:) = LI(i,j,:)/norm([LI(i,j,1),LI(i,j,2),LI(i,j,3)]);
    end
end


% Convert the image to grey, if not grey already
if size(RightImage,3)==3
    GRI = double(rgb2gray(RightImage));
else
    GRI = double(RightImage);
end
if size(LeftImage,3)==3
    GLI = double(rgb2gray(LeftImage));
else
    GLI = double(LeftImage);
end
GRI=imresize(GRI,scale); GLI=imresize(GLI,scale);

%hold on;
%figure;  imshow(uint8(GLI)); %hold off;hold on;
%figure;  imshow(uint8(GRI)); %hold off;hold on;

% ============================== %

%clearvars RightImage LeftImage

%%

% Initialize the Disparity Map
DisparityMap = zeros(size(GLI));

% Initialize the Cost Matrix (D)
D = zeros(size(GRI,1),size(GRI,2),size(GRI,2));
D(:,1,:) = repmat(occlusion*(1:length(D)),size(D,1),1);
D(:,:,1) = D(:,1,:);
% Initialize the Flag Matrix (A)
A = zeros(size(D));
Match = zeros(size(A,2),2);
%%
% ===== Scanning ====== %
% iterate through each row of the two images
%for y = 1:size(GRI,1),
%%
% ===== Matching/Occluding ====== %
% Iterate through X position of Left Picture
for i = 2:length(D),
    % Iterate through X position of Right Picture
    for j = 2:length(D),
        % Assign Min Cost in D and Flag in A
        [D(:,i,j),A(:,i,j)] = min([...
            (D(:,i-1,j-1)+(GLI(:,i)-GRI(:,j)).^2 ...
            + sum(abs(cross(LI(:,i,:),RI(:,j,:))),3)),...
            D(:,i-1,j)+occlusion,...
            D(:,i,j-1)+occlusion,...
            D(:,i-1,j-1)+2*occlusion]');
    end
end
% ============================== %
            %+ sum(abs(cross(LI(:,i,:),RI(:,j,:))),3)),...

%clearvars D occlusion

% ===== Interpreting Scan ====== %
% Prepare reconstruct of the optimum match

%ij = repmat(i*j,1,size(D,1)).*remapper;
% Initialize the Match-Index Array to max length
%%
ij=i;
for y=1:size(D,1),
% Back-track through the Flag Matrix
place=1;
while i>0 && j>0,
    switch A(y,i,j)
        case 1
            % pixel i of left scan-line matches 
            % pixel j of right scan-line
            Match(place,:)=[i,j];
            place=place+1;
            i=i-1; j=j-1;
            %disp('Match');
        case 2
            % pixel i of left scan-line has no match
            i=i-1;
            %disp('Left un-Match');
        case 3
            % pixel j of right scan-line has no match
            j=j-1;
            %disp('Right un-Match');
        case 4
            % pixel i of left scan-line has no match
            % and pixel j of right scan-line has no match
            i=i-1; j=j-1;
            %disp('Left&Right un-Match');
        otherwise
            % Just a default action
            %disp('Should not happen');
            i=i-1; j=j-1;
    end
end
%%
% ============================== %

%clearvars A j
%%
% ===== Storing Scan Info ====== %

% Set the second column of the Match matrix
% to be the inverse disparity of pixel i in
% the first column of the left scanline
% Match(:,2) = -1./abs(max(Match(:,1),Match(:,2))...
%     -min(Match(:,1),Match(:,2))-1);
% Match((Match(:,2)==-1),2) = -0.25;
Match(:,2)=abs(Match(:,2)-Match(:,1));
% Set end of index list
place = place-1;
% iterate through the non-occluded 
% pixels in the left scan-line
DisparityMap(y,Match(1:place,1)) = Match(1:place,2);
% ============================== %
i=ij;j=i;
%clearvars place i Match
%disp(y);
end
% ============================== %
%DM = interp2(DisparityMap,5,'cubic');
%DisparityMap=DM;
DisparityMap=DisparityMap/(max(max(DisparityMap)));

% ===== Saving/Outputting ====== %

% Convert me relateve stereo to a
% display-able setting (Produce)
Produce = mat2gray(DisparityMap,[0 1]);

%Produce = imresize(Produce,size(GRI,1)/size(Produce,1));
%clearvars y DisparityMap

figure;
imshow(Produce);
hold on;
%%
imwrite(Produce, Map, 'jpg');
% ============================== %

%clear all;

end