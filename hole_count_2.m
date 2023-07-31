function [hole_area,hole_area_hist]=hole_count_2(tif_filename)
close all;
% tif_filename;
hole_raw = imread(tif_filename) ;
h_fig = figure;
figure

imshow(hole_raw)
img_gray=hole_raw;
img_med = medfilt2(img_gray);
% imshow(img_med)
% figure(2)
% Perform image enhancement using contrast-limited adaptive histogram equalization (CLAHE)
img_clahe = adapthisteq(img_med);
imshow(img_clahe)

%  Binarize:
I_BW = imbinarize(img_clahe,0.40);

% for #13, set 0.50


% Remove those touching border
imshow(I_BW)

% waitforbuttonpress

% figure(3)
%% 
I_BW_m = ~medfilt2(I_BW,[3,3]); % Medium Filter, get rid of pepper noise
% imshow(I_BW_m)


% figure(4)
se = strel('cube',2);
I_BW_e=imerode(I_BW_m, se);
% imshow(I_BW_e)

% figure
vessels_filled = imfill(I_BW_e, 'holes');
imshow(vessels_filled)
title('vessels filled')

% figure(5)
BWnobord = imclearborder(vessels_filled,4);
% imshow(BWnobord)
% title('Cleared Border Image')





%%
figure
bw = BWnobord;


% L = watershed(I_BW);
% Lrgb = label2rgb(L);
% % figure(02)
% imshow(Lrgb)
%
%
% % figure(03)
% imshow(imfuse(bw,Lrgb))
% axis([10 175 15 155])
% % figure(04)

bw2 = ~bwareaopen(~bw, 10);
% imshow(bw2)

D = -bwdist(~bw2);
% imshow(D,[])

Ld = watershed(D);
% imshow(label2rgb(Ld))

bw2 = bw;
bw2(Ld == 0) = 0;
% imshow(bw2)

% figure
mask = imextendedmin(D,2);
% imshowpair(bw,mask,'blend')
%%
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = bw;
bw3(Ld2 == 0) = 0;
imshow(bw3)


% Define the output filename
output_filename = ['thresholded_' tif_filename];

% Save the image to a TIFF file
imwrite(bw3, output_filename);

%%

% Count the connected area
L = bwlabeln(bw3, 8);
S = regionprops(L, 'Area');
pos = ([S.Area] <= 10000) & ([S.Area] >= 5);  % To be set the area threshold
pos_ex = ~pos;
bw2 = ismember(L, find(pos));
bw2_ex = ismember(L, find(pos_ex));



S1 = [S.Area];
S1 = S1(pos);  % Final Area and number of connected regions
hole_area = S1;

N = length(S1);  % Number
disp('Holes number:')
disp(N);
edges = [0:50:10000];
figure
histogram(S1,edges)

[N2, edges] = histcounts(S1);
hole_area_hist = N2;


% Get the center of connected areas
C = regionprops(bw2, 'Centroid');  % to be processed
C1 = [C.Centroid];
C1 = reshape(C1, 2, length(C1)/2)';

% For exception
C_ex = regionprops(bw2_ex, 'Centroid');  % to be processed
C1_ex = [C_ex.Centroid];
C1_ex = reshape(C1_ex, 2, length(C1_ex)/2)';


% Mark the connected region on the orignal picture
figure(h_fig); hold on;
plot(C1(:,1), C1(:,2), 'r+', 'MarkerSize', 10);
plot(C1_ex(:,1), C1_ex(:,2), 'g+', 'MarkerSize', 10);
hold off;

xls_filename = [tif_filename 'holes_area.xlsx'];

% xls_filename='holes_area.xlsx';
% A={'holes_area';S1};

xlswrite(xls_filename,S1');

end