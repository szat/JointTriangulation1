% Copyright (c) Adrian Szatmari
% Author: Adrian Szatmari
% Date: 2017-11-30
% License: MIT, patent permitting
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% Implementation of Aronov, Seidel and Souvaine's On Compatible Triangulations of Simple Polygons
% Initial simple triangulations of shapes represented by c1 and c2 are mapped to the unit disk.
% On the unit disk they are intersected with one another and mapped back to c1 and c2. 
% The triangles on JT1 correspond to the triangles in JT2

% Load Shapes
I1 = imread('..\data\bw_fish4.png');
I2 = imread('..\data\bw_fish6.png');

BW1 = imbinarize(I1);
BW2 = imbinarize(I2);

% Build Contours
smoothing = 5;
sampling = 80;

C1 = bwboundaries(BW1,4); C1 = C1{1};
C1 = buildContour(BW1,C1,smoothing,'Sampling',sampling);

C2 = bwboundaries(BW2,4); C2 = C2{1};
C2 = buildContour(BW2,C2,smoothing,'Sampling',sampling);

%Check no self intersection
constrain1 = [1:length(C1);circshift(1:length(C1),-1)]';
constrain2 = [1:length(C2);circshift(1:length(C2),-1)]';
delaunayTriangulation(C1,constrain1);
delaunayTriangulation(C2,constrain2);

% jointTriangulation will fail if there are intersecting edges in either C1 or C2
[JT1, JT2, JC] = jointTriangulation(C1,C2);

% Visualize
figure
triplot(JT1)
set(gca,'Xdir','reverse')
view([90 -90])
title('Joint Triangulation on Shape 1');

figure
triplot(JC)
set(gca,'Xdir','reverse')
view([90 -90])
title('Joint Triangulation on the Unit Disk');

figure
triplot(JT2)
set(gca,'Xdir','reverse')
view([90 -90])
title('Joint Triangulation on Shape 2');