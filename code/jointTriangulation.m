function [JT1, JT2, JC] = jointTriangulation(c1, c2)
% Implementation of Aronov, Seidel and Souvaine's On Compatible Triangulations of Simple Polygons
% Initial simple triangulations of shapes represented by c1 and c2 are mapped to the unit disk.
% On the unit disk they are intersected with one another and mapped back to c1 and c2. 
% The triangles on JT1 correspond to the triangles in JT2
% DISCLAIMER: the method will fail if c1 or c2 have intersecting edges.

% Usage:   [JT1, JT2, JC] = jointTriangulation(c1, c2)
%
% Arguments:  
%               c1  - Closed contour m x 2
%               c2  - Closed contour m x 2, note: c1 and c2 have to be of the same length
%
% Returns: 
%               JT1 - Joint triangulation object on shape 1. A Matlab triangulation object. 
%               JT2 - Joint triangulation object on shape 2. A Matlab triangulation object. 
%                JC - Joint triangulation object on the unnit disk.  

% Example:
% I1 = imread('..\data\bw_fish5.png');
% I2 = imread('..\data\bw_fish8.png');
% BW1 = imbinarize(I1);
% BW2 = imbinarize(I2);
% smoothing = 5;
% sampling = 80;
% C1 = bwboundaries(BW1,4); C1 = C1{1};
% C1 = buildContour(BW1,C1,smoothing,'Sampling',sampling);
% C2 = bwboundaries(BW2,4); C2 = C2{1};
% C2 = buildContour(BW2,C2,smoothing,'Sampling',sampling); 
% [JT1, JT2, JC] = jointTriangulation(C1,C2);

% Reference: Aronov, Boris, Raimund Seidel, and Diane Souvaine. 
% "On compatible triangulations of simple polygons." 
% Computational Geometry 3.1 (1993): 27-35.

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

if(size(c1,1)~= size(c2,1))
   error('c1 and c2 have different lengths, c1 and c2 must correspond!'); 
end

boundary1 = [1:size(c1,1);circshift(1:size(c1,1),-1)]';
boundary2 = [1:size(c2,1);circshift(1:size(c2,1),-1)]';

%First two initial Delaunay Triangulations
DT1 = delaunayTriangulation(c1, boundary1);
inside1 = DT1.isInterior();
triangles1 = DT1.ConnectivityList(inside1,:);
edges1 = compute_edges(triangles1)';

DT2 = delaunayTriangulation(c2, boundary2);
inside2 = DT2.isInterior();
triangles2 = DT2.ConnectivityList(inside2,:);
edges2 = compute_edges(triangles2)';

%Map V1 and V2 onto the circle
circle = linspace(0,2*pi-2*pi/size(c1,1),size(c1,1));
VC = [cos(circle);sin(circle)]';

%Remove coincident edges, including circular border
edges1_inside = setdiff(edges1, edges2, 'rows');
edges2_inside = setdiff(edges2, edges1, 'rows');

%Compute point intesections
%In the disk
RC1 = [VC(edges1_inside(:,1),:), VC(edges1_inside(:,2),:)];
RC2 = [VC(edges2_inside(:,1),:), VC(edges2_inside(:,2),:)];
cap = lineSegmentIntersect(RC1,RC2);

%Remove points already on the circular border
epsilon = 0.0000000001;
cap.adj = cap.adj & ~(cap.d1to2 < epsilon | cap.d1to2 > 1-epsilon);

%Order the intersections on edges1
edges1_intersect = cell(size(cap.adj,1),1);
nb_intersect1 = 0;
for i = 1:length(cap.adj)
    p0 = VC(edges1_inside(i,1),:);
    edge = cap.adj(i,:);
    nb_intersect1 = nb_intersect1 + sum(edge); %nb of new points
    if(sum(edge) == 0)
       %error('Oops, sum(edge) == 0, lone edge (line 50).');
       continue;
    end
    edge_row = cap.row(i,:);
    edge_col = cap.col(i,:);
    vec = [edge_row(edge);edge_col(edge)]';
    distances = pdist2(p0, vec, 'Euclidean');
    [~, I] = sort(distances);
    other = find(edge); %finds the position of columns, i.e. of edges in graph 2
    edges1_intersect{i} = other(I);
end

%Give id numbers to new points
new_id = zeros(size(cap.adj));
tags = (1:nb_intersect1) + size(c1,1); %the number of initial points
new_id(cap.adj) = tags;

%Get the set of new edges from graph 1
new_edges1 = [];
for i = 1:length(edges1_intersect)
   current = edges1_inside(i,1);
   for j = 1:length(edges1_intersect{i})
       next = new_id(i,edges1_intersect{i}(j));
       new_edges1(end+1,:) = [current, next];
       current = next;
   end
   new_edges1(end+1,:) = [current, edges1_inside(i,2)];
end

%Order the intersections on edges2 (problem is symmetric)
cap.adj = cap.adj';
cap.row = cap.row';
cap.col = cap.col';
new_id = new_id';

edges2_intersect = cell(size(cap.adj,1),1);
nb_intersect2 = 0;
for i = 1:length(cap.adj)
    p0 = VC(edges2_inside(i,1),:);
    edge = cap.adj(i,:);
    nb_intersect2 = nb_intersect2 + sum(edge); %nb of new points
    if(sum(edge) == 0)
         %error('Oops, sum(edge) == 0, lone edge (line 93).');
         continue;
    end
    edge_row = cap.row(i,:);
    edge_col = cap.col(i,:);
    vec = [edge_row(edge);edge_col(edge)]';
    distances = pdist2(p0, vec, 'Euclidean');
    [~, I] = sort(distances);
    other = find(edge); %finds the position of columns, i.e. of edges in graph 2
    edges2_intersect{i} = other(I);
end

%Get the set of new edges from graph 2
new_edges2 = [];
for i = 1:length(edges2_intersect)
    current = edges2_inside(i,1);
    for j = 1:length(edges2_intersect{i})
       next = new_id(i,edges2_intersect{i}(j));
       new_edges2(end+1,:) = [current, next];
       current = next;
   end
   new_edges2(end+1,:) = [current, edges2_inside(i,2)];
end

cap.adj = cap.adj';
cap.row = cap.row';
cap.col = cap.col';
new_id = new_id';

%Find the exact position of the new intersection points
%On the circle
VC_new = zeros(nb_intersect1,2);
for i = 1:nb_intersect1
    spot = new_id == tags(i); %this could be done better?
    VC_new(i,:) = [cap.row(spot), cap.col(spot)];
end

constrain_joint = [new_edges1; new_edges2; setdiff(edges1, edges1_inside,'rows')];
V_joint = [VC; VC_new];

%%%PROJECT BACK 
%Find affine transformations per triangle
tforms1 = cell(length(triangles1),1);
for i = 1:length(triangles1)
     vo = c1(triangles1(i,:),:);
     vc = VC(triangles1(i,:),:);
     tforms1{i} = fitgeotrans(vc,vo,'affine'); 
end

tforms2 = cell(length(triangles2),1);
for i = 1:length(triangles2)
     vo = c2(triangles2(i,:),:);
     vc = VC(triangles2(i,:),:);
     tforms2{i} = fitgeotrans(vc,vo,'affine'); 
end

%Find in which triangles the intersection points are in the circle domain
find_triangles1 =  intriangleNBH(triangles1, VC, VC_new, 0.000000001);
find_triangles2 =  intriangleNBH(triangles2, VC, VC_new, 0.000000001);

%Map the intersection points back
V1_new = zeros(length(VC_new),2);
for i = 1:length(VC_new)
    tr = find(find_triangles1(i,:));
    if(isempty(tr))
       error('one point is unmatched!!!');
    end
    tf = tforms1{tr(1)};
    V1_new(i,:) = transformPointsForward(tf,[VC_new(i,1),VC_new(i,2)]); 
end

V2_new = zeros(length(VC_new),2);
for i = 1:length(VC_new)
    tr = find(find_triangles2(i,:));
    if(isempty(tr))
       error('one point is unmatched!!!');
    end
    tf = tforms2{tr(1)};
    V2_new(i,:) = transformPointsForward(tf,[VC_new(i,1),VC_new(i,2)]); 
end

%Get the total joint triangulation:
T_total = delaunayTriangulation(V_joint, constrain_joint);
triangles_total = T_total.ConnectivityList(:,:);
V1_joint = [c1; V1_new];
V2_joint = [c2; V2_new];

JT1 = triangulation(triangles_total, V1_joint);
JT2 = triangulation(triangles_total, V2_joint);
JC = triangulation(triangles_total, V_joint);
end

