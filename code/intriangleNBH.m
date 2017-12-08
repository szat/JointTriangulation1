function [TABLE] = intriangleNBH(triangles,vertices,queries,epsilon)
% Finds which query belongs to which triangle, within epsilon accuracy.
% Therefore a query might belong to more than one triangle, for instance if query is 
% on a vertex or on an edge of a triangle. 

% Usage: [TABLE] = intriangleNBH(triangles, vertices, queries, epsilon)
%
% Arguments:  
%     triangles  - m x 3, where triangles(i,:) are the vertex indices of triangle i.
%      vertices  - n x 2, where vertices(j,:) are the coordinates of vertex j. 
%       queries  - k z 2, where queries(k,:) are the coordinates of query k.
%       epsilon  - radius of the neighborhood in which to search for query. 
%
% Returns: 
%          TABLE - k x m, TABLE(i,j) == 1 if query i is in (or almost in) triangle j. 

% Copyright (c) Adrian Szatmari
% Author: Adrian Szatmari
% Date: 2017-12-06
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

e1 = [0,epsilon];
e2 = [0,-epsilon];
e3 = [epsilon,0];
e4 = [-epsilon,0];
Q1 = queries + e1;
Q2 = queries + e2;
Q3 = queries + e3;
Q4 = queries + e4;

T1 = intriangle(triangles, vertices, Q1);
T2 = intriangle(triangles, vertices, Q2);
T3 = intriangle(triangles, vertices, Q3);
T4 = intriangle(triangles, vertices, Q4);

TABLE = T1 | T2 | T3 | T4;
end

function [TABLE] = intriangle(triangles,vertices,queries)
nb_triangles = size(triangles,1);
nb_points = size(queries,1);

%INITIALIZATION
X1 = vertices(triangles(:,1),1)';
X2 = vertices(triangles(:,2),1)';
X3 = vertices(triangles(:,3),1)';
Y1 = vertices(triangles(:,1),2)';
Y2 = vertices(triangles(:,2),2)';
Y3 = vertices(triangles(:,3),2)';
X = queries(:,1);
Y = queries(:,2);

%SETTING MATRICES
XX1 = repmat(X1,nb_points,1);
XX2 = repmat(X2,nb_points,1);
XX3 = repmat(X3,nb_points,1);
YY1 = repmat(Y1,nb_points,1);
YY2 = repmat(Y2,nb_points,1);
YY3 = repmat(Y3,nb_points,1);
XX = repmat(X,1,nb_triangles);
YY = repmat(Y,1,nb_triangles);

%COMPUTATIONS, inspired from
%https://www.mathworks.com/matlabcentral/answers/277984-check-points-inside-triangle-or-on-edge-with-example

S = (X1 - X2).*(Y3-Y1) - (Y1-Y2).*(X3-X1);
S = repmat(S, nb_points,1);

DET1 = (XX3-XX).*(YY2-YY3) - (YY3-YY).*(XX2-XX3);
DET2 = (XX1-XX).*(YY3-YY1) - (YY1-YY).*(XX3-XX1);
DET3 = (XX2-XX).*(YY1-YY2) - (YY2-YY).*(XX1-XX2);

TABLE = (S.*DET1 >= 0) & (S.*DET2 >= 0) & (S.*DET3 >= 0);
end