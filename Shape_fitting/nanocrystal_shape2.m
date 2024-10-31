function vertices_output = nanocrystal_shape(plane_info,b)

%This code will compute the vertices that define the shape of a nanocrystal
%based on two inputs:
% plane_info: [h, k, l, multiplier], can contain multiple rows to define multiple planes
% b: Vectors that define the crystal lattice in [x1, y1, z1; x2, y2, z2; x3, y3, z3];

% Retrive plane properties from plane_info
g = plane_info(:,1:3);
% Retrieve multiplier from plane_info
d_spacings = ((plane_info(:,1).^2)/b(1,1).^2+(plane_info(:,2).^2)/b(2,2).^2+(plane_info(:,3).^2)/b(3,3).^2).^(-1/2)'; %in A
d_spacings = d_spacings*.1; %convert to nm
d = plane_info(:,4)'; %Distance in nm


%% Define reciprocal vectors
recip(1,:) = 2*pi*cross(b(2,:),b(3,:))./dot(b(1,:),cross(b(2,:),b(3,:)));
recip(2,:) = 2*pi*cross(b(3,:),b(1,:))./dot(b(2,:),cross(b(3,:),b(1,:)));
recip(3,:) = 2*pi*cross(b(1,:),b(2,:))./dot(b(3,:),cross(b(1,:),b(2,:)));

%Extend hkl along the vectors that define the lattice and apply the multiplier
%n_d = (g*b).*d'; %Extend based on spacing
n_d = g(:,1).*recip(1,:) +  g(:,2).*recip(2,:) +  g(:,3).*recip(3,:); %Change lengths
n_l = (n_d(:,1).^2 + n_d(:,2).^2 + n_d(:,3).^2).^(1/2); %length of planes
%%


n_norm = n_d./vecnorm(n_d')'; %Normal vectors for a*x + b*y+c*z = s
n_d_points = n_norm.*d'; %Points on plane

%Determine the normal vector (length of one)
%Determine how many planes you have, define variable for easier indexing
planes_index = 1:size(n_d,1);

%Obtain 's' value used to define plane for a*x + b*y+c*z = s
planes_s = d'; % s for a*x + b*y+c*z = s


% Find all plane combinations
planes_combinations = nchoosek(planes_index,2);
%determine which combinations are not parallel and save them
planes_non_par = [];
for i = 1:size(planes_combinations,1)
par_check(i) = abs(dot(n_norm(planes_combinations(i,1),:),n_norm(planes_combinations(i,2),:))./(norm(n_norm(planes_combinations(i,1),:))*norm(n_norm(planes_combinations(i,2),:)))); % Determine if planes are parallel
if round(par_check(i),4) ~= 1 %checks if planes are parallel
planes_non_par = [planes_non_par; planes_combinations(i,:)]; %Create list of non-parallel planes
end
end

% Find lines formed by plane pairs using cross product
for i = 1:size(planes_non_par,1)
[p, n] = intersect_planes(n_d_points(planes_non_par(i,1),:), n_norm(planes_non_par(i,1),:), n_d_points(planes_non_par(i,2),:), n_norm(planes_non_par(i,2),:));
lines_p(i,:) = round(p,10);
lines_n(i,:) = n;
end

%% Create plane-line combinations to determine vertices

%Use meshgrid to obtain all combinations
[L,P] = meshgrid(1:size(lines_n,1),planes_index);

lp_list = [L(:),P(:)]; %List of plane-line pairs

%Find line-plane pairs that are no parallel
lp_list_nonpar = [];
for i = 1:size(lp_list,1)
lp_par_check(i) = round(dot(lines_n(lp_list(i,1),:),n_norm(lp_list(i,2),:)),4);
if lp_par_check(i) ~= 0
lp_list_nonpar = [lp_list_nonpar; lp_list(i,:)];
end
end

%Find parametric value for line that intersects the plane, and save this value
for i = 1:size(lp_list_nonpar,1)
d_line(i,1) = dot(n_d_points(lp_list_nonpar(i,2),:)-lines_p(lp_list_nonpar(i,1),:), n_norm(lp_list_nonpar(i,2),:))/dot(lines_n(lp_list_nonpar(i,1),:),n_norm(lp_list_nonpar(i,2),:));
vertices_values(i,:) = lines_p(lp_list_nonpar(i,1),:)' +d_line(i,1)*lines_n(lp_list_nonpar(i,1),:)';
end
vertices_values = round(vertices_values,6); %Change this value to change accuracy, default 4
vertices_values = unique(vertices_values,'rows'); %Only keep unique vertices

% Compute s_values for vertices

%Use vertix and plane pairs to determine if vertices are within the
%nanocrystal, this assumes that the origin is always within the nanocrystal

[V,P] = meshgrid(1:size(vertices_values,1),planes_index);

vp_list = [V(:),P(:)];
v_list = [];
for i = 1:size(vp_list,1)
side_value(i,1) = round(dot(vertices_values(vp_list(i,1),:),n_norm(vp_list(i,2),:))-planes_s(vp_list(i,2)),4); %Don't set limit too high, 4 is good, 10 results in errors for non-integer planes
if side_value(i,1) > 0
v_list = [v_list;vp_list(i,1)]; %These are the vertices to exclude
end
end
vertices_final = [];

for f = 1:size(vertices_values,1)
if ismember(f,v_list)
else
vertices_final = [vertices_final; vertices_values(f,:)];
end
end

vertices_output = vertices_final; %Output the vertices that were found

end

% Function obtained from MATLAB community post to find line at intersection
function [p, n] = intersect_planes(p1, n1, p2, n2)
% p1 = point on plane 1
% n1 = normal vector of plane 1
% p2 = point on plane 2
% n2 = normal vector of plane 2
p0 = [0,0,0]; %Point closest to origin

M = [2 0 0 n1(1) n2(1)
     0 2 0 n1(2) n2(2)
     0 0 2 n1(3) n2(3)
     n1(1) n1(2) n1(3) 0 0
     n2(1) n2(2) n2(3) 0 0
    ];

b4 = p1(1).*n1(1) + p1(2).*n1(2) + p1(3).*n1(3);
b5 = p2(1).*n2(1) + p2(2).*n2(2) + p2(3).*n2(3);
b = [2*p0(1) ; 2*p0(2) ; 2*p0(3); b4 ; b5];

x = M\b;
p = x(1:3)';
n = cross(n1, n2);
end