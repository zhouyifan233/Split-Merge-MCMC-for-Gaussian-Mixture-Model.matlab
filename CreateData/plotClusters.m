function plotClusters(X, label)
global color

if isempty(color)
    color = rand(20, 3);
end

[dim, n] = size(X);
if nargin == 1
    label = ones(n,1);
end
assert(n == length(label));

% color = 'brgmcyk';

m = size(color, 1);
c = max(label);

figure(gcf);
clf;
hold on;
axis equal;
grid on;
if dim == 2
    view(2);
    for i = 1:c
        idc = label==i;
        scatter(X(1,idc),X(2,idc), 8, color(mod(i-1,m)+1, :), 'filled');
    end
elseif dim == 3
    view(3);
    for i = 1:c
        idc = label==i;
        scatter3(X(1,idc),X(2,idc),X(3,idc),36,color(mod(i-1,m)+1), :);
    end
else
    error('ERROR: data dimension should be 2D or 3D...');
end
