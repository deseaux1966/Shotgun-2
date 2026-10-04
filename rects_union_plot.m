function rects_union_plot(ax, X, Y, color)
% Draw the union of axis-aligned rectangles (columns of X,Y: 4 corners each) as
% one seamless polygon, avoiding anti-aliasing hairlines between abutting patches.
if isempty(X), return, end
e = 1e-9;
n = size(X,2);
ps(n) = polyshape();
for k = 1:n
    ps(k) = polyshape([min(X(:,k))-e max(X(:,k))+e max(X(:,k))+e min(X(:,k))-e], ...
                      [min(Y(:,k))-e min(Y(:,k))-e max(Y(:,k))+e max(Y(:,k))+e]);
end
P = union(ps);
plot(ax, P, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', 1);
end
