function d = df_cc_eudist(P,Q)
% Euclidean distance between the points in P and Q
% one point per row
% The points should have the same dimension.
% P can either be a single point or
% has to contains as many points as Q
% 
% Example:
% P = rand(10, 3);
% Q = rand(10, 3);
% D = df_cc_eudist(P,Q);
% D(1) will be norm(P(1,:) - Q(1,:))

% Need to to have the same dimensions
assert(size(P,2) == size(Q,2))
assert(size(P,1) == 1 || size(P,1) == size(Q,1));

if numel(Q) == 0
    d = [];
    return
end

if size(P,1)==1
    P = repmat(P, [size(Q,1),1]);
end

d = P-Q;
d = d.^2;
d = sum(d,2);
d = d.^(1/2);

end
