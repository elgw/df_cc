% Demo of df_cc
% Two channels, a constant shift

% Number of dots to use
ndots = 10;
% The domain/image size
dom = [512, 512, 60];

displacement = 17; % displacement, pixels
pixelsize = [130, 130, 200];

%% Reference channel
C1 = rand(ndots, 3);
C1(:,1) = C1(:,1)*dom(1);
C1(:,2) = C1(:,2)*dom(2);
C1(:,3) = C1(:,3)*dom(3);

%% The other channel, where dots are shifted
delta = randn(1,3);
delta = displacement*delta/norm(delta);
C2 = C1;
for kk = 1:3
    C2(:, kk) = C2(:, kk) + delta(kk);
end

%% Prepare input for df_cc_create
D = {};
D{1,2} = [C1, C2];
D{2,1} = [C2, C1];

settings = df_cc_create('getDefaults');
settings.filename = 'demo1.cc'; % For correction coefficients 
settings.pixelsize = pixelsize;
df_cc_create('dots', D, ...
    'channels', {'red', 'green'}, ...
    'settings', settings);

% Now correct the dots in the second channel to match the first
C2_corr = df_cc_apply_dots('dots', C2, ...
    'from', 'green', ...
    'to', 'red', ...
    'ccFile', 'demo1.cc');


% 3D scatter
% sc3 = @(X, y) scatter3(X(:,1), X(:,2), X(:,3), y);
% 2D scatter
sc3 = @(X, y) scatter(X(:,1), X(:,2), y);

figure,
sc3(C1, 'ro')
hold on
sc3(C2, 'go')
sc3(C2_corr, 'gx')

for kk = 1:size(C2,1)
    a = C2(kk,:); b = C2_corr(kk,:);
    plot([a(1) b(1)], [a(2), b(2)], 'k');
end

legend({'Reference channel', 'Other channel', 'Corrected'})
axis equal
title('demo1, a shift')
