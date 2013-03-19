%
% Plotting 1d advection NODAL DG w/ netcdf
% -----------------------------------------

%% Read in data
ncfilename = 'dg1d_mfill.nc';

x = nc_varget(ncfilename,'x');
t = nc_varget(ncfilename,'time');
rhoq = nc_varget(ncfilename,'rhoq');
nodes = nc_varget(ncfilename,'nodes');
init = nc_varget(ncfilename,'init');
cputime = nc_varget(ncfilename,'ctime');


ncfilename = 'dg1d_nofil.nc';

x2 = nc_varget(ncfilename,'x');
t2 = nc_varget(ncfilename,'time');
rhoq2 = nc_varget(ncfilename,'rhoq');
nodes2 = nc_varget(ncfilename,'nodes');
init2 = nc_varget(ncfilename,'init');
cputime2 = nc_varget(ncfilename,'ctime');
%% Plot ICs

figure('Color', 'white');
subplot(1,1,1);
hh1 = plot(x, init, 'LineWidth', 1);
xlabel('x'); ylabel('phi(x)');
axis([0 1 -.2 1.2]);
ht = title(sprintf('Time: %0.2f sec', t(1)));

% Plot animation
% Get figure size
pos = get(gcf, 'Position');
width = pos(3); height = pos(4);

% Preallocate data (for storing frame data)
mov = zeros(height, width, 1, length(t), 'uint8');

% Loop through by changing XData and YData
%

for id = 1:length(t)%/10
    ind = id;%*10;
    % Update graphics data. This is more efficient than recreating plots.
    set(hh1(1), 'XData', x          , 'YData', rhoq(:,ind));
    set(ht, 'String', sprintf('Time: %0.2f sec', t(ind)));

    % Get frame as an image
    f = getframe(gcf);

    % Create a colormap for the first frame. For the rest of the frames,
    % use the same colormap
    if id == 1
        [mov(:,:,1,id), map] = rgb2ind(f.cdata, 256, 'nodither');
    else
        mov(:,:,1,id) = rgb2ind(f.cdata, map, 'nodither');
    end
end


% Create animated GIF
imwrite(mov, map, 'animation.gif', 'DelayTime', 0, 'LoopCount', inf);


%% Plot exact solution and last frame
%g = @(x) sin(6*pi*x) + sin(8*pi*x);

xx = linspace(0,1,101);
g = zeros(101,1);
g(xx>0.25 & xx<0.75) = 1;
%}
figure(2)
%fplot(g,[0 1],'-.')
plot(xx,g,'-.')
hold on
ii = 2;
plot(x,rhoq(:,ii),'r')
plot(x,rhoq2(:,ii),'g')
hold off
axis([0 1 -.2 1.2])
%% Compute error
%{
num_elem = 8;
ecent = zeros(1,num_elem);
dx = 1.0/num_elem;
ecent(1) = dx/2.0;
for i = 2:num_elem
    ecent(i) = ecent(i-1) + dx;
end
err = zeros(1,num_elem);

for i = 1:num_elem
    s = @(x) dgsoln(x,A,i,nodes,length(t),dx,ecent);
    q = @(x) (g(x) - soln(x)).^2;
    err(i) = quad(q,ecent(i)-dx/2.0,ecent(i)+dx/2.0,1.0e-10);
end

error = sqrt(sum(err));
%}
%steps = length(rhoq,2);
nx = length(x);
g = zeros(nx,1);
g(x>0.25 & x<0.75) = 1;

err = zeros(1,nx);
err2 = err;

err = (rhoq(:,end)-g).^2;
err2 = (rhoq2(:,end)-g).^2;

error = sqrt((1.0/nx)*sum(err))
error2 = sqrt((1.0/nx)*sum(err2))

%% Compute min value

for ii = 1:length(t)
    % ii
    tmp = squeeze(rhoq(:,ii));
    m = min(tmp(:));
    m
end

%% Plot elements
nelem = 16;
dxel = 1.0/nelem;
ecent = zeros(nelem,1);
ecent(1) = dxel/2.0;
for ii=2:nelem
    ecent(ii) = ecent(ii-1)+dxel;
end
hold on
for i = 1:nelem
     line([ecent(i)+dxel/2 ecent(i)+dxel/2],[0 1],'LineWidth',3)
 end
hold off