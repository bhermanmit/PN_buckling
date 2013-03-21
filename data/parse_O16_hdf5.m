function parse_O16_hdf5(hdfile,name,thermal,inelastic)

% read in capture cross section
sigc = h5read(hdfile,'/MF3/MT102/P0/dil0');
sigc = double(sigc.sigma);

% add other disapperance to capture
sigc = add_to_vector(sigc,hdfile,'/MF3/MT107/P0/dil0');

% get number of groups
ng = length(sigc);

% read in sigma scattering 0 cross section info
sigs0 = zeros(ng);
sigs0 = overwrite_matrix(sigs0,hdfile,'/MF6/MT2/P0/dil0');
if thermal ~= 0; sigs0 = overwrite_matrix(sigs0,hdfile,horzcat('/MF6/MT',num2str(thermal),'/P0/dil0')); end;
sigs1 = zeros(ng);
sigs1 = overwrite_matrix(sigs1,hdfile,'/MF6/MT2/P1/dil0');
if thermal ~= 0; sigs1 = overwrite_matrix(sigs1,hdfile,horzcat('/MF6/MT',num2str(thermal),'/P1/dil0')); end;
if inelastic == 1
    for i = 51:55
        sigs0 = add_matrix(sigs0,hdfile,horzcat('/MF6/MT',num2str(i),'/P0/dil0'));
        sigs1 = add_matrix(sigs1,hdfile,horzcat('/MF6/MT',num2str(i),'/P1/dil0'));
    end
end
% compute total cross section
sigt = zeros(ng,1);
for g = 1:ng
    sigt(g) = sigc(g) + sum(sigs0(g,:));
end

% write output file
save(name,'sigt','sigc','sigs0','sigs1');

end

function xsmat = overwrite_matrix(xsmat,hdfile,root)

% extract data from hdf5 file
data = h5read(hdfile,horzcat(root,'/meta'));
grp = data.group;
ijj = data.ijj;
njj = data.njj;
xs = double(h5read(hdfile,horzcat(root,'/sigma')));

% overall counter
kount = 1;

% begin loop around groups
for g = 1:length(grp)
    if g == 1 && grp(g) == 1 && ijj(g) == 70
        kount = kount + njj(g);
        continue
    end
    xsmat(grp(g),ijj(g):ijj(g)+njj(g)-1) = xs(kount:kount+njj(g)-1);
    kount = kount + njj(g);
    
end

end

function xsmat = add_matrix(xsmat,hdfile,root)

% extract data from hdf5 file
data = h5read(hdfile,horzcat(root,'/meta'));
grp = data.group;
ijj = data.ijj;
njj = data.njj;
xs = double(h5read(hdfile,horzcat(root,'/sigma')));

% overall counter
kount = 1;

% begin loop around groups
for g = 1:length(grp)
    if g == 1 && grp(g) == 1 && ijj(g) == 70
        kount = kount + njj(g);
        continue
    end
    xsmat(grp(g),ijj(g):ijj(g)+njj(g)-1) = ...
    xsmat(grp(g),ijj(g):ijj(g)+njj(g)-1) + xs(kount:kount+njj(g)-1);
    kount = kount + njj(g);
    
end

end

function xsvec = add_to_vector(xsvec,hdfile,root)

% extract data from hdf5 file
data = h5read(hdfile,horzcat(root));
grp = data.group;
xs = double(data.sigma);

% begin loop around groups
for g = 1:length(grp)
    xsvec(grp(g)) = ...
    xsvec(grp(g)) + xs(grp(g));    
end

end
