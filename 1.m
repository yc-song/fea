%%
model = createpde('structural','static-planestress');
radius = 20.0;
width = 50.0;
totalLength = 4*width;
R1 = [3 4 0  totalLength ...
           totalLength 0 ...
          0 0 width width]'; 
C1 = [1 0 0 radius 0 0 0 0 0 0]';
gdm = [R1 C1];
ns = char('R1','C1');
g = decsg(gdm,'R1-C1',ns');
geometryFromEdges(model,g);
figure(1);
pdegplot(model,'EdgeLabel','on');
axis equal
title 'Geometry with Edge Labels';
structuralProperties(model,'YoungsModulus',200E3,'PoissonsRatio',0.25);
structuralBC(model,'Edge',[3,4],'Constraint','symmetric'); %symmetric하게 y,x에 대하여 고정
structuralBoundaryLoad(model,'Edge',1,'SurfaceTraction',[100;0]); %surface traction
%%
generateMesh(model,'Hmax',radius/6); % Max. element edge size. geometry order-quadratic (6개의 node를 가진 triangle)
figure(2);
pdemesh(model); % plot
R=solve(model);
%%
figure(3);
pdeplot(model,'XYData',R.Strain.sxx,'ColorMap','jet') %stress 대신 displacement, von mises등 가능
axis equal
title 'Normal Stress Along x-Direction';
