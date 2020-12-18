%역학과 설계 프로젝트- 기계항공공학부 2017-11416 송종현
%https://github.com/yc-song/fea_practice
%% problem #1
clc; clear; fclose all;
% 배열형 데이터 설정/Material Property 설정
y = zeros(9,1); x = zeros(9,1);
E = 7E4; v = 0.25;
C = E/(1-v^2)*[1 v 0; v 1 0; 0 0 (1-v)/2];
St_avg = 25;
i=1;
% 모듈 시작
if ~exist('run_data','var')
    for radius = 1:1:9
        [node, element] = geoNmesh(radius);
        % Boundary Condition(Constraint)
        Con_ux = 2*find(node(1,:) > -1E-9 & node(1,:) < 1E-9)-1;
        Con_uy = 2*find(node(2,:) > -1E-9 & node(2,:) < 1E-9);
        Con_u  = [Con_ux Con_uy];
        ndof = setdiff(1:2*length(node),Con_u);
        % Loading Condition
        f = [0; 25];
        u_L = find(node(2,:) > 28-1E-9 & node(2,:) < 28+1E-9);
        % Global K Matrix
        K = Kmat(node, element, C);
        % Global F Vector
        F = force(node, element, f, u_L);
        % Calculate Global U
        U = zeros(2*length(node),1);
        U(ndof) = K(ndof,ndof)\F(ndof);
        % Calculate Stress Vector/Stress Concentration Factor
        St_G = stress(node, element, U, C);
        St_max = max(St_G);
        y(i) = St_max / St_avg;
        x(i) = radius/(20-2*radius);
        i=i+1;
    end
    figure(1);
    plot(x,y,'-o','MarkerIndices',1:1:length(y));
    title('The stress concentration factor against \rho /d')
    xlabel('\rho / d')
    ylabel('Stress Concentration Factor')
end

function [KG] = Kmat(node, element, C)
IK = []; JK = []; K = [];
for i = 1 : size(element,2)
  % Construct the Local Ke
  X = [node(1,element(1,i)) node(1,element(2,i)) node(1,element(3,i))]';
  Y = [node(2,element(1,i)) node(2,element(2,i)) node(2,element(3,i))]';
  S = [1 X(1) Y(1); 1 X(2) Y(2); 1 X(3) Y(3)];
  A = (X(1)*(Y(2)-Y(3))+X(2)*(Y(3)-Y(1))+X(3)*(Y(1)-Y(2)))/2;
  dhdx = [[0 1 0]*(S\[1 0 0]'); [0 1 0]*(S\[0 1 0]'); [0 1 0]*(S\[0 0 1]')];
  dhdy = [[0 0 1]*(S\[1 0 0]'); [0 0 1]*(S\[0 1 0]'); [0 0 1]*(S\[0 0 1]')];
  Be = [dhdx(1) 0 dhdx(2) 0 dhdx(3) 0;
    0 dhdy(1) 0 dhdy(2) 0 dhdy(3);
    dhdy(1) dhdx(1) dhdy(2) dhdx(2) dhdy(3) dhdx(3)];
  Ke = Be'*C*Be*A;
  % Change to Global Nod num.
  IK = [IK; repmat([2*element(1,i)-1 2*element(1,i) 2*element(2,i)-1 2*element(2,i)...
    2*element(3,i)-1 2*element(3,i)]',6,1)];
  JK = [JK repmat([2*element(1,i)-1 2*element(1,i) 2*element(2,i)-1 2*element(2,i)...
    2*element(3,i)-1 2*element(3,i)],6,1)];
  K = [K; Ke(:)];
end
KG = sparse(IK(:), JK(:), K(:));
end

function [St_G] = stress(node, element, U, C)
syms x y
St_G=[];
for i = 1 : size(element,2)
    % Stiffness Matrix와 같은 방식으로 접근
  X = [node(1,element(1,i)) node(1,element(2,i)) node(1,element(3,i))]';
  Y = [node(2,element(1,i)) node(2,element(2,i)) node(2,element(3,i))]';
  S = [1 X(1) Y(1); 1 X(2) Y(2); 1 X(3) Y(3)];
  dhdx = [[0 1 0]*(S\[1 0 0]'); [0 1 0]*(S\[0 1 0]'); [0 1 0]*(S\[0 0 1]')];
  dhdy = [[0 0 1]*(S\[1 0 0]'); [0 0 1]*(S\[0 1 0]'); [0 0 1]*(S\[0 0 1]')];
  %{B}
  Be = [dhdx(1) 0 dhdx(2) 0 dhdx(3) 0; 0 dhdy(1) 0 dhdy(2) 0 dhdy(3); dhdy(1) dhdx(1) dhdy(2) dhdx(2) dhdy(3) dhdx(3)];
  U_Local = [U(2*element(1,i)-1,1) U(2*element(1,i),1) U(2*element(2,i)-1,1) U(2*element(2,i),1) U(2*element(3,i)-1,1) U(2*element(3,i),1)];
  s_Local = C*Be*U_Local';
  St_G = [St_G;s_Local(2,1)];
 
end
end

function [node, element] = geoNmesh(radius)
model = createpde('structural','static-planestress');
width = 10.0;
totalLength = 28.0;
R1 = [3 4 0  width ...
        width 0 ...
         0 0 totalLength totalLength]'; 
C1 = [1 0 0 radius 0 0 0 0 0 0]';
gdm = [R1 C1];
sf = 'R1-C1';
ns = char('R1','C1');
g = decsg(gdm,sf,ns');
geometryFromEdges(model,g);
structuralProperties(model,'YoungsModulus',7E4, 'PoissonsRatio',0.25);
structuralBC(model,'Edge',[3,4],'Constraint','symmetric'); %symmetric하게 y,x에 대하여 고정
structuralBoundaryLoad(model,'Edge',2,'SurfaceTraction',[0;25]); %surface traction
generateMesh(model,'Hmax',radius/6,'Hmin',0.1,'GeometricOrder','linear');
node = model.Mesh.Nodes;
element = model.Mesh.Elements;
end

function [FG] = force(node, element, f, u_L)
syms x y
F = []; IF = [];
for i = 1 : size(element,2)
  X = [node(1,element(1,i)) node(1,element(2,i)) node(1,element(3,i))]';
  Y = [node(2,element(1,i)) node(2,element(2,i)) node(2,element(3,i))]';
  S = [1 X(1) Y(1); 1 X(2) Y(2); 1 X(3) Y(3)];
  h = [[1 x y]*(S\[1 0 0]'); [1 x y]*(S\[0 1 0]'); [1 x y]*(S\[0 0 1]')];
  % {F}
  [D,ia] = intersect(element(:,i),u_L');
  if length(ia) == 2
    fe = (eval(int(subs([h(1) 0 h(2) 0 h(3) 0; 0 h(1) 0 h(2) 0 h(3)],...
      y,28),x,min(node(1,D)), max(node(1,D))))'*f);
  else
    fe = zeros(6,1);
  end
  IF = [IF; [2*element(1,i)-1 2*element(1,i) 2*element(2,i)-1 2*element(2,i) 2*element(3,i)-1 2*element(3,i)]'];
  F = [F; fe];
end
% Construct Global Force Vector
FG = sparse(IF(:),1,F(:));
end
