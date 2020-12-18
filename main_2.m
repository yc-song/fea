%% Problem 2
%역학과 설계 프로젝트- 기계항공공학부 2017-11416 송종현
%https://github.com/yc-song/fea_practice
clc; clear; fclose all;
% 배열형 데이터 설정/Material Property 설정

% 모듈 시작
if ~exist('run_data','var')
    E = 5E4; v = 0.25;
    C = E/((1+v)*(1-2*v))*[1-v v v 0; v 1-v v 0; v v 1-v 0; 0 0 0 (1-2*v)/2];
    % 수치해 구하기
    [node, element] = geoNmesh;
    % Boundary Condition(Constraint)
    Con_u = 2*find(node(2,:) > -1E-9 & node(2,:) < 1E-9);
    ndof = setdiff(1:2*length(node),Con_u);
    % Loading Condition
    f = [100; 0];
    u_L = find(node(1,:) > 50-1E-9 & node(1,:) < 50+1E-9);
    % Global K Matrix
    K = Kmat(node, element, C);
    % Global F Vector
    F = force(node, element, f, u_L);
    % Calculate Global U
    U = zeros(2*length(node),1);
    U(ndof) = K(ndof,ndof)\F(ndof);
    % Calculate Stress Vector
   St_r = stress_rr(node, element, U, C);
   St_t = stress_tt(node,element,U,C);


    % 분석해 구하기
    y_ana_rr=zeros(51,1);
    y_ana_tt=zeros(51,1);
    x1=zeros(50,1);
    for j=1:51
        x1(j)=49+j;
        
    end
    for i = 1:51
        y_ana_rr(i)=50^2*100/(100^2-50^2)*(1-100^2/(49+i)^2);
        y_ana_tt(i)=50^2*100/(100^2-50^2)*(1+100^2/(49+i)^2);
    end
        x=zeros(1,size(element,2));
    for i = 1:size(element,2)
        x(i)=mean([node(1,element(1,i)), node(1,element(2,i)), node(1,element(3,i))]);
    end
%플랏팅
    figure(1);   
    scatter(x,St_r)
    hold on
    plot(x1,y_ana_rr)
    title('FEM v. Analytical Solution (The Radial Stress)')
    xlabel('length(mm)');
    ylabel('Stress (N/mm^2)');
    legend('FEM','Analytical')
    hold off
    figure(2);   
    scatter(x,St_t)
    hold on
    plot(x1,y_ana_tt)
    title('FEM v. Analytical Solution (The Hoop Stress)')
    xlabel('length(mm)');
    ylabel('Stress (N/mm^2)');
    legend('FEM','Analytical')
    hold off
    

end

function [KG] = Kmat(node, element, C)
IK = []; JK = []; K = [];
for i = 1 : size(element,2)
  % Construct the Local Ke
  X = [node(1,element(1,i)) node(1,element(2,i)) node(1,element(3,i))]';
  Y = [node(2,element(1,i)) node(2,element(2,i)) node(2,element(3,i))]';
  S = [1 X(1) Y(1); 1 X(2) Y(2); 1 X(3) Y(3)];
  r_avg = mean(X);
  z_avg = mean(Y);
  alpha = [[1 0 0]*(S\[1 0 0]'); [1 0 0]*(S\[0 1 0]'); [1 0 0]*(S\[0 0 1]')];
  A = (X(1)*(Y(2)-Y(3))+X(2)*(Y(3)-Y(1))+X(3)*(Y(1)-Y(2)))/2;
  dhdx = [[0 1 0]*(S\[1 0 0]'); [0 1 0]*(S\[0 1 0]'); [0 1 0]*(S\[0 0 1]')];
  dhdy = [[0 0 1]*(S\[1 0 0]'); [0 0 1]*(S\[0 1 0]'); [0 0 1]*(S\[0 0 1]')];
  Be = [dhdx(1) 0 dhdx(2) 0 dhdx(3) 0;
    0 dhdy(1) 0 dhdy(2) 0 dhdy(3);
    (alpha(1)/r_avg+dhdx(1)+dhdy(1)*z_avg/r_avg) 0 (alpha(2)/r_avg+dhdx(2)+dhdy(2)*z_avg/r_avg) 0 (alpha(3)/r_avg+dhdx(3)+dhdy(3)*z_avg/r_avg) 0;
    dhdy(1) dhdx(1) dhdy(2) dhdx(2) dhdy(3) dhdx(3)];
  Ke =  2*pi*Be'*C*Be*r_avg*A;
  % Change to Global Nod num.
  IK = [IK; repmat([2*element(1,i)-1 2*element(1,i) 2*element(2,i)-1 2*element(2,i)...
    2*element(3,i)-1 2*element(3,i)]',6,1)];
  JK = [JK repmat([2*element(1,i)-1 2*element(1,i) 2*element(2,i)-1 2*element(2,i)...
    2*element(3,i)-1 2*element(3,i)],6,1)];
  K = [K; Ke(:)];
end
KG = sparse(IK(:), JK(:), K(:));
end

function [St_r] = stress_rr(node, element, U, C)
syms x y
St_r = [];

for i = 1 : size(element,2)
    % Stiffness Matrix와 같은 방식으로 접근
  X = [node(1,element(1,i)) node(1,element(2,i)) node(1,element(3,i))]';
  Y = [node(2,element(1,i)) node(2,element(2,i)) node(2,element(3,i))]';
  S = [1 X(1) Y(1); 1 X(2) Y(2); 1 X(3) Y(3)];
  r_avg = mean(X);
  z_avg = mean(Y);
  dhdx = [[0 1 0]*(S\[1 0 0]'); [0 1 0]*(S\[0 1 0]'); [0 1 0]*(S\[0 0 1]')];
  dhdy = [[0 0 1]*(S\[1 0 0]'); [0 0 1]*(S\[0 1 0]'); [0 0 1]*(S\[0 0 1]')];
  alpha = [[1 0 0]*(S\[1 0 0]'); [1 0 0]*(S\[0 1 0]'); [1 0 0]*(S\[0 0 1]')];
  %{B}
 Be = [dhdx(1) 0 dhdx(2) 0 dhdx(3) 0; 0 dhdy(1) 0 dhdy(2) 0 dhdy(3); alpha(1)/r_avg+dhdx(1)+dhdy(1)*z_avg/r_avg 0 alpha(2)/r_avg+dhdx(2)+dhdy(2)*z_avg/r_avg 0 alpha(3)/r_avg+dhdx(3)+dhdy(3)*z_avg/r_avg 0; dhdy(1) dhdx(1) dhdy(2) dhdx(2) dhdy(3) dhdx(3)];
 U_Local = [U(2*element(1,i)-1,1) U(2*element(1,i),1) U(2*element(2,i)-1,1) U(2*element(2,i),1) U(2*element(3,i)-1,1) U(2*element(3,i),1)];
 s_Local = C*Be*U_Local'; 
 St_r = [St_r;s_Local(1,1)];
end
end

function [St_t] = stress_tt(node, element, U, C)
syms x y
St_t=[];
for i = 1 : size(element,2)
    % Stiffness Matrix와 같은 방식으로 접근
  X = [node(1,element(1,i)) node(1,element(2,i)) node(1,element(3,i))]';
  Y = [node(2,element(1,i)) node(2,element(2,i)) node(2,element(3,i))]';
  S = [1 X(1) Y(1); 1 X(2) Y(2); 1 X(3) Y(3)];
  r_avg = mean(X);
  z_avg = mean(Y);
  dhdx = [[0 1 0]*(S\[1 0 0]'); [0 1 0]*(S\[0 1 0]'); [0 1 0]*(S\[0 0 1]')];
  dhdy = [[0 0 1]*(S\[1 0 0]'); [0 0 1]*(S\[0 1 0]'); [0 0 1]*(S\[0 0 1]')];
  alpha = [[1 0 0]*(S\[1 0 0]'); [1 0 0]*(S\[0 1 0]'); [1 0 0]*(S\[0 0 1]')];
  %{B}
 Be = [dhdx(1) 0 dhdx(2) 0 dhdx(3) 0; 0 dhdy(1) 0 dhdy(2) 0 dhdy(3); alpha(1)/r_avg+dhdx(1)+dhdy(1)*z_avg/r_avg 0 alpha(2)/r_avg+dhdx(2)+dhdy(2)*z_avg/r_avg 0 alpha(3)/r_avg+dhdx(3)+dhdy(3)*z_avg/r_avg 0; dhdy(1) dhdx(1) dhdy(2) dhdx(2) dhdy(3) dhdx(3)];
U_Local = [U(2*element(1,i)-1,1) U(2*element(1,i),1) U(2*element(2,i)-1,1) U(2*element(2,i),1) U(2*element(3,i)-1,1) U(2*element(3,i),1)];
s_Local = C*Be*U_Local'; 
St_t = [St_t;s_Local(3,1)];
end
end

function [node, element] = geoNmesh()
symModel = createpde('structural','static-axisymmetric');
R = 100.0;
r = 50.0;
z = 1000;
R2 = [3,4, [r,R,R,r,z,z,0,0]]';
gdm = [R2];
ns = char('R2');
g = decsg(gdm,'R2',ns');
geometryFromEdges(symModel,g);
structuralProperties(symModel,'YoungsModulus',5E4, 'PoissonsRatio',0.25);
structuralBC(symModel,'Edge',[1 3],'Constraint','symmetric'); %symmetric하게 y,x에 대하여 고정
structuralBoundaryLoad(symModel,'Edge',4,'SurfaceTraction',[100;0]); %surface traction
generateMesh(symModel,'Hmax',5,'GeometricOrder','linear');
node = symModel.Mesh.Nodes;
element = symModel.Mesh.Elements;
end

function [FG] = force(nod, ele, f, u_L)
syms x y
F = []; IF = [];
for i = 1 : size(ele,2)
  X = [nod(1,ele(1,i)) nod(1,ele(2,i)) nod(1,ele(3,i))]';
  Y = [nod(2,ele(1,i)) nod(2,ele(2,i)) nod(2,ele(3,i))]';
  S = [1 X(1) Y(1); 1 X(2) Y(2); 1 X(3) Y(3)];
  h = [[1 x y]*(S\[1 0 0]'); [1 x y]*(S\[0 1 0]'); [1 x y]*(S\[0 0 1]')];
  % {F}
  [D,ia] = intersect(ele(:,i),u_L');
  if length(ia) == 2
    fe = (2*pi*50)*(eval(int(subs([h(1) 0 h(2) 0 h(3) 0; 0 h(1) 0 h(2) 0 h(3)],x,50),y,min(nod(2,D)), max(nod(2,D))))'*f);
  else
    fe = zeros(6,1);
  end
  IF = [IF; [2*ele(1,i)-1 2*ele(1,i) 2*ele(2,i)-1 2*ele(2,i)...
    2*ele(3,i)-1 2*ele(3,i)]'];
  F = [F; fe];
end
% Construct Global Force Vector
FG = sparse(IF(:),1,F(:));
end
