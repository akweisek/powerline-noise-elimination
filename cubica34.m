%function  [R,y]=cubica34(x)
function  [y]=cubica34(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CubICA (IMPROVED CUMULANT BASED ICA-ALGORITHM)
%
% This algorithm performes ICA by diagonalization of third- and
% fourth-order cumulants simultaneously.
%
%  [R,y]=cubica34(x)
%
% - x is and NxP matrix of observations 
%     (N: Number of components; P: Number of datapoints(samplepoints)) 
% - R is an NxN matrix such that u=R*x, and u has 
%   (approximately) independent components.
% - y is an NxP matrix of independent components
%  
% This algorithm does exactly (1+round(sqrt(N)) sweeps.
% 
% Ref: T. Blaschke and L. Wiskott, "An Improved Cumulant Based
% Method for Independent Component Analysis", Proc. ICANN-2002,
% Madrid, Spain, Aug. 27-30.
%
% questions, remarks, improvements, problems to: t.blaschke@biologie.hu-berlin.de.
%
% Copyright : Tobias Blaschke, t.blaschke@biologie.hu-berlin.de.
%
% 2002-02-22
%
%
% Last change:2003-05-19 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  [N,P]=size(x);
  Q=eye(N);
  resolution=0.001;
  
  % centering and whitening
  
 % fprintf('\ncentering and whitening!\n\n');
  
  x=x-mean(x,2)*ones(1,P);
  [V,D]=eig(x*x'/P);
  W=diag(real(diag(D).^(-0.5)))*V';
  y=W*x;
  
 % fprintf('rotating\n');
  
  % start rotating
  
  for t=1:(1+round(sqrt(N))), 
    for i=1:N-1,
      for j=i+1:N,
		
	%calculating the new cumulants
	
	u=y([i j],:);
		
	sq=u.^2;
	
	sq1=sq(1,:);
	sq2=sq(2,:);
	u1=u(1,:)';
	u2=u(2,:)';
	
	C111=sq1*u1/P;
	C112=sq1*u2/P;
	C122=sq2*u1/P;
	C222=sq2*u2/P;
	
	C1111=sq1*sq1'/P-3;
	C1112=(sq1.*u1')*u2/P;
	C1122=sq1*sq2'/P-1;
	C1222=(sq2.*u2')*u1/P;
	C2222=sq2*sq2'/P-3;
	
	% coefficients
	
	c_34=(1/6)*(1/8)*(3*(C111^2+C222^2)-9*(C112^2+C122^2)-6*(C111*C122+C112*C222));
	
	c_44=(1/24)*(1/16)*(7*(C1111^2+C2222^2)-16*(C1112^2+C1222^2)-12*(C1111*C1122+C1122*C2222)-36*C1122^2-32*C1112*C1222-2*C1111*C2222);
	
	s_34=(1/6)*(1/4)*(6*(C111*C112-C122*C222));
	
	s_44=(1/24)*(1/32)*(56*(C1111*C1112-C1222*C2222)+48*(C1112*C1122-C1122*C1222)+8*(C1111*C1222-C1112*C2222));
	
	c_48=(1/24)*(1/64)*(1*(C1111^2+C2222^2)-16*(C1112^2+C1222^2)-12*(C1111*C1122+C1122*C2222)+36*C1122^2+32*C1112*C1222+2*C1111*C2222);
	
	s_48=(1/24)*(1/64)*(8*(C1111*C1112-C1222*C2222)-48*(C1112*C1122-C1122*C1222)-8*(C1111*C1222-C1112*C2222));
	
	phi_4=-atan2(s_34+s_44,c_34+c_44);
	phi_8=-atan2(s_48,c_48);
	
	B_4=sqrt((c_34+c_44)^2+(s_34+s_44)^2);
	B_8=sqrt(c_48^2+s_48^2);
	
	%calculating the angle
	
	approx=-phi_4/4-(pi/2)*fix(-phi_4/pi);
	
	intervall=(approx-pi/8):resolution:(approx+pi/8);
	
	psi_34=B_8*cos(8*intervall+phi_8)+B_4*cos(4*intervall+phi_4);
	
	[~,index]=max(psi_34);
	
	phi_max=intervall(index);
	
	% a different way to calculate the angle is via the matlab
        % function fminbnd. The command would look like:
	%fun=[num2str(B_8),'*(-1)*cos(8*x+',num2str(phi_8),')-',num2str(B_4),
	%'*cos(4*x+',num2str(phi_4),')'];
	%phi_max=fminbnd(fun,approx-pi/8,approx+pi/8);

	
	%Givens-rotation-matrix Q_ij
	
	Q_ij=eye(N);

	c=cos(phi_max);
	s=sin(phi_max);
	
	Q_ij(i,j)=s;
	Q_ij(j,i)=-s;
	Q_ij(i,i)=c;
	Q_ij(j,j)=c;
	
	Q=Q_ij*Q;

	% rotating y
	
	y([i j],:)=[c s;-s c]*u;
	
      end %j
    end %i
  end %t
    
  R=Q*W;
  
  return