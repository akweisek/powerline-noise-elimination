%%   Removal of Powerline Noise via Blind Source Separation and Wavelet Analysis
%    Copyright (C) 2015  Samuel Akwei-Sekyere
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

%%
%
%Input : Data = Single row of neural data+AC noise, Sample_rate= sampling rate (in Hz) of recorded data, freq=the AC frequency you want to remove, thresh=the threshold within AC to look into window= window size (in milliseconds eg. 100000), wavlib = optional, pre-made wavelet library
%Output: Signal = raw neural data extracted , AC_noise = AC_noise extracted
%
%
% To create the wavelet library (optional, because it will still create it for you anyway);
%    wavelet_library=wav_m(scales,translations,time);   .... You may use the code below to know how to choose the scales and translations.
%
% This code shows how the algorithm works; it might not be the most computationally efficient way of running the code.

%

function [Signal,AC_noise]=power_noise_remove(Data,Sample_rate,freq,thresh,window,wavlib)
%
%window in ms
win_width=round(window./1000*Sample_rate);

Signal=[];
AC_noise=[];

%%%%%%%%%%%%%%%%%%%%%
%Initializing ... will need this later
time=(0:1:(win_width-1))/Sample_rate;
time=time-(mean(time));
scales=30:0.5:80;
translations=time;
real_morlet_wavelet=real_morlet(translations,freq,2);

if nargin < 6
    wavlib=wav_m(scales,translations,time);   %obtaining the wavelet library
end

%%
h=waitbar(0,'Please wait for the denoising process... \newline \newlineCopyright (C) 2015  Samuel Akwei-Sekyere');
s=length(Data);

b=1;
%%%%%%%%%%%%%%%%%%%%%%%%
while b < length(Data)
    if length(Data(b:end))>=win_width
        [signal_n,ac] = extract_single_in(Data(b:b+win_width-1),wavlib,scales,freq,thresh,real_morlet_wavelet,translations);  %Calculation for each window
        AC_noise=horzcat(AC_noise,ac(:,:));
        Signal=horzcat(Signal,signal_n);
               
        
    else
        time=(0:1:length(Data(b:end))-1)/Sample_rate;
        time=time-(mean(time));
        scales=30:0.5:80;
        translations=time;
        wavlib=wav_m(scales,translations,time);   %obtaining the wavelet library
        real_morlet_wavelet=real_morlet(translations,freq,2);
        
        [signal,ac] = extract_single_in(Data(b:end),wavlib,scales,freq,thresh,real_morlet_wavelet,translations);
        ac=ac(:,:);
        
        AC_noise=horzcat(AC_noise,ac(:,:));
        Signal=horzcat(Signal,signal);
        
        
    end

   b=b+win_width;
   
   waitbar(b/s)
end

close(h)
   
   waitbar(1,'Completed! \newline \newlineCopyright (C) 2015  Samuel Akwei-Sekyere')

%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%
function [signal,ac_noise,Decomp]=extract_single_in(Data,wavlib,scales,freq,thresh,real_morlet_wavelet,translation)


%%Decomposition of the Data by EEMD

Decomp=eemd(Data,0.05,100,10);   %Call EEMD function ; number of IMFs = 10
%%Detecting Noise(s)

signal_decomp=[];
noise_decomp=[];



%%

for i=1:length(Decomp(:,1))
	[Xtic]=Scalogram_Prob(Decomp(i,:),wavlib,scales,freq,thresh);
	
    if strcmp(Xtic,'signal')
    		signal_decomp=vertcat(signal_decomp,Decomp(i,:));	
    elseif strcmp(Xtic,'ac_noise')
            noise_decomp=vertcat(noise_decomp,Decomp(i,:));
    else
            continue
    end
 end


if ~isempty(noise_decomp)
    
    [amplitude_correction]=find_closest_amplitude_correction(Data,wavlib,scales,freq);
    [~,ac_noise_decomp_n]=cubica34(noise_decomp);
    [~,~,scalogram]=find_closest_ac(ac_noise_decomp_n,wavlib,scales,freq);
    
    for g=1:length(Data)
       ac_noise(1,g)=sum(scalogram.*real_morlet((translation-translation(g)),freq,2));
    end
    
    c=(sum(abs(real_morlet_wavelet))).^2;
    %a=((2.^(3/4)));
    ac_noise=(1./c).*(rms(amplitude_correction)./rms(scalogram)).*(ac_noise);   %Starting point...
    
	%amps=0.01:0.001:3;
    
    %least_sq=zeros(length(amps),length(ac_noise));
    %cost=zeros(1,length(amps));
    %for r = 1:length(amps)
    %    least_sq(r,:) = amps(r).*ac_noise;
    %    cost(1,r)=sum((least_sq(r,:)-Data).^2);
    %end
    
    %[~,indices]=min(cost);
    %factor=amps(1,indices);
	
	%Alternatively, after solving the equation
	factor=sum(ac_noise.*Data)/sum(ac_noise.^2);
    
    ac_noise=factor.*ac_noise;
    
else
    ac_noise=zeros(size(Decomp(1,:)));
end



signal=Data-ac_noise;


%%%%%%%%%%%%%%---------------------------------------------------------------------------------

function [Xtic]=Scalogram_Prob(Data,wavlib,scales,freq,thresh)

Scalogram=abs(cwt_extract(Data,wavlib));

Scalogram=abs(Scalogram-mean2(Scalogram));
Scalogram=sum((Scalogram'));

[~,ind]=max(Scalogram);

find_freq=scales(ind);

if freq-thresh<find_freq && find_freq<freq+thresh
    Xtic='ac_noise';
else
    Xtic='signal';
end
%%%%%%%%%%%%----------------------------------------------------------------------------------

function [ac_noise,index_choose,scalogram_reconstruct]=find_closest_ac(input,wavlib,scales,freq)

for p=1:length(input(:,1))
    scalogram=cwt_extract(input(p,:),wavlib);
    scalogram_reconstruct=scalogram;
    scalogram=sum(((abs(abs(scalogram)-mean2(abs(scalogram))))'));
    [~,ind]=max(scalogram);
    freq_prob(p,:)=scales(ind);
end
to_sort=freq_prob-freq;
[~,index]=sort(abs(to_sort));

ac_noise=input(index(1),:);
index_choose=find(scales==freq);
scalogram_reconstruct=scalogram_reconstruct(index_choose,:);

%%%%%%%%%%%%%%--------------------------------------------------------------------------------------------
function [scalogram_reconstruct]=find_closest_amplitude_correction(input,wavlib,scales,freq)

scalogram_reconstruct=cwt_extract(input,wavlib);
index_choose=find(scales==freq);
scalogram_reconstruct=scalogram_reconstruct(index_choose,:);


%%
function  [R,y]=cubica34(x)
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