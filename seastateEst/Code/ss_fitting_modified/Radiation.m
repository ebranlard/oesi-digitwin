function  ssRad = Radiation(gp,Rad,A,Ainf,B,HM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function to compute the SS Matrices                        
%                                                                 
%  Tiago Duarte, Marco Alves
%
% Tiago Duarte
% Instituto Superio Tecnico - IST Lisbon, Portugal
%
% in collaboration with:
% WavEC Offshore Renewables (www.wavec.org)
% NREL (www.wind.nrel.gov)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Re-Write matrices in the dimensional form
ldof=length(gp.DoF); dof=zeros(1); ii=0;

for d=1:1:ldof
    if gp.DoF(d)==1, ii=ii+1;
%        if d<=ldof/gp.nd
%            disp(identify(d));%Function written down below to identify each node
%        end
       dof(ii)=d; %Vector with the number of each DoF: surge=1 heave=3 ...
    end
end

% Dimensionalize the matrices

kmax=length(A(:,1))/(ldof+1);%number of frequency
AM(kmax,ii,ii)=0; %Initialize matrix
AMinf=AM; %Initialize infinite Added Mass
D=AM;%Initialize matrix
K=AM;%Initialize matrix
Dd=AM;%Initialize matrix
AMd=AM;%Initialize matrix
AMinfd=AM;%Initialize matrix
w=zeros(1,kmax);%Initialize matrix

for l=1:1:ii
    kl=dof(l);
    for j=1:1:ii
        kj=dof(j);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NOTE: For more than 6DoF is advisable  %
        % to check the dimensionless method      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x=rem(l,6); if x==0; x=6; end
        y=rem(j,6); if y==0; y=6; end
        switch (x & y)
               case (x<=3 & y<=3), n=3;
               case (x<=3 & y>3)|(x>3 & y<=3), n=4;
               case (x>3 & y>3), n=5;
        end
        for k=1:kmax % for each frequency
            w(k)=A((k-1)*(ldof+1)+1,1);
            DimAM=1000*(A(1,2))^n;
            DimD=1000*(A(1,2))^n*w(k);
            AMinf(l,j)=(Ainf(kl,2*kj));
            AM(k,l,j)=(A((kl+1)+((k-1)*(ldof+1)),2*kj));
            D(k,l,j)=(B((kl+1)+((k-1)*(ldof+1)),2*kj));
            K(k,l,j)=(D(k,l,j)*DimD)+1i*w(k)*(AM(k,l,j)-AMinf(l,j))*DimAM;
            Dd(k,l,j)=D(k,l,j)*DimD;
            AMd(k,l,j)=AM(k,l,j)*DimAM;
            AMinfd(l,j)=AMinf(l,j)*DimAM;
        end
    end
end

%Interpolate the matrices to get a evenly spaced K and D matrices

expnt=8; %Interpolated sample number {2^expnt}                         
wi = 0.0:w(k)/(2^expnt-1):w(k); 
np=length(wi);
Ks(np,ii,ii)=0;
Ks(:,:,:)= interp1(w,K,wi,'linear','extrap');
Ds(:,:,:)= interp1(w,Dd,wi,'linear','extrap');

DKs= reshape(max(D),ii,ii);
DKs=max(diag(DKs));

%Get the significant Entries of K
for s=1:1:ii
    for r=s:1:ii 
        if r~=s
           if max(abs(D(:,s,r)))<=0.0005*DKs, %TD: Should this parameter be user defined?
              Ks(:,s,r)=0.0; 
           end
        else
           if max(abs(D(:,s,s)))<=0.0001*DKs,%TD: Should this parameter be user defined?
              Ks(:,s,s)=0.0; 
              disp([identify(dof(s)),' - Unimportant DoF']);
           end 
        end
    end
end
I=find(wi>=Rad.twr(1) & wi<=Rad.twr(2)); %most relevant freq. range
wwac=1-Rad.wwf;                          %freq. weights in the a&c intervals (marginal ranges)
wwb=(kmax/length(I)-1)*Rad.wwf+1;        %freq. weights in the b interval (most relevant range)
wiw=[wwac*ones(1,I(1)-1), wwb*ones(1,length(I)),...
     wwac*ones(1,length(wi)-length(I)-I(1)+1)];

switch Rad.tfi % Fitting method to use
       case 1, [ssRad] = FreqID(Rad,Ks,HM,wi,wiw,dof,ii); %Freq. Identification
       case 2, [ssRad] = FDI(Rad,Ks,AMd,AMinfd,Dd,wi,w,dof,ii); %FDI Toolbox
       case 3, [ssRad] = TimeID(Rad,Ks,wi,wiw,dof,ii,Ds,DKs); % TD LS Method
       case 4, [ssRad] = TimeID2(Rad,Ds,DKs,wi,dof,ii); %TD Realization Theory
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%% PARAMETRIC MODEL -Frequency domain Identification- %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ssRad] = FreqID(Rad,Ks,HM,wi,wiw,dof,ii)
fprintf('Using identification in the frequency domain.\n') %Method used

Kss=zeros(length(Ks),1);  %Initializations
KS=zeros(length(Ks),1);
ssRad(:).D=[];
N=zeros(1,1);
q=1;

for s=1:1:ii % for each entry of the matrix K(w)
    for r=s:1:ii % 
        if max(abs(Ks(:,s,r)))>0.0 %If it is significant
           Rsquare=0; m=2; n=1;
           Kss(:,q)=Ks(:,s,r);
           
           while Rsquare<Rad.fit %LS fit
                 [b,a] = invfreqs(Kss(:,q),wi,n,m,wiw);%,'trace');%,'trace'
                 Ksa=freqs(b,a,wi).';
                 res=Kss(:,q)-Ksa;                     %Residuals
                 avrg=mean(Kss(:,q).*wiw.');           %Weighted average value
                 sse=sum(wiw.'.*(res.*conj(res)));     %Weighted summed square of residuals 
                 sst=sum(wiw.'.*((Kss(:,q)-avrg).*...  %Sum of squares about the mean
                         conj(Kss(:,q)-avrg)));  
                 if sst==0, Rsquare=1-sse; else Rsquare=1-sse/sst; end 
                 m=m+1; %TD: Increase both n and m to keep the relative degree 1
                 n=n+1;
            end 
            KS(:,q)=Ksa;
            
            [A B C D]=tf2ss(b,a); %Get the equivalent ss matrices
            [~,p,~] = tf2zp(b,a);  %Zeros Poles and Gain (z,p,k)
            x=find(real(p)>0.0);
            if isempty(x)==0 
                if max(real(p(x))./abs(p(x)))<0.05
                  [T,E] = eig(A); x=find(diag(E)>0.0);
                  E(x,x)=-real(E(x,x))+1i*imag(E(x,x));
                  A=real(T*E/T);
                  % TD: Remove the free floating stability!
%                   if s==r 
%                      mode=dof(s);
%                      stability(mode,HM,A,B,C);
%                   end
               else
                   error('Radiation dynamic system UNSTABLE');
                   
               end
            else 
                  % TD: Remove the free floating stability from this
                  % toolbox
%                if s==r 
%                   mode=dof(s);
%                   stability(mode,HM,A,B,C);
%                end
            end
            %Assemble output structure
            ssRad(q).A=A; ssRad(q).B=B;
            ssRad(q).C=C; ssRad(q).D=D;
            if r==s, N(s,r)=dof(s); else N(s,r)=-1; end  %Real DoFs
            q=q+1;
        end
    end 
end
%% Fit Visualisation %%%
if Rad.ppmf==1
   Func='Transfer Function'; %Identification based on a TF 
   [i,j] = find(N.');
   mRd.idx=[i,j];            %Indexes of the radiation models 
   mRd.n=nnz(N);             %Number of radiation models 
   mRd.Ap=KS;                %Models aproximation
   mRd.Or=Kss;               %Real radiation models
   mRd.it=wi;                %Frequency interval 
   visualisation(mRd,Func);  %Fit Visualisation (FDI)
end
 %% New (reduced) convolution array %%% 
 %Comented out TD: In order to garantee that the output matrices are always
 %refering to 6 DOF's functions!
 
idx=find(diag(N)==0);
for j=1:length(idx)        % Loop to verify if there is a coss-copling % 
%     if nnz(N(:,idx(j)))>0  % term higher than the diagonal value       %
%        idx(j)=0;
%     end
end
N = removerows(N.',nonzeros(idx));
N = removerows(N.',nonzeros(idx));
[i,j,s] = find(N.'); s(s<0)=0;
for u=1:length(i)
    ssRad(u).ij=[j(u),i(u),s(u)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% PARAMETRIC MODEL -FDI Toolbox Identification- %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ssRad]=FDI(Rad,Ks,AM,AMinf,D,wi,w,dof,ii)

    fprintf('Using the FDI Toolbox v1.2, for identification in the frequency domain.\n') %Method being used
    addpath('FDI_Toolbox_v1.2'); %Add the FDI folder
    
    %Input option structure
    FDIopt.OrdMax=20; % - Maximum order of transfer function to be considered. Typical value 20.
    FDIopt.AinfFlag=1;% - if set to 1, the algorithm uses Ainf (3D hydrodynamic data), %if set to 0, the algorithm estimates Ainf (2D hydrodynamic data);
    FDIopt.Method=2; % - There are 3 parameter estimation methods (1,2,3). See help of fit_siso_fresp. Recomended method 2 (best trade off between accuracy and speed)
    FDIopt.Iterations=20;% - Related to parameter estimation methods 2 and 3. See help of fit_siso_fresp. Typical value 20.
    FDIopt.PlotFlag=0;% - If set to 1 the function plots the results of each iteration of the automatic order detection in a separate figure.
    FDIopt.LogLin=1;%  - logarithmic or linear frequency scale for plotting.
    FDIopt.wsFactor=0.1;% - Sample faster than the Hydrodynamic code for plotting. Typical value 0.1.
    FDIopt.wminFactor=0.1;% - The minimum frequency to be used in the plot is FDIopt.wminFactor*Wmin, where %Wmin is the minimum frequency of the dataset used for identification.%Typical value 0.1.
    FDIopt.wmaxFactor=5; %- the maximum frequency to be used in the plot is FDIopt.wmaxFactor*Wmax, where Wmax is the maximum frequency of the dataset used for identification. Typical value 5.
    
    %Initializations
    Kss=zeros(length(Ks),1);  
    KS=zeros(length(Ks),1);
    ssRad(:).D=[];
    N=zeros(1,1);
    q=1;
    W=wi;
    Ainf=AMinf;
    
    %Interpolation of the Added mass and damping matrices
    np=length(wi);
    AMM(np,ii,ii)=0;
    AMM(:,:,:)= interp1(w,AM,wi,'linear','extrap');
    
    BMM(np,ii,ii)=0;
    BMM(:,:,:)= interp1(w,D,wi,'linear','extrap');
    
    for s=1:1:ii % for each...
        for r=s:1:ii % entry of the matrix ..
            if max(abs(Ks(:,s,r)))>0.0
                
                Kss(:,q)=Ks(:,s,r);%Real model
                Dof=[s,r];
                
                %Call the FDI toolbox
                [KradNum,KradDen,Ainf_hat]=FDIRadMod(W',AMM(:,s,r),Ainf(s,r,1),...
                    BMM(:,s,r),FDIopt,Dof,Rad.twr(1),Rad.twr(2),Rad.fit);
                
                b=KradNum;
                a=KradDen;

                KS(:,q)=freqs(b,a,wi);
                [A B C D]=tf2ss(b,a);
                [~,p,~] = tf2zp(b,a);  %Zeros Poles and Gain (z,p,k)
                
                x=find(real(p)>0.0);
                if isempty(x)==0 
                    if max(real(p(x))./abs(p(x)))<0.05
                      [T,E] = eig(A); x=find(diag(E)>0.0);
                      E(x,x)=-real(E(x,x))+1i*imag(E(x,x));
                      A=real(T*E/T);
                      if s==r 
                         mode=dof(s);
                         stability(mode,HM,A,B,C);
                      end
                   else
                       error('Radiation dynamic system UNSTABLE');
                   end

                end
                ssRad(q).A=A; ssRad(q).B=B;
                ssRad(q).C=C; ssRad(q).D=D;
                if r==s, N(s,r)=dof(s); else N(s,r)=-1; end  %Real DoFs
                q=q+1;
            end
        end
    end
    
    close 102 %Window
    %% Fit Visualisation %%%
if Rad.ppmf==1
   Func='Transfer Function'; %Identification based on a TF 
   [i,j] = find(N.');
   mRd.idx=[i,j];            %Indexes of the radiation models 
   mRd.n=nnz(N);             %Number of radiation models 
   mRd.Ap=KS;                %Models aproximation
   mRd.Or=Kss;               %Real radiation models
   mRd.it=wi;                %Frequency interval 
   visualisation(mRd,Func);  %Fit Visualisation (FDI)
end
    
%% New (reduced) convolution array %%% 

 %Comented out TD: In order to garantee that the output matrices are always
 %refering to 6 DOF's functions!
% idx=find(diag(N)==0);
% for j=1:length(idx)        % Loop to verify if there is a coss-copling % 
%     if nnz(N(:,idx(j)))>0  % term higher than the diagonal value       %
%        idx(j)=0;
%     end
% end
% N = removerows(N.',nonzeros(idx));
% N = removerows(N.',nonzeros(idx));
[i,j,s] = find(N.'); s(s<0)=0;
for u=1:length(i)
    ssRad(u).ij=[j(u),i(u),s(u)];
end

    rmpath('FDI_Toolbox_v1.2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%% PARAMETRIC MODEL -Time domain Least Squares Identification- %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ssRad] = TimeID(Rad,Ks,wi,wiw,dof,ii,D,DKs)
fprintf('Using time domain identification with prony.m\n') %Method Used


%TD: Please uncomment this section to use the IFFT method
% Kss=zeros(length(Ks),1);   %Initializations
% Kts=zeros(length(Ks),1);
% KT=zeros(length(Ks),1);
% ssRad(:).D=[];
% N=zeros(1,1);
% q=1;
% fmax=wi(end)/(2*pi);   %Freq max [Hz] or Nyquist frequency
% dT=1/(2*fmax);         %Sampling period
% dw=wi(3)-wi(2);
% Ts=pi/(wi(end)+dw);
% wzeros = wiw'; 
% wzeros(wiw>0)=1; 

% %%% Using the B(w) integration instead of FFT!!!!!!!!!!!!!!!!!!!!!!!
% TD 23/10/12

%Initializations
Kss=zeros(length(Ks),1);   
ssRad(:).D=[];
N=zeros(1,1);
q=1; %Index for the number of significant entries of Kt
fmax=wi(end)/(2*pi);   %Freq max [Hz] or Nyquist frequency
dT=0.1;%1/(2*fmax);         %Sampling period
dw=wi(3)-wi(2); %Frequency sample
Ts=pi/(wi(end)+dw);
wzeros = wiw'; 
wzeros(wiw>0)=1; 
stime=0:dT:200;%dT*(0:length(wi)-1);%Time vector
Bw=D;%Damping matrix

% Remove the unimportant terms of the D matrix
for s=1:1:ii
    for r=s:1:ii 
        if r~=s
           if max(abs(D(:,s,r)))<=0.05*DKs,
              %if max(abs(Ks(:,s,r)))<=0.05*DKs, 
              Bw(:,s,r)=0.0; 
           end
        else
           if max(abs(D(:,s,s)))<=0.0001*DKs,
              %if max(abs(Ks(:,s,s)))<=0.0001*DKs, 
              Bw(:,s,s)=0.0; 
%               disp([identify(dof(s)),' - Unimportant DoF']);
           end 
        end
    end
end


for s=1:1:ii %for each entry of the matrix
    for r=s:1:ii
        if max(abs(Ks(:,s,r)))>0.0 
            % %%% Using the B(w) integration instead of FFT!!!!!!!!!!!!!!!!!!!!!!!
            % TD 23/10/12
                     
           %Impluse response function integration
           for t=1:length(stime) %for each time step
               COS=0;
               for i=2:length(wi)-1 %for each frequency
                    COS=COS+2*Bw(i,s,r).*cos(wi(i)*stime(t));
               end               
               Ktemp(t,1)=dw/pi*(COS+Bw(1,s,r)+Bw(end,s,r)*cos(wi(end)*stime(t)));               
           end
           
           Kts(:,q)=Ktemp;%Store the time decay function
            
           % Please uncomment this to use the IFFT Method 
           
           %kt0 = (ifft([0;Ks(2:end,s,r);0;flipud(conj(Ks(2:end,s,r)))],'symmetric')); 
%            kstemp = Ks(:,s,r);%.*wzeros;
%            Ktemp =(ifft([0;kstemp(2:end);0;flipud(conj(kstemp(2:end)))],'symmetric'));
%            np = round(length(Ktemp)/2);
            np=length(Ktemp);
%            Kts(:,q) = Ktemp(1:np); 
%            stime=dT*(0:length(Kts(:,q))-1).';
%            Kss(:,q) = Ks(:,s,r);

            %Initializations
           Rsquare=0; m=-1; n=0;

           while (Rsquare<Rad.fit)
                 m=m+1; %Increase both to have a 1 relative degree
                 n=n+1;
                 
                 [b,a]=prony(Kts(:,q),m,n);
                 
                 Kta=impz(b,a,np);
                 sse=(Ktemp(1:np)-Kta)'*(Ktemp(1:np)-Kta); %Weighted summed square of residuals 
                 sst=((Ktemp(1:np)-mean(Ktemp(1:np)))'*(Ktemp(1:np)-mean(Ktemp(1:np)))); %Sum of squares about the mean
                 if sst==0, 
                    Rsquare=1-sse;
                 else
                    Rsquare=1-sse/sst;
                 end
           end
%            sys=tf(b,a,dT);
%            [A,B,C,D,Ts]=ssdata(d2c(sys));
%            
           
           
           [Ad Bd Cd Dd]=tf2ss(b,a);                  %Discrete domain A,B,C,D state-space arrays
           sys=(d2c(ss(Ad,Bd,Cd*dT,0,dT),'tustin')); %Continuous domain A,B,C,D state-space arrays
           [A B C D]=ssdata(sys);
           
           Rsquare=0;
           n=1;
%            while Rsquare<Rad.fit*0.9
%                     
%                     n=n+1;%Increase the order                  
%                 
%                     [AM,BM,CM,DM,TOTBND,HSV] =balmr(A,B,C,0,1,n);
%                     [Y(:,1),T]=impulse(ss(AM,BM,CM,DM,0),stime);
%                     
%                     res=Kts(:,q)-Y;                     %Residuals
%                     avrg=mean(Kts(:,q));             %Sum of squares about the mean
%                     sst=sum(((Kts(:,q)-avrg).*conj(Kts(:,q)-avrg)));  %Weighted average value
% %                     sse=sum(wiw.'.*(res.*conj(res)));     %Weighted summed square of residuals 
%                     Rsquare=1-sum(res.^2)/sst; %Sum of squares about the mean
%                     
% %                     if 
% % %                  if sse==0, Rsquare=1-sse; else Rsquare=1-sse/sst; end  
% %                  if m==n, m=m+1; elseif n<m, n=n+1; end  %n<=m for LTI systems 
%                     if n==length(A)
%                         Rsquare=1;
%                         disp('No state was removed')
%                     end
%            end 
% 
%            Rsquare
           
           KT(:,q)=impulse(ss(A,B,C,0),stime);
%            hankelsv(ss(A,B,C,0));
%            title(['K_{',num2str(s),num2str(r),'}'])
           
           %Stability of the system
           [~,p,~] = zpkdata(ss(A,B,C,0));
           unstable=0;
           for j=1:length(p)
               if real(p{j})>0 
                  unstable=1; 
               end 
           end       
           if unstable==1, error('Radiation dynamic system UNSTABLE'); end
           
           %Assemble the output
          ssRad(q).A=A; ssRad(q).B=B;
          ssRad(q).C=C; ssRad(q).D=D;
          if r==s, N(s,r)=dof(s); else N(s,r)=-1; end  %Real DoFs          
          q=q+1;
        end 
    end 
end
%% Fit Visualisation %%%
if Rad.ppmf==1
   Func='Impulse Response Function'; %Identification besed on a IRF 
   [i,j] = find(N.');
   mRd.idx=[i,j];                 %Indexes of the radiation models 
   mRd.n=nnz(N);                  %Number of radiation models 
   mRd.Ap=KT;                     %Models aproximation
   mRd.Or=Kts;                    %Real radiation models
   mRd.it=stime;                  %Simulation time interval 
   visualisation(mRd,Func);       %Fit Visualisation (TDI)
end
%% New (reduced) convolution array %%% 
 %Comented out TD: In order to garantee that the output matrices are always
 %refering to 6 DOF's functions!
idx=find(diag(N)==0);
for j=1:length(idx)        % Loop to verify if there is a coss-copling % 
%     if nnz(N(:,idx(j)))>0 % term higher than the diagonal value       %
%        idx(j)=0;
%     end
end
N = removerows(N.',nonzeros(idx));
N = removerows(N.',nonzeros(idx));
[i,j,s] = find(N.'); s(s<0)=0;
for u=1:length(i)
    ssRad(u).ij=[j(u),i(u),s(u)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% PARAMETRIC MODEL -Time domain - Realization Theory!!!- %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ssRad] = TimeID2(Rad,D,DKs,wi,dof,ii)
fprintf('Using time domain identification, with K(t) and imp2ss.m.\n')%Method chosen

%Initializations
ssRad(:).D=[];
N=zeros(1,1);
q=1; %Index for the number of significant entries of Kt
dT=0.1;%1/(2*fmax);         %Sampling period
dw=wi(3)-wi(2); %Frequency sample
stime=0:dT:200;%dT*(0:length(wi)-1);%Time vector
Bw=D;%Damping matrix

% Remove the unimportant terms of the D matrix
for s=1:1:ii
    for r=s:1:ii 
        if r~=s
           if max(abs(D(:,s,r)))<=0.05*DKs,
             Bw(:,s,r)=0.0; 
           end
        else
           if max(abs(D(:,s,s)))<=0.0001*DKs,
              Bw(:,s,s)=0.0; 
           end 
        end
    end
end


for s=1:1:ii %for each entry of the matrix
    for r=s:1:ii
        if max(abs(Bw(:,s,r)))>0.0 
           
           %Impluse response function integration
           for t=1:length(stime) %for each time step
               COS=0;
               for i=2:length(wi)-1 %for each frequency
                    COS=COS+2*Bw(i,s,r).*cos(wi(i)*stime(t));
               end
               Ktemp(t,1)=dw/pi*(COS+Bw(1,s,r)+Bw(end,s,r)*cos(wi(end)*stime(t)));             
           end
           
           Kt(:,q)=Ktemp;%Store the time decay function
           
           %Fitting
           [A,B,C,D,totbnd,hsv] = imp2ss(Ktemp,dT,1,1,Rad.fit);
           [y,t]=impulse(ss(A,B,C*dT,D*dT,0),stime(end));
           
           %Reduction of the number of states
           a=0;
           if Rad.fmt==1 %If Manual order reduction is selected
                while a==0 %Is the user satisfied with the goodness of the fit
               
                       [AM,BM,CM,DM,TOTBND,HSV] =balmr(A,B,C*dT,0,3); %System must be casual (D=0)
                       [Y,T]=impulse(ss(AM,BM,CM,DM,0),stime(end));  

                       figure('NumberTitle','off','Name','Impulse Response Function')
                       pp=[3.0 2.0 4 3.5]; set(gcf,'paperposition',pp,'position',pp*100)
                       plot(stime,Ktemp,t,y,T,Y)
                       legend('K(t)','~K(t)-original (i=200)',strcat('~K(t)-reduced (i=',num2str(length(AM)),')'))
                       xlabel('Time [s]')
                       if s==r
                           if s>3
                               ylabel('[kg.m^2/s^2]')
                           else
                               ylabel('[kg/s^2]')
                           end
                       else
                           ylabel('[kg.m/s^2]')
                       end
                       title(strcat('K_{',num2str(s),num2str(r),'}(t)'))
                       
                   a=input('Is the state reduction good (1), or do you want to re-do it (0)?');    
                end
           else % Automatic order detection
               n=1;
               Rsquare=0;
               while Rsquare<Rad.fit
                    
                    n=n+1;%Increase the order
                
                    [AM,BM,CM,DM,TOTBND,HSV] =balmr(A,B,C*dT,0,1,n);%System must be casual (D=0)
                    [Y(:,1),T]=impulse(ss(AM,BM,CM,DM,0),stime);
                    
                    res=Kt(:,q)-Y;                     %Residuals
                    avrg=mean(Kt(:,q));             %Sum of squares about the mean
                    sst=sum(((Kt(:,q)-avrg).*conj(Kt(:,q)-avrg)));  %Weighted average value
                    Rsquare=1-sum(res.^2)/sst %Sum of squares about the mean
               end 
                
             
               if Rad.ppmf==1 %Figure plot
                  figure('NumberTitle','off','Name','Impulse Response Function')
                  pp=[3.0 2.0 4 3.5]; set(gcf,'paperposition',pp,'position',pp*100)
                       plot(stime,Ktemp,t,y,T,Y)
                       legend('K(t)','~K(t)-original (i=200)',strcat('~K(t)-reduced (i=',num2str(length(AM)),')'))
                       xlabel('Time [s]')
                       if s==r
                           if s>3
                               ylabel('[kg.m^2/s^2]')
                           else
                               ylabel('[kg/s^2]')
                           end
                       else
                           ylabel('[kg.m/s^2]')
                       end
                       title(strcat('K_{',num2str(s),num2str(r),'}(t)')) 
               end
           end
 
           %Stability of the system
           [~,p,~] = zpkdata(ss(AM,BM,CM,DM));
           unstable=0;
           for j=1:length(p)
               if real(p{j})>0 
                  unstable=1; 
               end 
           end       
           if unstable==1, error('Radiation dynamic system UNSTABLE'); end
          
           %Output assembly
          ssRad(q).A=AM; ssRad(q).B=BM;
          ssRad(q).C=CM; ssRad(q).D=DM;
          if r==s, N(s,r)=dof(s); else N(s,r)=-1; end  %Real DoFs          
          q=q+1;
        end 
    end
    
end

%% New (reduced) convolution array %%%
 %Comented out TD: In order to garantee that the output matrices are always
 %refering to 6 DOF's functions!
idx=find(diag(N)==0);
for j=1:length(idx)        % Loop to verify if there is a coss-copling % 
%     if nnz(N(:,idx(j)))>0 % term higher than the diagonal value       %
%        idx(j)=0;
%     end
end
N = removerows(N.',nonzeros(idx));
N = removerows(N.',nonzeros(idx));
[i,j,s] = find(N.'); s(s<0)=0;
for u=1:length(i)
    ssRad(u).ij=[j(u),i(u),s(u)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%% Function to Identify each mode %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function md=identify(d)
switch (d)
       case 1, md='Surge';
       case 2, md='Sway';
       case 3, md='Heave';
       case 4, md='Roll';
       case 5, md='Pitch';
       case 6, md='Yaw';
       otherwise, md=([num2str(d),'DoF']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%% Function to display the Parametric Model Fit (TDI or FDI) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
function visualisation(mRd,Func) 
for fig=1:mRd.n
    etr=[num2str(mRd.idx(fig,1)),'_',num2str(mRd.idx(fig,2))];
    figure('NumberTitle','off','Name',Func)
    pp=[3.0 2.0 4 3.5]; set(gcf,'paperposition',pp,'position',pp*100)
    plot(mRd.it,real(mRd.Ap(:,fig)),'g','linewidth',2), hold on
    plot(mRd.it,real(mRd.Or(:,fig)),'b','linewidth',2)
    if length(Func)==17 
       title(['K_',etr],'fontsize',10,'FontWeight','bold')
       xlabel('\omega [rad/s]'); ylabel('K(s)'); 
       %legend('\simK(s)','K(s)',1)
    elseif length(Func)==25    
       title(['K_',etr],'fontsize',12,'FontWeight','bold')
       xlabel('Time [s]'); ylabel('K(t)[kN]'); 
       legend('\simK(t)','K(t)',1)
    end
    set(get(gca,'YLabel'),'fontsize',9,'fontweight','b')
    set(get(gca,'XLabel'),'fontsize',9,'fontweight','b')
    grid on
end