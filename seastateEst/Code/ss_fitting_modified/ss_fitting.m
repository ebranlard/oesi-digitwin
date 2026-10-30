function ss_fitting(Rad_in)
%% State Space fitting of the convolution integral
% This routine fits a state-space model for the convolution term of wave
% radiation forces on a floating body, based on the coefficients of WAMIT. 
% All the inputs should be define in the "opt_file". 
% To run the code please type in the comand window: 
% >> ss_fitting('SS_Fitting_Options.inp')
% Please refer to the Theory and UserManual for more information.
%
% This routine was developed by:
% Tiago Duarte
% tduarte@hidro1.ist.utl.pt
% Instituto Superio Tecnico - IST Lisbon, Portugal
%
% in collaboration with:
% WavEC Offshore Renewables (www.wavec.org)
% NREL (www.wind.nrel.gov)
% V1.00.01

%% Import options from python

% Assemble structures
gp = struct();
gp.DoF = double(Rad_in.DoF); 

ULEN = 1; % Reference length of the WAMIT files

try
    Rad = struct();
    Rad.FileName = char(Rad_in.FileName);
    Rad.twr = double(Rad_in.twr);
    Rad.wwf = double(Rad_in.wwf);
    Rad.tfi = double(Rad_in.tfi);
    Rad.fit = double(Rad_in.fit);
    Rad.ppmf = double(Rad_in.ppmf);
    Rad.fmt = double(Rad_in.fmt);
catch
    disp('Error parsing the input variables passed from Python.')
    return
end
if length(Rad.twr)~=2
    disp('Frequency range vector must contain minimum and maximum value.')
    return
elseif Rad.twr(1)<0
    disp('Minimum frequency of the frequency range should be equal or larger than zero.')
    return
elseif Rad.twr(2)<Rad.twr(1)
    disp('The frequency range should be defined as:')
    disp('[Fmin Fmax]')
    disp('Where:')
    disp('Fmin<Fmax and Fmin>=0')
    return
end

if isempty(Rad.wwf)==1
    disp('Frequencies weighting factores must be between 0 and 1')
    return 
elseif (Rad.wwf<0 || Rad.wwf>1) 
    disp('Frequencies weighting factores must be between 0 and 1')
    return 
end

if (isempty(Rad.tfi)|| Rad.tfi<1 || Rad.tfi>5 )
    
    disp('The method chosen need to be 1 to 4:')
    disp('    1 - FREQ Method')
    disp('    2 - FDI Toolbox')
    disp('    3 - LS Time-Domain Method')
    disp('    4 - Realization Theory')
    return
end

if (isempty(Rad.fit)||Rad.fit<0 || Rad.fit>=1 )
    
    disp('The goodness of the fit must be between 0 and 1')
    return
end

if (isempty(Rad.ppmf) || (Rad.ppmf~=0 && Rad.ppmf~=1) )
    
    disp('The plot flag must be 0 or 1')
    return
end

if (isempty(Rad.fmt) || (Rad.fmt~=0 && Rad.fmt~=1) )
    
    disp('The order reduction flag must be 0 or 1')
    return
end
%% Reads WAMIT .1 file
% Get the AddedMass, INFAddedMass and Damping Matrices
A=[];
B=[];
   fid1 = fopen(strcat(Rad.FileName,'.1'));
   if fid1==-1
        disp(['The WAMIT file "' Rad.FileName '" was not found.'])
        disp(' The file name can contain relative or absolute path, and should not contain the extension ".1".')
        disp(' Ex.: "HydroData/Spar" or "C:/HydroData/Spar".' )
        return
   end
   OPT1 = textscan(fid1,'%f %n %n %f %f','commentstyle','W');
   fclose(fid1); 
   ctr=OPT1{1}(1); %-1 if 0 freq, 0 if Inf freq
   if ctr~=0, 
       n=1; 
       while OPT1{1}(n)==OPT1{1}(n+1), 
           n=n+1; 
       end
   end
   if ctr==0, %Inf freq
       n=length(OPT1{1}); 
   end
   m=max(OPT1{2}); frm='%3.0f\t%2.6E\t ';
   AM(m,2*m)=0; DM(m,2*m)=0;   
   for s=2:m, frm=horzcat(frm,'%3.0f\t%2.6E\t'); end
   frm=horzcat(frm,'\n');
   for ij=1:1:m, AM(:,2*ij-1)=(10+ij:10:10*m+ij); end
   if ctr<0,
       OPT1{5}(1:n)=0.0; 
   end
   if ctr==0, fid1a=fopen('INFAddedMass.log','a+'); else

       for ij=1:1:m, DP(:,2*ij-1)=(10+ij:10:10*m+ij); end
   end
   for r=0:n:length(OPT1{1})-n
       for k=1+r:1:n+r  
           b=horzcat(OPT1{2}(k),OPT1{3}(k));
           AM(b(2),2*b(1))=OPT1{4}(k);
           AM(b(1),2*b(2))=OPT1{4}(k);
           if ctr~=0.0
              DP(b(2),2*b(1))=OPT1{5}(k);
              DP(b(1),2*b(2))=OPT1{5}(k);
           end
           d=mod(k,m);
           if d==0
              switch ctr & r
                  %ver melhor a freq inf & zero
                  case (-1 & 0), FREQ=0.0; frms=' %3.3f\t%3.3f\n';
                  case (0 & 0), FREQ='INF FREQ'; frms=' %8s\t%3.3f\n';
                  otherwise, FREQ=2*pi/OPT1{1}(k); frms=' %3.3f\t%3.3f\n';
              end
           end
       end
   if isfinite(FREQ)
        A=[A; FREQ ULEN zeros(1,11); AM zeros(6,1)];
        B=[B; FREQ ULEN zeros(1,11); DP zeros(6,1)];
   else
        Ainf=AM;
   end
   end    

ldof=length(Ainf(:,1)); %Nº of DoF in the wamit output array (gener. modes)
ii=length(gp.DoF)-ldof; %Assigned modes vector need to have length = ldof
if sum(gp.DoF)==0
    disp('Error in reading the DoF vector.')
    disp('At least one DoF must be active (non-zero).')
    return
elseif ii~=0 
    disp('Error in reading the DoF vector.')
    disp('The length of the vector must match the number of Degrees of Freedom of the WAMIT file.')
    return
end


for d=1:1:ldof
    if (gp.DoF(d)~=0 & gp.DoF(d)~=1), n='stop'; return, else n='go'; end
end
if strcmp('stop', n)==1, disp('DoF Input Error'), return, end 

%% Calls the Radiation function to fit the State-Space Model

HM=[];

ssRad = Radiation(gp,Rad,A,Ainf,B,HM);  

gp.dof=nonzeros(vertcat(ssRad.ij)*[0 0 1].');

%% Assembly of the arrays (A,B,C) of the state-space representation
spdof=zeros(1,6); %Vector containing states per dof
ii=6;%length(gp.dof); TD: Keep the B and C matrices with NX6 and 6xN
for k=1:length(ssRad)
    
    a=ssRad(k).ij;
    if a(1)~=a(2) & a(1)<a(2)% Its a cross term!
        %Find the other cross term
        kk=0;
        kkk=k;
        while kk==0
            if ssRad(kkk).ij(1)==a(2)
                kk=1;
                kkk=kkk+1;
            else
                kkk=kkk+1;
            end
        end
        % Add a new entry
        ssRad(kkk+1:end+1)=ssRad(kkk:end);
        ssRad(kkk)=ssRad(k);
        ssRad(kkk).ij=[a(2) a(1)];
    end
end

for k=1:1:length(ssRad)
    a=ssRad(k).ij;
    if k==1, %first mode
        AA=ssRad(1).A;
        spdof(1,a(1))=spdof(1,a(1))+length(ssRad(1).A);
    else %other modes       
       AA=blkdiag(AA,ssRad(k).A);
       spdof(1,a(1))=spdof(1,a(1))+length(ssRad(k).A);
    end
    
end
CC(ii,length(AA))=0;       
BB(length(AA),ii)=0;  

for k=1:1:length(ssRad)
    a=ssRad(k).ij;
    lc(k)=length(ssRad(k).C);
    lb(k)=length(ssRad(k).B);
    if k==1
       CC(1,1:lc(1))=-ssRad(k).C; %C array of the rad. state-space repr.
       c(1)=lc(1);
       BB(1:lb(1),1)=ssRad(k).B;  %B array of the radiation state-space repr.
       e(1)=lb(1);
    else 
       e(k)=e(k-1)+lb(k); 
       BB(e(k-1)+1:e(k),a(1))=ssRad(k).B; 
       c(k)=c(k-1)+lc(k); 
       CC(a(2),c(k-1)+1:c(k))=-ssRad(k).C;
    end
end

%% Check stability and other proprieties of the system

sys=ss(AA,BB,CC,0,0);

if isstable(sys)~=1 %Stability
   error('Resulting SS system unstable! Try with a lower R^2 value or different method!')
end

if isproper(sys)~=1 %Proper
    warning('The system is not proper!')
end

%General calculations
[mag,phase,wout] = bode(sys,0:0.1:5); %Frequency response (tf)
jj=0;
for ii=1:6
    jj=jj+sum(mag(ii,ii,:)<0);
end

if jj>0 %Check for passivity
    warning('The system is not passive! Choose other R^2 value or different method!')
end

if sum(sum(mag(:,:,1)))~=0 %Low Frequency limit
    warning('Low frequency limit not fulfilled! Only guarantee in the FDI method!')
end

%% Save the Matrix AA BB and CC in a new ss file

formata='%6.6e '; %Format of AA
for ii=1:size(AA,2)-1
    formata=[formata '%6.6e '];
end
formata=[formata '\r\n']; 

formatb='%6.6e ';%Format of BB
for ii=1:size(BB,2)-1
    formatb=[formatb '%6.6e '];
end
formatb=[formatb '\r\n'];

formatc='%6.6e ';%Format of CC
for ii=1:size(CC,2)-1
    formatc=[formatc '%6.6e '];
end
formatc=[formatc '\r\n'];

% Writes output file
fid=fopen(strcat(Rad.FileName, '.ss'),'w+');%'Verification_tests/marin_semi_097/marin_semi_097_Time.ss','w+');

fprintf(fid,'%s\r\n',['SS_Fitting v1.00.01: State-Spaces Matrices obtained using ', identify_method(Rad.tfi),' method, on ',datestr(now)]);
fprintf(fid,'%s\r\n',[num2str(gp.DoF) '    %Enabled DoFs']);
fprintf(fid,'%s\r\n',[num2str(size(AA,2)) '                  %Radiation states']);
fprintf(fid,'%s\r\n',[num2str(spdof) '    %Radiation states per DOFs']);

fprintf(fid,formata,AA');

fgetl(fid);%skip one line

fgetl(fid);%skip one line

fprintf(fid,formatb,BB');

fgetl(fid);%skip one line

fgetl(fid);%skip one line
fprintf(fid,formatc,CC');

fclose(fid);

% End of program
disp(['Fitted a SS model with ',num2str(length(AA)),' states.'])
disp(['Results saved in ', Rad.FileName,'.ss'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%% Function to Identify each method %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function md=identify_method(d)
switch (d)
       case 1, md='Freq_ID';
       case 2, md='FDI Toolbox';
       case 3, md='TD-LS';
       case 4, md='TD-Realization Theory';
end
