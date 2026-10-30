function [KradNum,KradDen,Ainf_hat]=FDIRadMod(W,A,Ainf,B,FDIopt,Dof,wmin,wmax,r2Thres)
%Frequency-domain identification of radiation function models of marine
%structures. This function identifies a SISO transfer function
%corresponding to the coupling specified.
%
%Use: [KradNum,KradDen,Ainf_hat]=FDIRadMod(W,A,Ainf,B,FDIopt,Dof)
%
%W - is the vector of frequencies at which A(w) and B(w) are compued.
%
%A - is the vector of frequency dependant added mass coeffiricents A(w),
%
%Ainf -  is the infinite frequency added mass coefficient,
%
%B - is the vector of potential damping coefficients B(w),
%
%FDIopt - is a structure with the following fields,
% 
%FDIopt.OrdMax - Maximum order of transfer function to be considered.
%Typical value 20.
%
%FDIopt.AinfFlag - if set to 1, the algorithm uses Ainf (3D hydrodynamic data), 
%if set to 0, the algorithm estimates Ainf (2D hydrodynamic data);
%
%FDIopt.Method - There are 3 parameter estimation methods (1,2,3). See help
%of fit_siso_fresp. Recomended method 2 (best trade off between accuracy and speed)
%
%FDIopt.Iterations - Related to parameter estimation methods 2 and 3. See help of
%fit_siso_fresp. Typical value 20.
%     
%FDIopt.PlotFlag - If set to 1 the function plots the results of each iteration
%of the automatic order detection in a separate figure.
%
%FDIopt.LogLin  - logarithmic or linear frequency scale for plotting.
%
%FDIopt.wsFactor - Sample faster than the Hydrodynamic code for plotting. Typical value
%0.1.
%
%FDIopt.wminFactor - The minimum frequency to be used in the plot is FDIopt.wminFactor*Wmin, where
%Wmin is the minimum frequency of the dataset used for identification.
%Typical value 0.1.
%
%FDIopt.wmaxFactor - the maximum frequency to be used in the plot is FDIopt.wmaxFactor*Wmax, where
%Wmax is the maximum frequency of the dataset used for identification.
%Typical value 5.
%
%Dof [i j] - vector with the coupling to be indentified. i=1,2,..,6. j=1,2,...,6. 
%
%
%Description:
%
%The function first calls EditAB.m to prepare the data for identification. 
%Then, depending on the option FDIopt.AinfFlag, the function calls the 
%appropriate computation routine. 
%
%The function also makes an automatic order estimate by 
%increasing the order of the approximation and computing the coefficient of 
%determination related to the fitting of both added mass and damping. When 
%both these coefficients reach a value greater or equal to 0.99, the 
%function stops increasing the order, and the re-constructued added mass 
%and damping are plotted together with the non-parametric data used for 
%identification. At this point, the function prompts the user to either 
%adjust the order of the approximation manually via a keyboard input or  
%leave the model as it is and exit the function. The user can make as many 
%changes in order as required, and every time there is a change in the 
%order, the model is re-estimted and the data re-plotted.  
%
% Created by Tristan Perez (trisan.perez@ntnu.no)
% Date 2009/9/1, Trondheim, Norway.
% Revision:
%
% Copyright (C) 2009 Thor I. Fossen and Tristan Perez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% E-mail: contact@marinecontrol.org
% URL:    <http://www.marinecontrol.org>

%% Prepapre data for identification: Select frequency range and eliminate
[A,B,W]=EditAB(A,B,W,wmin,wmax);

%% Compute the retardation function Freq response K(jw)
Kw = B+complex(0,W).*(A-Ainf*ones(size(A)));

%% Frequency domain identification
MethOpt =[FDIopt.Method;FDIopt.Iterations]; 
PlotOpt =[FDIopt.PlotFlag;FDIopt.LogLin;FDIopt.wsFactor; ...
            FDIopt.wminFactor;FDIopt.wmaxFactor];

%% Use Ainf as part of the dataset        
if FDIopt.AinfFlag==1,
    %Auto order estimates
     FDIopt.Ord = 2; %Initial order of the approximation
     [KsNum,KsDen]=ident_retardation_FD(W,Kw,FDIopt.Ord,MethOpt,PlotOpt);
     %
     Kw_hatFD=freqs(KsNum,KsDen,W);
     Brecfd = real(Kw_hatFD);
     Arecfd = imag(Kw_hatFD)./W+Ainf*ones(size(W));
     
     SSEB = (B-Brecfd)'*(B-Brecfd);
     SSTB =(B-mean(B)*ones(size(B)))'*(B-mean(B)*ones(size(B)));
     R2B = 1 - SSEB/SSTB;
     
     SSEA = (A-Arecfd)'*(A-Arecfd);
     SSTA =(A-mean(A)*ones(size(A)))'*(A-mean(A)*ones(size(A)));
     R2A = 1 - SSEA/SSTA;
     
%      r2Thres = 0.97;
     while  (R2B<r2Thres) & (R2A<r2Thres),
     FDIopt.Ord = FDIopt.Ord +1;
     if FDIopt.Ord > FDIopt.OrdMax,
             break
     end
     [KsNum,KsDen]=ident_retardation_FD(W,Kw,FDIopt.Ord,MethOpt,PlotOpt);
     %Compute coeff of determination
     Kw_hatFD=freqs(KsNum,KsDen,W);
     Brecfd = real(Kw_hatFD);
     Arecfd = imag(Kw_hatFD)./W+Ainf*ones(size(W));
     
     SSEB = (B-Brecfd)'*(B-Brecfd);
     SSTB =(B-mean(B)*ones(size(B)))'*(B-mean(B)*ones(size(B)));
     R2B = 1 - SSEB/SSTB;
     
     SSEA = (A-Arecfd)'*(A-Arecfd);
     SSTA =(A-mean(A)*ones(size(A)))'*(A-mean(A)*ones(size(A)));
     R2A = 1 - SSEA/SSTA;
     end
    
    %Manual order selection
    for l=1:1000,
    %Call radiation identification function
    [KsNum,KsDen]=ident_retardation_FD(W,Kw,FDIopt.Ord,MethOpt,PlotOpt); 
    %Compute retardation frequency response
     ws=FDIopt.wsFactor*(W(2)-W(1));
     wmax=W(length(W));
     wmin=W(1);
     Wext2 =[FDIopt.wminFactor*wmin:ws:FDIopt.wmaxFactor*wmax]'; 
     Kw_hatFD=freqs(KsNum,KsDen,Wext2);
    
    %plot retardation Freq. responses 
    Xpl=[0.01 10];
    figure(102) 
    subplot(221)
    semilogx(W,20*log10(abs(Kw)),'o-r',Wext2,20*log10(abs(Kw_hatFD)),'--b','LineWidth',2)
    xlim(Xpl)
    grid on
    legend('K(jw)',['K_{hat}(jw) order ',num2str(FDIopt.Ord)])
    ylabel('|K(jw)| [dB]') %N=Kg m/s^2 = (Kg/s) m/s  
    xlabel('Frequency [rad/s]')
    title(['Convolution Model DoF ',num2str(Dof(1)),num2str(Dof(2))])
    subplot(223)
    semilogx(W,(180/pi)*angle(Kw),'o-r',Wext2,(180/pi)*angle(Kw_hatFD),'--b','LineWidth',2)
    xlim(Xpl)
    grid on
    ylabel('Phase K(jw) [deg]') %N=Kg m/s^2 = (Kg/s) m/s  
    xlabel('Frequency [rad/s]')

    %% Reconstruction of added mass and damping
    %Compute coeff of determination
     Kw_hatFDa=freqs(KsNum,KsDen,W);
     Brecfd = real(Kw_hatFDa);
     Arecfd = imag(Kw_hatFDa)./W+Ainf*ones(size(W));
     
     SSEB = (B-Brecfd)'*(B-Brecfd);
     SSTB =(B-mean(B)*ones(size(B)))'*(B-mean(B)*ones(size(B)));
     R2B = 1 - SSEB/SSTB;
     
     SSEA = (A-Arecfd)'*(A-Arecfd);
     SSTA =(A-mean(A)*ones(size(A)))'*(A-mean(A)*ones(size(A)));
     R2A = 1 - SSEA/SSTA;
    
    
     %Construct B and A t plot
     Brecfd = real(Kw_hatFD);
     Arecfd = imag(Kw_hatFD)./Wext2+Ainf*ones(size(Wext2));
     
    figure(102) 
    subplot(224)
    semilogx(W,B,'o-r',Wext2,Brecfd,'-b','LineWidth',2)
    xlim(Xpl)
    grid on
    legend('B',['Best FD ident, order ', num2str(FDIopt.Ord)])
    ylabel('Damping') %N=Kg m/s^2 = (Kg/s) m/s  
    xlabel('Frequency [rad/s]')
    title(['Potential Damping DoF ',num2str(Dof(1)),num2str(Dof(2))])
    %
    subplot(222)
    semilogx(W,A,'o-r',Wext2,Arecfd,[Wext2(1) Wext2(length(Wext2))],[Ainf Ainf],'--b','LineWidth',2)
    xlim(Xpl)
    legend('A',['Aest FD indet, order ', num2str(FDIopt.Ord)],'Ainf')
    ylabel('Added Mass')
    xlabel('Frequency [rad/s]')
    title(['Added Mass DoF ',num2str(Dof(1)),num2str(Dof(2))])
    grid on

    %% Ask for a new iteration
%     display('---------------------------------------------------------------------')
  %  KK=input(['Current Order ',num2str((FDIopt.Ord)),', (Change or 0-to exit):']);
KK=0;
    display(' ')
    if KK==0,
        break
    end
    FDIopt.Ord=KK;
    end
end

Ainf_hat = Ainf;

%% Do not use Ainf as part of the dataset
if FDIopt.AinfFlag==0,
    % Compute Complex Coefficients 
    Ac = (B./complex(0,W)+A);

    % Frequency domain identification
    MethOpt =[FDIopt.Method;FDIopt.Iterations]; 
    PlotOpt =[FDIopt.PlotFlag;FDIopt.LogLin;FDIopt.wsFactor; ...
                FDIopt.wminFactor;FDIopt.wmaxFactor];

        
    %Auto order estimates
    FDIopt.Ord = 2;
    [AsNum,AsDen]=ident_retardation_FDna(W,Ac,FDIopt.Ord,MethOpt,PlotOpt);
    %
    Ainf_hat=AsNum(1);
    P = [AsNum(:)-AsDen(:)*Ainf_hat ;0];
    Q = AsDen(:);
    
    KsNum =P(3:FDIopt.Ord+2,1)'; %Estimated fluid memory model TF.
    KsDen =Q';                   %Estimated fluid memory model TF.
    %
    %
     Kw_hatFD=freqs(KsNum,KsDen,W);
     Brecfd = real(Kw_hatFD);
     Arecfd = imag(Kw_hatFD)./W+Ainf*ones(size(W));
     
     SSEB = (B-Brecfd)'*(B-Brecfd);
     SSTB =(B-mean(B)*ones(size(B)))'*(B-mean(B)*ones(size(B)));
     R2B = 1 - SSEB/SSTB;
     
     SSEA = (A-Arecfd)'*(A-Arecfd);
     SSTA =(A-mean(A)*ones(size(A)))'*(A-mean(A)*ones(size(A)));
     R2A = 1 - SSEA/SSTA;
     
     r2Thres = 0.99;
     
     while  (R2B<r2Thres) & (R2A<r2Thres),
         FDIopt.Ord = FDIopt.Ord +1;
         if FDIopt.Ord > FDIopt.OrdMax,
             break
         end
         [AsNum,AsDen]=ident_retardation_FDna(W,Ac,FDIopt.Ord,MethOpt,PlotOpt);
         %
         Ainf_hat=AsNum(1);
         P = [AsNum(:)-AsDen(:)*Ainf_hat ;0];
         Q = AsDen(:);
         KsNum =P(3:FDIopt.Ord+2,1)'; 
         KsDen =Q'; 
         
         %Compute coeff of determination
         Kw_hatFD=freqs(KsNum,KsDen,W);
         Brecfd = real(Kw_hatFD);
         Arecfd = imag(Kw_hatFD)./W+Ainf*ones(size(W));

         SSEB = (B-Brecfd)'*(B-Brecfd);
         SSTB =(B-mean(B)*ones(size(B)))'*(B-mean(B)*ones(size(B)));
         R2B = 1 - SSEB/SSTB;

         SSEA = (A-Arecfd)'*(A-Arecfd);
         SSTA =(A-mean(A)*ones(size(A)))'*(A-mean(A)*ones(size(A)));
         R2A = 1 - SSEA/SSTA;
     end
    
    for l=1:1000,

        [AsNum,AsDen]=ident_retardation_FDna(W,Ac,FDIopt.Ord,MethOpt,PlotOpt);
        %% Compute retardation frequency response
         ws=FDIopt.wsFactor*(W(2)-W(1));
         wmax=W(length(W));
         wmin=W(1);
         Wext2 =[FDIopt.wminFactor*wmin:ws:FDIopt.wmaxFactor*wmax]';
         Ac_hatFD=freqs(AsNum,AsDen,Wext2);
         
        %% Estimated added mass and fluid memory model from As  
        %As is the identified parametric model of Ac
        Ainf_hat=AsNum(1);
        P = [AsNum(:)-AsDen(:)*Ainf_hat ;0];
        Q = AsDen(:);
        KsNum =P(3:FDIopt.Ord+2,1)'; %Estimated fluid memory model TF.
        KsDen =Q'; %Estimated fluid memory model TF.

        
        %Compute retardation frequency response
         ws=FDIopt.wsFactor*(W(2)-W(1));
             wmax=W(length(W));
             wmin=W(1);
             Wext2 =[FDIopt.wminFactor*wmin:ws:FDIopt.wmaxFactor*wmax]';
         Kw_hatFD=freqs(KsNum,KsDen,Wext2);
         
        %plot fluid memory model frequency response 
        Xpl=[min(Wext2) max(Wext2)];

        figure(102) 
        subplot(211)
        semilogx(Wext2,20*log10(abs(Kw_hatFD)),'--b','LineWidth',2)
        xlim(Xpl)
        grid on
        legend(['K hat(jw) order ',num2str(FDIopt.Ord)])
        ylabel('|K(jw)| [dB]') %N=Kg m/s^2 = (Kg/s) m/s  
        xlabel('Frequency [rad/s]')
        title('Fluid-Memory Model Frequency Response')
        subplot(212)
        semilogx(Wext2,(180/pi)*angle(Kw_hatFD),'--b','LineWidth',2)
        xlim(Xpl)
        grid on
        ylabel('Phase K(jw) [deg]') %N=Kg m/s^2 = (Kg/s) m/s  
        xlabel('Frequency [rad/s]')
        %end

        %% Reconstruction of Added Mass and Damping

        Bhat=real(Kw_hatFD);
        Ahat =imag(Kw_hatFD)./Wext2+Ainf_hat*ones(size(Wext2));

        figure(103) 
        subplot(221)
        semilogx(W,(abs(Ac)),'o-r',Wext2,(abs(Ac_hatFD)),'--b','LineWidth',2)
        xlim(Xpl)
        grid on
        title(['DoF ',num2str(Dof(1)),num2str(Dof(2))])
        legend('Ac(jw)',['Achat(jw) order ',num2str(FDIopt.Ord)])
        ylabel('|Ac(jw)|') %N=Kg m/s^2 = (Kg/s) m/s  
        xlabel('Freq. [rad/s]')
        subplot(223)
        semilogx(W,(180/pi)*angle(Ac),'o-r',Wext2,(180/pi)*angle(Ac_hatFD),'--b','LineWidth',2)
        xlim(Xpl)
        grid on
        ylabel('Phase Ac(jw) [deg]') %N=Kg m/s^2 = (Kg/s) m/s  
        xlabel('Frequency [rad/s]')
        %
        subplot(222)
        semilogx(W,A,'o-r',Wext2,Ahat,'--b','LineWidth',2)
        xlim(Xpl)
        grid on
        legend('A',['Aest'])
        %end
        title(['DoF ',num2str(Dof(1)),num2str(Dof(2))])
        ylabel('Added Mass')
        xlabel('Frequency [rad/s]')
        subplot(224)
        semilogx(W,B,'o-r',Wext2,Bhat,'--b','LineWidth',2)
        xlim(Xpl)
        grid on
        legend('B',['Best'])
        ylabel('Damping')
        xlabel('Frequency [rad/s]')
        %% Ask for a new iteration
%         display('---------------------------------------------------------------------')
        KK=input(['Current Order ',num2str((FDIopt.Ord)),', (Change or 0-to exit):']);
        display(' ')
        if KK==0,
            break
        end
        FDIopt.Ord=KK;
    end
end
%% Return Result
KradNum=KsNum;
KradDen=KsDen;

