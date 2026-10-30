function out=MakeStable(p)
%This function checks if a polynomial has roots in the 
%closed right half plane, and if so, the unstable roots are 
%reflected about the imaginary axis.
r=roots(p);
for k=1:length(r)
        if  real(r(k)) > 0,
            r(k) = -r(k)
        end
end
out=poly(r);