%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Two-dimensional Fourier transform                                      %
%                                                                         %
%                                                                         %
%  (c) Christopher W. Roy, 2018-12-04                                     %
%  fetal.xcmr@gmail.com                                                   %
%  https://github.com/cwroy/Fetal-XCMR/                                   %
%  Roy, C.W. et al. Fetal XCMR: a numerical phantom for fetal             %
%  cardiovascular magnetic resonance imaging. Journal of Cardiovascular   %
%  Magnetic Resonance 21, 29 (2019).                                      %
%  https://doi.org/10.1186/s12968-019-0539-2                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function x = fft2c(x)

fctr = size(x,1)*size(x,2);
x = 1/sqrt(fctr)*fftshift(fft2(ifftshift(x)));

end