function [hlp] = firtest(htf,N)
%%% fir linear phase %%%
% do the FFT trick
       % FFT length (even!)
H1 = fft(htf,2*N);
hlp = real(ifft(abs(H1)));
hlp=circshift(hlp,N);
end