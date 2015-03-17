%%   Morlet Wavelet
%    Copyright (C) 2014  Samuel Akwei-Sekyere
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
%
%Inputs ; t represents the time stamp
%         centf represents the center frequency
%         bandw represents the frequency bandwidth
%
%Outputs ; comp_morl represents the complex morlet wavelet
%
%Notes: f(a)=f(c)/(a*s) ; where f(a) is the pseudo-frequency
%                      corresponding to scale a
%                      f(c) is the center frequency of the wavelet
%                      a is the scale and s is the sampling period
%                      (ie., s=1/sampling rate)

function comp_morl=complex_morlet(t,centf,bandw)

c=(1+exp(-((centf).^2))-(2.*exp(-((3/4).*(centf).^2)))).^(-1/2);
comp_morl=c.*(1/sqrt((pi.*bandw))).*exp(2.*1i.*pi.*centf.*t).*exp(-t.^2./bandw);  %!Normalized