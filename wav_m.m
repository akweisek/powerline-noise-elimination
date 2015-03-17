%%   Wavelet Library Creation
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

%Inputs: scales -- row vector of scales to use
%        translations -- row vector of translations to use
%        wavelet -- which wavelet to use
%
%Output: wavelet_lib -- the library of wavelets used


function [wavelet_lib]=wav_m(scales,translations,time)


wavelet_lib=zeros(length(scales),length(time),length(translations));
h=waitbar(0,'Wavelet Library... \newline \newlineCopyright (C) 2015  Samuel Akwei-Sekyere');
p=0;
q=length(scales)*length(translations);
for i=1:length(scales)
	for j=1:length(translations)
		wavelet_lib(i,:,j)=real_morlet(time-translations(j),scales(i),2);    %Creating a library using Morlet wavelet
        p=p+1;
        waitbar(p/q)
    end
end
close(h)
