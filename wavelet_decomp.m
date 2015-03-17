%%   Wavelet Decomposition
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

function noise_decomp=wavelet_decomp(Data,scales,translation,wavlib)
scalogram=cwt_extract(Data,wavlib);

    for h=1:length(scales)
        for g=1:length(Data)
            noise_decomp(h,g)=sum(scalogram(h,:).*real_morlet((translation-translation(g)),scales(h),2));
        end
    end