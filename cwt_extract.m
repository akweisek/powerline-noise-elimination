%%   Wavelet Transformation (part of the pseudo-convolution)
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
%Inputs: data -- the data ; it should be a single row vector
%        wav_m-- a matrix of wavelet at various scales and translations, use wav_m.m to generate this data (mxnxp matrix) m=scale,n=time,p=translation
%		 sample_rate -- the sampling rate of the data
%
%Output: CWT_data -- wavelet transformation of the data

function [CWT_data]=cwt_extract(data,wav_m)

%Alternatively,
%data_new=permute(repmat(repmat(data,length(wav_m(:,1,1)),1),1,1,length(wav_m(1,1,:))),[1,2,3]);
%data_new=data_new.*wav_m;
%CWT_data=permute(sum(permute(data_new,[2,1,3])),[2,3,1]);

%Alternatively,
%CWT_data=zeros(length(wav_m(:,1,1)),length(wav_m(1,1,:)));
%for i=1:length(wav_m(:,1,1))    %scale, a
%	for j=1:length(wav_m(1,1,:))    %translation, b
%		CWT_data(i,j)=sum(wav_m(i,:,j).*data);   %multiply the mother wavelet at different scales and translations by the data and sum the result
%    end                                                                  
%end

CWT_data=bsxfun(@times,wav_m,data);
CWT_data=permute(sum(permute(CWT_data,[2,1,3])),[2,3,1]);