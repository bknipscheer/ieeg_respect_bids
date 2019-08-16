%Parse annotation of the form:
%KEYWORD{ChannelName,ChannelSubset,(*ChannelName),(*ChannelSubset)} where
%ChannelName= i.e. Gr1 ChannelSubset= i.e. Gr[1:5]
%
%input_str - string to parse 
%chs       - cell contaings the current channels in use 
%ch2Use    - logical array containing which channel (of ch) is present on
%            the input_str string to parse
%TODO: checking for exception 


%     Copyright (C) 2019 Matteo Demuru
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.


function [ch2use]=parse_annotation(input_str,chs)


if(~contains(input_str,';'))
    error('Wrong annotation format: %s',input_str)
end

C=strsplit(input_str,';');
C=C(~cellfun(@isempty,C));
ch2use=zeros(size(chs));

for i=2:numel(C)
    curr_str=C{i};
    if(~isempty(curr_str))
        if(contains(curr_str,'[')) %it is a ChannelSubset
            ch_subset_str=parse_ch_subset(curr_str,chs);
            ch2use_temp = zeros(size(chs));
            for chan = 1:size(ch_subset_str,1)
                ch2use_temp(:,chan) = strcmpi(chs,ch_subset_str{chan});
            end
        else % ChannelName
             ch2use_temp = strcmpi(chs,curr_str);
        end

    else
        error('Wrong annotation format: %s',input_str)
    end
    ch2use = ch2use + sum(ch2use_temp,2);
end



