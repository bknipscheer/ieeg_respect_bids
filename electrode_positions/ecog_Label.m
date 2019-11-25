function ecog_Label(els,varargin)
%     Copyright (C) 2009  K.J. Miller
%                   2019  D. van Blooijs

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
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% adaptations:
% Sept 2019: shows electrode label instead of number

if isempty(varargin)
    name = num2cell(1:size(els,1));
    msize = 20; % marker size
    fsize = 10; % text font size
else
    if length(varargin)>=1
        name = varargin{1};
        msize = 20;
        fsize = 10;
    end
    if length(varargin)>=2
        msize = varargin{2};
        fsize = 10;
    end
    if length(varargin)>=3
        msize = varargin{2};
        fsize = varargin{3};
    end
    
end

hold on, plot3(els(:,1),els(:,2),els(:,3),'.','MarkerSize',msize,'Color',[.99 .99 .99])
for k=1:length(els(:,1))
    text(els(k,1)*1.05,els(k,2)*1.05,els(k,3)*1.05,name{k},'FontSize',fsize,'HorizontalAlignment','center','VerticalAlignment','middle')
end



