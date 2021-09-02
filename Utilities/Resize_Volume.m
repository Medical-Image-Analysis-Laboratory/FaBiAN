%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This script shrinks or expands a multi-dimensional image for a         %
%  desired size.                                                          %
%                                                                         %
%                                                                         %
%  Christopher W. Roy, 2018-12-04                                         %
%  fetal.xcmr@gmail.com                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [I,y,x] = Resize_Volume(I,sz)

if size(I,1)~=sz(1)||size(I,2)~=sz(2)||size(I,3)~=sz(3)
    
    if ~exist('sz','var')
        sz=round(size(I,1)/2);
    end
    if length(sz)==1
        sz=[sz,sz];
    end
    
    if size(I,1)>=sz(1)
        y=(size(I,1)-sz(1))/2;
        y=ceil(y+1:y+sz(1));
        I=I(y,:,:,:);
    else
        p1=sz(1)-size(I,1);
        I=cat(1,zeros(floor(0.5*p1),size(I,2),size(I,3),size(I,4)),I,zeros(ceil(0.5*p1),size(I,2),size(I,3),size(I,4)));
    end
    
    
    if size(I,2)>=sz(2)
        x=(size(I,2)-sz(2))/2;
        x=ceil(x+1:x+sz(2));
        I=I(:,x,:,:);
    else
        p2=sz(2)-size(I,2);
        I=cat(2,zeros(size(I,1),floor(0.5*p2),size(I,3),size(I,4)),I,zeros(size(I,1),ceil(0.5*p2),size(I,3),size(I,4)));
    end
    
    if size(I,3)>1
        if size(I,3)>=sz(3)
            x=(size(I,3)-sz(3))/2;
            x=ceil(x+1:x+sz(3));
            I=I(:,:,x,:);
        else
            p3=sz(3)-size(I,3);
            I=cat(3,zeros(size(I,1),size(I,2),floor(0.5*p3),size(I,4)),I,zeros(size(I,1),size(I,2),ceil(0.5*p3),size(I,4)));
        end
    end
end

end