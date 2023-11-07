function [m,pts]=ws(d,all,fig)
if nargin<2; all=false; end
if nargin<3; fig=false; end
% find all line intersections and count multiplicities
% the computations are done on the triangle 
%  T=[0 0;1 0;0 1]'
% that is roughly
%  |\
%  |_\
% and intersections are computed only in the region
% U :=U1 cup U2
% U1:= x<=1/3 and y<=x
% U2:= 1/3<=x<=1/2 and y<=1-2 x
%
% Lines are described as a vector of coefficients:
% [a; b; c] represents the zeros of the equation a*x+b*y+c=0
% so that the intersection between two lines is, in projective coordinates,
% the cross product of the two lines' coefficient-vectors.
% 
% All computations are done in integers by manually taking care of the
% denominator

% To list all lines, first we build the list them between the vertical
% and orizontal edge and then we rotate the set twice by the rotation
xc=(1:d-1)';yc=1:d;cc=-xc.*yc;
VH=[reshape(repmat(xc*d,numel(yc),1),[],1), reshape(repmat(yc*d,numel(xc),1),[],1), cc(:) ]';
R=@(L)[-1 -1 1; 1 0 0 ; 0 0 1]'*L;
VD=R(VH);
HD=R(VD);
% Then we restrict all sets of lines to only those passing through the
% region U described above
g1= ((VH(1,:)+VH(2,:))>=-3*VH(3,:) ) & (VH(1,:)~=VH(2,:) | VH(2,:)~=-VH(3,:) );
L1=VH(:,g1);
% v-edge -> d-edge
g2=(VD(1,:)+VD(2,:))<=-3*VD(3,:) & (VD(3,:)~=0 | VD(1,:)~=0);
L2=VD(:,g2);
% h-edge -> d-edge
g3= 2*HD(3,:)<=-HD(1,:);
L3=HD(:,g3);
% Having the set of all lines
L=cat(2,L1,L2,L3);
% we find all intersection points inside the region U
int=ws_intersection(L,L);
good=int(1,:)>0 & int(2,:)>0 & int(3,:)>0;
good=good & ( (3*int(1,:)<=int(3,:) & int(2,:)<=int(1,:))...
    |(2*int(1,:)<=int(3,:) & int(2,:)<= int(3,:)-2*int(1,:)));
int=int(:,good);
% take the unique set and count multiplicities taking into account that 
% if k lines contain a point p then  
[pts,~,po]=unique(int','rows');
l=floor(3*d/2);
comb=(1:l).*(2:l+1);
combInv(comb)=(2:l+1)';
m=combInv(accumarray(po,1));
[m,ord]=sort(m,'descend');
pts=pts(ord,:)';
if ~all
pts=pts(:,m==m(1));
m=m(1);
end
if fig
figure;
hold on
ws_plot(L);
scatter(pts(1,:)./pts(3,:), pts(2,:)./pts(3,:));
hold off;
end

end