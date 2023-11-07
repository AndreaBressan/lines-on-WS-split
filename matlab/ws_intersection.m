function [p_pts] = ws_intersection(L1,L2)
% intersect all pairs of lines
L1r=repmat(reshape(L1,3,1,[]),1,size(L2,2),1);
L2r=repmat(reshape(L2,3,[],1),1,1,size(L1,2));
p_pts=cross(L1r,L2r,1);
f1=gcd(p_pts(1,:),p_pts(2,:));
f2=gcd(f1,p_pts(3,:));
p_pts=p_pts(:,:)./f2;
p_pts=p_pts.*sign(p_pts(3,:));
end