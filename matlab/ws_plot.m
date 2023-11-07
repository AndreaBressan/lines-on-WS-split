function ws_plot(L)
L=L';
xs=-L(:,3)./L(:,1); xgood=(0<xs & xs<=1);
ys=-L(:,3)./L(:,2); ygood=(0<=ys & ys<1);
ds=(-L(:,3)-L(:,2))./(L(:,1)-L(:,2)); dgood=(0<=ds & ds<1);

e1=xgood&ygood;
x1=[zeros(sum(e1),1),xs(e1)]';
y1=[ys(e1),zeros(sum(e1),1)]';
e2=(xgood&dgood) & ~e1;
x2=[ds(e2),xs(e2)]';
y2=[1-ds(e2),zeros(sum(e2),1)]';
e3=(ygood&dgood) & ~(e1 |e2);
x3=[zeros(sum(e3),1),ds(e3)]';
y3=[ys(e3),1-ds(e3)]';
plot([x1 x2 x3],[y1 y2 y3]);
end