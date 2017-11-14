%% Capillary bridge


%% Initial info
discr=0.5/71.5;
%R=0.0025; %radius of sphere in m
gamma=0.07; %  superficial tension (N/m)
profile_up=250; % upper bridge profile coordinate y (better ~10-20 pix lower value) 250
profile_dn=650; % lower bridge profile coordinate y, (10-20 pix upper value)510
% profile_up and profile_dn should be outside of the profile
profile_accuracy_up=3;% accuracy of upper bridge profile (distance sphere- bridge profile)
profile_accuracy_dn=1;% accuracy of lower bridge profile (distance sphere- bridge profile)
circ_up=100; %beginning of searching of upper circle
circ_dn=800; %end of searching of lower circle
pw=1; %bridge approx at contact points: 1- polynomial 6th degr, 2-polynomial 2nd degr, 3- linear
nfl=[1:3];       % image file numbers 
dir=['']; % directory
file=['Image00100.jpg']; % original file
file1=['binary100.tif']; % binarized file

res(1,1:29)={'file', 'dist', 'y*=rg', 'hf angle', 'cont ang', 'condiam_up',...
    'condiam_dn', 'curvat', 'deltaP', 'FdeltaP', 'FST', 'Fcap', 'profile',...
    'min_crit', 'max_crit', 'Vol_int', 'Vol_nodoid_c', 'Vol_nodoid_s',...
    'Vol_Guldin', 'ERR_max', 'dt_ul', 'dt_ur', 'dt_dl', 'dt_dr', 'th_up', 'th_dn', 'VGul_left', 'VGul_right', 'Rayon'}; %head of results table
Vol_integral=0;
H=0;
ERR_max=0;
prof='not_found';


%% File
%s0=input('Single image? (y/n)','s');
s0='y';
if(s0=='y'),
    nfl=1;
end

for ifl=1:length(nfl),
    s2=num2str(nfl(ifl));
    
    if(s0=='n'),
        if(nfl(ifl)<10),
    file=[file(1:end-5) s2 file(end-3:end)];  % original image file name
    file1=[file1(1:end-7) s2 file1(end-5:end)];  % transformed to binary image   
    elseif(nfl(ifl)<100),
    file=[file(1:end-6) s2 file(end-3:end)];  % original image file name
    file1=[file1(1:end-8) s2 file1(end-5:end)];  % transformed to binary image  
        elseif(nfl(ifl)<1000),
          file=[file(1:end-7) s2 file(end-3:end)];  % original image file name
    file1=[file1(1:end-9) s2 file1(end-5:end)];  % transformed to binary image  
        elseif(nfl(ifl)>=1000),
          file=[file(1:end-8) s2 file(end-3:end)];  % original image file name
    file1=[file1(1:end-10) s2 file1(end-5:end)];  % transformed to binary image  
        end
    end

tiffobj=Tiff([dir file1],'r');
imagedata=tiffobj.read();
tiffobj.close();
[nx, ny]=size(imagedata); % image size in pixels
[hx, hy]=find(imagedata~=0); % search of profiles (points differents than 0)


%% Searching of spheres profile

ar_c1=[]; %array of upper sphere points (x,y)
ar_c2=[]; %array of lower sphere points (x,y)
row=1; %row of array
for i=circ_up:profile_up % points of upper sphere
    h=find(hx==i); % choose of row 
    hm=find(hy(h)==min(hy(h))); % left side of upper sphere for row i
    hM=find(hy(h)==max(hy(h))); % right side of upper sphere
    
    if(~isempty(hm) & ~isempty(hM)) % if both sides exist,
        ar_c1(2*row-1,1)=hy(h(hm)); % adding coordinates of left and right profile
        ar_c1(2*row-1,2)=i;
        ar_c1(2*row,1)=hy(h(hM));
        ar_c1(2*row,2)=i;
    end
     row=row+1; % next row
end

row=1;

for i=profile_dn:circ_dn % the same fo  r lower sphere
    h=find(hx==i); 
    hm=find(hy(h)==min(hy(h)));
    hM=find(hy(h)==max(hy(h)));
    
    if(~isempty(hm) & ~isempty(hM))  
        ar_c2(2*row-1,1)=hy(h(hm));
        ar_c2(2*row-1,2)=i;
        ar_c2(2*row,1)=hy(h(hM));
        ar_c2(2*row,2)=i;
    end
    row=row+1;
end


%% Approximation of spheres: centers and radii
%--------------------------------------------------------------------------
% 
%     Algebraic circle fit by Taubin
%      G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
%                  Space Curves Defined By Implicit Equations, With
%                  Applications To Edge And Range Image Segmentation",
%      IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this is a version optimized for stability, not for speed
%
% http://people.cas.uab.edu/~mosya/cl/MATLABcircle.html
%--------------------------------------------------------------------------

centroid = mean(ar_c1);   % the centroid of the data set
X = ar_c1(:,1) - centroid(1);  %  centering data
Y = ar_c1(:,2) - centroid(2);  %  centering data
Z = X.*X + Y.*Y;
Zmean = mean(Z);
Z0 = (Z-Zmean)/(2*sqrt(Zmean));
ZXY = [Z0 X Y];
[U,S,V]=svd(ZXY,0);
A = V(:,3);
A(1) = A(1)/(2*sqrt(Zmean));
A = [A ; -Zmean*A(1)];
Par1 = [-(A(2:3))'/A(1)/2+centroid , sqrt(A(2)*A(2)+A(3)*A(3)-4*A(1)*A(4))/abs(A(1))/2];
% Par1(1), Par1(2)= x,y coordinates of center, Par1(3) = radius, upper sphere

centroid = mean(ar_c2);   % the centroid of the data set
X = ar_c2(:,1) - centroid(1);  %  centering data
Y = ar_c2(:,2) - centroid(2);  %  centering data
Z = X.*X + Y.*Y;
Zmean = mean(Z);
Z0 = (Z-Zmean)/(2*sqrt(Zmean));
ZXY = [Z0 X Y];
[U,S,V]=svd(ZXY,0);
A = V(:,3);
A(1) = A(1)/(2*sqrt(Zmean));
A = [A ; -Zmean*A(1)];
Par2 = [-(A(2:3))'/A(1)/2+centroid , sqrt(A(2)*A(2)+A(3)*A(3)-4*A(1)*A(4))/abs(A(1))/2];
% Par2(1), Par2(2)= x,y coordinates of center, Par2(3) = radius, lower sphere

dist=sqrt((Par2(1)-Par1(1))^2+(Par2(2)-Par1(2))^2)-Par2(3)-Par1(3); % distance between spheres
r=(Par1(3)+Par2(3))/2; % mean radii

%pixel to mm
%discr=R/r;
dist_m=dist*discr;
r_m=r*discr;
r1_m=Par1(3)*discr;
r2_m=Par2(3)*discr;

r1 = Par1(3);
r2 = Par2(3);



%% Plotting of both spheres
f1=figure;
A1=imread([dir file]);
colormap(gray);
imagesc(A1); 
axis equal; 
axis tight;
hold on;
plot(Par1(1), Par1(2), 'bx'); % center of upper sphere
plot(Par2(1), Par2(2), 'bx'); % center of lower sphere
viscircles([Par1(1), Par1(2)], Par1(3), 'EdgeColor','r', 'LineWidth',1, 'LineStyle',':'); % upper sphere profile
viscircles([Par2(1), Par2(2)], Par2(3), 'EdgeColor','r', 'LineWidth',1,'LineStyle',':'); % lower sphere profile


%% Capillary bridge

xl=[]; % x points of left profile
xr=[]; % x points of right profile
yl=[]; % y points of left profile
yr=[]; % y points of right profile

for i=profile_up:profile_dn % capillary bridge
    H=find(hx==i); % row x
    Hm=find(hy(H)==min(hy(H))); % minimum (left side) of given row
    HM=find(hy(H)==max(hy(H))); % maximum (right side) of given row
    if(~isempty(hm) & ~isempty(hM)),
        xl=[xl hy(H(Hm))]; % x points of left profile
        yl=[yl i]; % y points of left profile
        xr=[xr hy(H(HM))]; % x points of right profile
        yr=[yr i]; % y points of right profile
    end
end


% sam profil : dla wyliczenia na podstawie profilu
% Par1(1) wspolrzedna x srodka kulki 1, Par1(2) wspolrzedna y,
% Par1(3) promien kulki 1
% Par2(1) wspolrzedna x srodka kulki 1, Par2(2) wspolrzedna y,
% Par2(3)promien kulki 1




hx1=find((xl-Par1(1)).^2+(yl-Par1(2)).^2-(Par1(3)+profile_accuracy_up)^2>0 & ...
    (xl-Par2(1)).^2+(yl-Par2(2)).^2-(Par2(3)+profile_accuracy_dn).^2>0); % search of left bridge profile (out of spheres) 
hx2=find((xr-Par1(1)).^2+(yr-Par1(2)).^2-(Par1(3)+profile_accuracy_up).^2>0 & ...
    (xr-Par2(1)).^2+(yr-Par2(2)).^2-(Par2(3)+profile_accuracy_dn).^2>0); % search of right bridge profile (out of spheres)

X1=smooth(xl(hx1),0.25,'moving'); % smoothing of left profile
X2=smooth(xr(hx2),0.25,'moving'); % smoothing of right profile

[yy i1 i2]=intersect(yl(hx1),yr(hx2)); % find common set of y-coordinates for left and right profiles 
hs_cvx=find(X2(i2)-X1(i1)==min(X2(i2)-X1(i1))); % find the points with the minimal distance between the right and left bridge profile sides
y_str_cvx=min(X2(i2)-X1(i1))/2;   % neck radius as a minimal distance between two sides of the bridge profile (convex bridge)
hs_ccv=find(X2(i2)-X1(i1)==max(X2(i2)-X1(i1)));
y_str_ccv=max(X2(i2)-X1(i1))/2; % neck radius as a maximal distance between two sides of the bridge profile (concave bridge)


%y_nk=mean(yy(hs));



%% Contact points

pt1=[xl(hx1(1)) yl(hx1(1))]; % up left
pt2=[xr(hx2(1)) yr(hx2(1))]; % up right
pt3=[xl(hx1(end)) yl(hx1(end))]; % low left
pt4=[xr(hx2(end)) yr(hx2(end))]; % low right


% plot of contact points
plot(pt1(1),pt1(2),'rx','MarkerSize', 10, 'LineWidth', 2);
plot(pt2(1),pt2(2),'rx','MarkerSize', 10, 'LineWidth', 2);
plot(pt3(1),pt3(2),'rx','MarkerSize', 10, 'LineWidth', 2);
plot(pt4(1),pt4(2),'rx','MarkerSize', 10, 'LineWidth', 2);


% half-filling angle delta (dt)

condiam_up=sqrt((pt2(1)-pt1(1))^2+(pt2(2)-pt1(2))^2); % upper contact diameter
condiam_dn=sqrt((pt4(1)-pt3(1))^2+(pt4(2)-pt3(2))^2); % lower contact diameter
condiam_up_m=condiam_up*discr; 
condiam_dn_m=condiam_dn*discr;
% dt_up=asin(condiam_up/(2*Par1(3))); % upper half-filling angle, mean
% dt_dn=asin(condiam_dn/(2*Par2(3))); % lower half-filing angle, mean
% dt=mean([dt_up dt_dn])/pi*180; % mean half-filling angle

dt1=atan((Par1(1)-pt1(1))/(pt1(2)-Par1(2))); % left upper delta
dt2=atan((pt2(1)-Par1(1))/(pt2(2)-Par1(2))); % right upper delta
dt3=atan((Par2(1)-pt3(1))/(Par2(2)-pt3(2))); % left lower delta
dt4=atan((pt4(1)-Par2(1))/(Par2(2)-pt4(2))); % right lower delta
dt_up=mean([dt1 dt2]);
dt_dn=mean([dt3 dt4]);
d1 = dt_up; % contact angle for the first sphere
d2 = dt_dn; % contact angle for the second sphere
dt=mean([dt1 dt2 dt3 dt4]); % mean delta
dt_deg=dt/pi*180; % delta mean in degrees
d1_deg=radtodeg(d1)
d2_deg=radtodeg(d2)

% % Gorge radius, contact diameter
% if (y_str_cvx<condiam_up/2 | y_str_cvx<condiam_dn/2) & (y_str_ccv<condiam_up/2 | y_str_ccv<condiam_dn/2)
%     y_str=y_str_cvx;
%     x_nk=mean(X1(i1(hs_cvx)))+y_str; % -||- x-axis wspolrzedna x srodka mostka = wspolrz lewa dla minimum profilu + prom gorge
%     ynk_l=mean(yl(find(X1==max(X1)))); % coordinate y of left maximum (at gorge)
%     ynk_r=mean(yr(find(X2==min(X2)))); % coordinate y of right minimum (at gorge)
%     %y_nk=mean(yl(hs_cvx)); % wspolrzedna y srodka = wartosc y dla minimum profilu
% elseif (y_str_ccv>condiam_up/2 | y_str_ccv>condiam_dn/2) & (y_str_cvx>condiam_up/2 | y_str_cvx>condiam_dn/2)
%     y_str=y_str_ccv;
%     x_nk=mean(X1(i1(hs_ccv)))+y_str; % wspolrzedna x srodka mostka = wspolrz lewa dla minimum profilu + prom gorge
%     ynk_l=mean(yl(find(X1==min(X1)))); % coordinate y of left maximum (at gorge)
%     ynk_r=mean(yr(find(X2==max(X2)))); % coordinate y of right minimum (at gorge)
%     %y_nk=mean(yl(hs_ccv)); % wspolrzedna y srodka = wartosc y dla minimum profilu
% else
%     y_str=y_str_cvx;
% end

if (mean(X1(i1(hs_cvx))) > mean(pt1(1), pt3(1))) & (mean(X2(i2(hs_cvx))) < mean(pt2(1), pt4(1)))
   pont='convex';
   y_str=y_str_cvx;
   x_nk=mean(X1(i1(hs_cvx)))+y_str; % -||- x-axis wspolrzedna x srodka mostka = wspolrz lewa dla minimum profilu + prom gorge
   %ynk_l=mean(yl(find(X1==max(X1)))); % coordinate y of left maximum (at gorge)
   %ynk_r=mean(yr(find(X2==min(X2)))); % coordinate y of right minimum (at gorge)
   %y_nk=mean(yl(hs_cvx)); % wspolrzedna y srodka = wartosc y dla minimum profilu
elseif (mean(X1(i1(hs_ccv))) < mean(pt1(1), pt3(1))) & (mean(X2(i2(hs_ccv))) > mean(pt2(1), pt4(1)))
    pont='concave';
    y_str=y_str_ccv;
    x_nk=mean(X1(i1(hs_ccv)))+y_str; % wspolrzedna x srodka mostka = wspolrz lewa dla minimum profilu + prom gorge
%     ynk_l=mean(yl(find(X1==min(X1)))); % coordinate y of left maximum (at gorge)
%     ynk_r=mean(yr(find(X2==max(X2)))); % coordinate y of right minimum (at gorge)
    %y_nk=mean(yl(hs_ccv)); % wspolrzedna y srodka = wartosc y dla minimum profilu
else
    pont='cylindre';
    y_str=y_str_cvx;
    
end







y_nk=((Par1(2)+Par1(3))+(Par2(2)-Par2(3)))/2;  % position of the neck center along y-axis => in the middle between particles
%y_nk=mean([ynk_l, ynk_r]);


plot(x_nk, y_nk, 'rx'); % plot of bridge "center" x_nk, y_nk



%% Contact angle theta
% approximation by 2nd polynom of all the points of the left (p1) and right 
% (p2) sides of the bridge profiles: x(y)=p(1)*y^2+p(2)*y+p(3)



%plotting approximations
% th=0.3101;
% xc=y_nk-pt1(2);
% a1=(xc*cot(dt+th)-2*(r*sin(dt)-y_str))/(xc^3);
% b1=(3*(r*sin(dt)-y_str)-xc*cot(dt+th))/(xc^2);
% c1=y_str;
% 
% plot(x_nk-((a1*(y_nk-yl(hx1)).^3+b1*(y_nk-yl(hx1)).^2+c1)),yl(hx1),'b-','linewidth',2);
% 
% xc=y_nk-pt2(2);
% a1=(xc*cot(dt+th)-2*(r*sin(dt)-y_str))/(xc^3);
% b1=(3*(r*sin(dt)-y_str)-xc*cot(dt+th))/(xc^2);
% c1=y_str;
% 
% plot(x_nk+((a1*(y_nk-yr(hx1)).^3+b1*(y_nk-yr(hx1)).^2+c1)),yr(hx1),'b-','linewidth',2);

% polynomial fit 6th degree
if(pw==1),
    
[P1,S1,mu1]=polyfit(yl(hx1),xl(hx1),6);
[P2,S2,mu2]=polyfit(yr(hx2),xr(hx2),6);
plot(polyval(P1,yl(hx1),S1,mu1),yl(hx1),'g-','linewidth',1);
plot(polyval(P2,yr(hx2),S2,mu2),yr(hx2),'g-','linewidth',1);    
      
    %function retval = polyfit_convert(p2, x)
% convert 3 return parameter polynomial fit vector to the type returned if
% one return value is requested:
%
% p1 = polyfit(x,y,n);
% [p2,S,mu] = polyfit(x,y,n);
% p3 = polyfit_convert(p2);
% => p1 == p3
n = numel(P1)-1;
m = mean(yl(hx1));
s = std(yl(hx1));

rv1 = zeros(size(P1));
for i = [0:n];
for j = [0:i];
rv1(n+1-j) = rv1(n+1-j) + P1(n+1-i)*nchoosek(i, j)*(-m)^(i-j)/s^i;
end
end


n = numel(P2)-1;
m = mean(yr(hx2));
s = std(yr(hx2));

rv2= zeros(size(P2));
for i = [0:n];
for j = [0:i];
rv2(n+1-j) = rv2(n+1-j) + P2(n+1-i)*nchoosek(i, j)*(-m)^(i-j)/s^i;
end
end

    tb1=[-(6*rv1(1)*pt1(2)^5+5*rv1(2)*pt1(2)^4+4*rv1(3)*pt1(2)^3+3*rv1(4)*pt1(2)^2+2*rv1(5)*pt1(2)+rv1(6)),-1];
    tb2=[-(6*rv2(1)*pt2(2)^5+5*rv2(2)*pt2(2)^4+4*rv2(3)*pt2(2)^3+3*rv2(4)*pt2(2)^2+2*rv2(5)*pt2(2)+rv2(6)),-1];
    tb3=[(6*rv1(1)*pt3(2)^5+5*rv1(2)*pt3(2)^4+4*rv1(3)*pt3(2)^3+3*rv1(4)*pt3(2)^2+2*rv1(5)*pt3(2)+rv1(6)),1];
    tb4=[(6*rv2(1)*pt4(2)^5+5*rv2(2)*pt4(2)^4+4*rv2(3)*pt4(2)^3+3*rv2(4)*pt4(2)^2+2*rv2(5)*pt4(2)+rv2(6)),1];

    % polynomial fit 2nd degree
end
    
if(pw==2),
    

    
    %approximation by a second order polynom: x(y)=p(1)*y^2+p(2)*y +p(3)
    p1=polyfit(yl(hx1),xl(hx1),2);
    p2=polyfit(yr(hx2),xr(hx2),2);
    
    %plotting approximations
    plot(p1(1)*yl(hx1).^2+p1(2)*yl(hx1)+p1(3),yl(hx1),'b-','linewidth',2);
    plot(p2(1)*yr(hx2).^2+p2(2)*yr(hx2)+p2(3),yr(hx2),'b-','linewidth',2);
    
    % directional vector tb for each profile end at triple point pt: 
    % tb_x=(p(1)*pt_y+p(2))*dy, tb_y=dy; with dy idem. as for linear approximation
    
   
    tb1=[-(2*p1(1)*pt1(2)+p1(2)) -1]; % upper left
    tb2=[-(2*p2(1)*pt2(2)+p2(2)) -1]; % upper right
    tb3=[(2*p1(1)*pt3(2)+p1(2)) 1];   %lower left
    tb4=[(2*p2(1)*pt4(2)+p2(2)) 1];   % lower right
    
end

% linear fit
if(pw==3),

    nip1=fix(length(hx1)/40);
    nip2=fix(length(hx2)/40);
    
    % linear approximations x(y)=p(1)*y+p(2)
    p1=polyfit(yl(hx1(1:nip1)),xl(hx1(1:nip1)),1);
    p2=polyfit(yr(hx2(1:nip2)),xr(hx2(1:nip2)),1);
    p3=polyfit(yl(hx1(end-nip1:end)),xl(hx1(end-nip1:end)),1);
    p4=polyfit(yr(hx2(end-nip2:end)),xr(hx2(end-nip2:end)),1);
    
    % plotting liner approximation
plot(p1(1)*yl(hx1(1:nip1))+p1(2),yl(hx1(1:nip1)),'c-','linewidth',2);
plot(p2(1)*yr(hx2(1:nip2))+p2(2),yr(hx2(1:nip2)),'c-','linewidth',2);
plot(p3(1)*yl(hx1(end-nip1:end))+p3(2),yl(hx1(end-nip1:end)),'c-','linewidth',2);
plot(p4(1)*yr(hx2(end-nip2:end))+p4(2),yr(hx2(end-nip2:end)),'c-','linewidth',2);
 
    % directional vector tb (non normalized): tb_x=x'(y)*dy=p(1)*dy, tb_y=dy;
    % take dy=+/-1 to get a directional vector pointing from the triple point "outside the profile";
    % Remember: the y-axis is reversed, thus -1 for upper and 1 for lower
    tb1=[-(p1(1)) -1]; % upper left
    tb2=[-(p2(1)) -1]; % upper right
    tb3=[(p3(1)) 1];   % lower left
    tb4=[(p4(1)) 1];   % lower right
end


% normalisation of the profile directional vectors
tb1=tb1/norm(tb1);
tb2=tb2/norm(tb2);
tb3=tb3/norm(tb3);
tb4=tb4/norm(tb4);

% directional vectors of the particle surface profile at the triple points; 
% sign is chosen to be pointed as for tb
% tg1=[(pt1(2)-Par1(2)), -(pt1(1)-Par1(1))];
% tg1=tg1/norm(tg1);
% tg2=[(pt2(2)-Par1(2)), pt2(1)-Par1(1)];
% tg2=tg2/norm(tg2);
% tg3=[-(pt3(2)-Par2(2)), -(pt3(1)-Par2(1))];
% tg3=tg3/norm(tg3);
% tg4=[-(pt4(2)-Par2(2)), pt4(1)-Par2(1)];
% tg4=tg4/norm(tg4);

tg1=[-cos(dt1) -sin(dt1)];
tg2=[cos(dt2) -sin(dt2)];
tg3=[-cos(dt3) sin(dt3)];
tg4=[cos(dt4) sin(dt4)];


th1=acos(dot(tb1,tg1));   % upper left
th2=acos(dot(tb2,tg2));   % upper right
th3=acos(dot(tb3,tg3));   % lower left
th4=acos(dot(tb4,tg4));   % lower right
th=mean([th1 th2 th3 th4]); % resulting theta = average over all thetas
theta1=mean([th1 th2]);
theta2=mean([th3 th4]);
th_deg=th/pi*180; %theta in deg
th_up=mean([th1 th2]);
th_dn=mean([th3 th4]);
theta1_deg = radtodeg(theta1);
theta2_deg = radtodeg(theta2);


%% Profile shape estimation

crit_min=r*sin(dt)*sin(dt+th); % lower criterion for convex ndoid: r sin(delta)sin(delta+theta)
crit_max=r*sin(dt); % upper criterion for convex ndoid: r sin(delta)
crit_ccv=r*sin(dt)/sin(dt+th); %criterion for concave profiles
VC1=pi/3*Par1(3)^3*(1-cos(dt_up))^2*(2+cos(dt_up)); %volume of upper spherical cap
VC2=pi/3*Par2(3)^3*(1-cos(dt_dn))^2*(2+cos(dt_dn)); %volume of lower spherical cap

%% ===> NODOID CASE
% Note: x-coordinates here corresponds to y-coordinates in the article
if(y_str>crit_min&crit_max>y_str), %y_str>crit_min&crit_max>y_str
    prof='nodoid';
    denom_cal_nodo = r1*sin(d1)*sin(d1+theta1)-r2*sin(d2)*sin(d2+theta2);
    a  = 0.5*(-r1^2*sin(d1)^2 + r2^2*sin(d2)^2)/denom_cal_nodo;
    b2 = r1*r2*sin(d1)*sin(d2)*(-r1*sin(d1)*sin(d2+theta2) + r2*sin(d2)*sin(d1+theta1))/denom_cal_nodo
    e  = sqrt(a^2+b2)/a;
    tau  = acos(e*(b2-r^2*sin(dt)^2)/(b2+r^2*sin(dt)^2));
    tau1 = acos(e*(b2-r1^2*sin(d1)^2)/(b2+r1^2*sin(d1)^2));
    tau2 = acos(e*(b2-r2^2*sin(d2)^2)/(b2+r2^2*sin(d2)^2));
    H  = 1/(a*discr);
    t  = [-tau1:(tau1+tau2)/200:tau2];
    
    for i=1:length(t)
        % integral part of the nodoid y-coordinate
        ynodoid=@(t0,e)sin(t0).^2./sqrt(e^2-cos(t0).^2);   % simplified form
        %ynodoid=@(t0,e)cos(t0)./((e+cos(t0))*sqrt(e^2-cos(t0).^2));
        % calculation of the integral
        Y_int=integral(@(t0)ynodoid(t0,e),0,t(i),'RelTol',0,'AbsTol',1e-12);
              
        % nodoid y-coordinates
        Ynod(i)=-a*Y_int+a*sin(t(i))*sqrt((e-cos(t(i)))/(e+cos(t(i))));
        %Ynod(i)=(b2/a)*Y_int;
       
        % nodoid x-coordinates
        Xnod(i)=sqrt(b2)*sqrt((e-cos(t(i)))/(e+cos(t(i))));
    end
    
    % plotting the nodoid portions (left and right, +/- Xnod) centered in the
    % neck center
    plot(-Xnod+x_nk,Ynod+y_nk,'r-','linewidth',2);
    plot(Xnod+x_nk,Ynod+y_nk,'r-','linewidth',2);
   
    % volume of liquid by integration of Ynod(x)^2 dx - volume of spherical
    % caps
    funvol=@(tt,e) ((e-cos(tt)).*cos(tt))./((e+cos(tt)).^2.*sqrt(e^2-cos(tt).^2));
    vol_init=integral(@(tt)funvol(tt,e),0,tau,'RelTol',0,'AbsTol',1e-12); % volume integral part
    Vol_integral=2*pi*(b2^2/a)*vol_init-VC1-VC2; % resulting volume of liquid
    
    % --> Estimation of maximal error of the profile approximation by nodoid
    y0=[ min(Ynod) : max(Ynod)];
    X0=interp1(Ynod,Xnod,y0,'linear');
    
    Y1=interp1(yl(hx1)-y_nk,xl(hx1)-x_nk,y0,'linear','extrap');
    Y2=interp1(yr(hx2)-y_nk,xr(hx2)-x_nk,y0,'linear','extrap');
    
    ERR_max=max(max(abs(1+Y1./y0)),max(abs(1-Y2./y0)));
    
end


%% ====> UNDULOID CASE
% Note: x-coordinates here corresponds to y-coordinates in the article
if(y_str<crit_min) %y_str<crit_min
    prof='unduloid';
    denom_cal_undu =  r1*sin(d1)*sin(d1+theta1)-r2*sin(d2)*sin(d2+theta2);
    a=0.5*(r1^2*sin(d1)^2-r2^2*sin(d2)^2)/denom_cal_undu;
    b2=r1*r2*sin(d1)*sin(d2)*(r1*sin(d1)*sin(d2+theta2)-r2*sin(d2)*sin(d1+theta1))/denom_cal_undu;
    e=sqrt(a^2-b2)/a;
    tau=acos(1/e*((b2-r^2*sin(dt)^2)/(b2+r^2*sin(dt)^2)));
    tau1=acos(1/e*((b2-r1^2*sin(d1)^2)/(b2+r1^2*sin(d1)^2)));
    tau2=acos(1/e*((b2-r2^2*sin(d2)^2)/(b2+r2^2*sin(d2)^2)));
    H=-1/abs(a*discr);
    t=[-tau1:(tau1+tau2)/200:tau2];
    for i=1:length(t)
        % integral part of the unduloid y-coordinate
        yunduloid=@(t0,e) 1./(1+e*cos(t0))./sqrt(1-e^2*cos(t0).^2);
        % calculation of the integral
        Y_int=integral(@(t0)yunduloid(t0,e),0,t(i),'RelTol',0.1,'AbsTol',1e-12); % 0.1 and 1e-12
        % unduloid y-coordinates
        Yund(i)=b2/a*Y_int;
        % unduloid x-coordinates
        Xund(i)=sqrt(b2)*sqrt((1-e*cos(t(i)))/(1+e*cos(t(i))));
    end
    % plotting the unduloid portions (left and right, +/- Xund) centered in the
    % neck center
    plot(Xund+x_nk,Yund+y_nk,'m-','linewidth',2);
    plot(-Xund+x_nk,Yund+y_nk,'m-','linewidth',2);
    
    % volume of liquid by integration of Yund(x)^2 dx - volume of spherical
    % caps
    funvolu=@(tt,e) (1-e*cos(tt))./((1+e*cos(tt)).^2.*sqrt(1-e^2*cos(tt).^2));
    vol_initu=integral(@(tt)funvolu(tt,e),0,tau,'RelTol',0,'AbsTol',1e-12);
    Vol_integral=2*pi*(b2^2/a)*vol_initu-VC1-VC2; % resulting volume
    
    
    % Estimation of the profile approximation by the unduloid
    
    y0=[min(Yund):max(Yund)];
    X0=interp1(Yund,Xund,y0,'nearest'); % original value 'linear'
    
    Y1=interp1(yl(hx1)-y_nk,xl(hx1)-x_nk,y0,'linear','extrap');
    Y2=interp1(yr(hx2)-y_nk,xr(hx2)-x_nk,y0,'linear','extrap');
    
    ERR_max=max(max(abs(1+Y1./y0)),max(abs(1-Y2./y0)));
    
end


%% === > CIRCULAR CYLINDER CASE

if(y_str==crit_max&y_str>crit_min)
    prof='circular';
    dt=asin(y_str/r);
    H=1/(y_str*discr*sin(dt));
    plot([-r*sin(dt)+Par1(1) -r*sin(dt)+Par2(1)],[r*cos(dt)+Par1(2) -r*cos(dt)+Par2(2)],'m');
    plot([r*sin(dt)+Par1(1) r*sin(dt)+Par2(2)],[r*cos(dt)+Par1(2) -r*cos(dt)+Par2(2)],'m');
    
end

%% ===> CATENOID CASE
if(y_str==crit_min)
    prof='catenoid';
    y0=-r*sin(dt):0.1:r*sin(dt);
    x0=r*sin(dt)*sin(dt+th)*cosh(y0/(r*sin(dt)*sin(dt+th)));
    H=0;
    plot(x_nk-x0,y0+y_nk,'c','linewidth',2);
    plot(x_nk+x0,y0+y_nk,'c','linewidth',2);
    
end


%% ===> CONCAVE NODOID
if (crit_max<y_str & y_str<crit_ccv)
    prof='Concave nodoid';
    
    a=0.5*(r^2*sin(dt)^2-y_str^2)/(y_str-r*sin(dt)*sin(dt+th));
    b2=y_str*r*sin(dt)*((r*sin(dt)-y_str*sin(dt+th))/(y_str-r*sin(dt)*sin(dt+th)));
    e=sqrt(a^2+b2)/a;
    tau=(acos(e*(b2-r^2*sin(dt)^2)/(b2+r^2*sin(dt)^2)))-pi;
    H=1/(a*discr);
    t=[pi-tau:tau/20:pi+tau];
    
    for i=1:length(t)
        % integral part of the nodoid y-coordinate
        %ynodoid=@(t0,e)sin(t0).^2./sqrt(e^2-cos(t0).^2);   % simplified form
        ynodoid=@(t0,e)cos(t0)./((e+cos(t0)).*sqrt(e^2-cos(t0).^2))
        % calculation of the integral
        Y_int=integral(@(t0)ynodoid(t0,e),0,t(i),'RelTol',0,'AbsTol',1e-12);
              
        % nodoid y-coordinates
        %Ynod(i)=-a*Y_int+a*sin(t(i))*sqrt((e-cos(t(i)))/(e+cos(t(i))));
        Ynod(i)=pi*(b2^2/a)*Y_int;
       
        % nodoid x-coordinates
        Xnod(i)=sqrt(b2)*sqrt((e-cos(t(i)))/(e+cos(t(i))));
    end
    
    % plotting the nodoid portions (left and right, +/- Xnod) centered in the
    % neck center
    plot(-Xnod+x_nk,Ynod+y_nk,'r-','linewidth',1);
    plot(Xnod+x_nk,Ynod+y_nk,'r-','linewidth',1);
   
    % volume of liquid by integration of Ynod(x)^2 dx - volume of spherical
    % caps
    funvol=@(tt,e)(e-cos(tt)./e+cos(tt)).*(cos(tt)./((e+cos(tt).*sqrt(e^2-cos(tt).^2))));;
    vol_init=integral(@(tt)funvol(tt,e),pi-tau,pi+tau,'RelTol',0,'AbsTol',1e-12); % volume integral part
    Vol_integral=pi*(b2^2/a)*vol_init-VC1-VC2; % resulting volume of liquid
    
    % --> Estimation of maximal error of the profile approximation by nodoid
    y0=[ min(Ynod) : max(Ynod)];
    X0=interp1(Ynod,Xnod,y0,'linear');
    
    Y1=interp1(yl(hx1)-y_nk,xl(hx1)-x_nk,y0,'linear','extrap');
    Y2=interp1(yr(hx2)-y_nk,xr(hx2)-x_nk,y0,'linear','extrap');
    
    ERR_max=max(max(abs(1+Y1./y0)),max(abs(1-Y2./y0)));
end
   
%% ===> CONCAVE UNDULOID
if (y_str>crit_ccv)
    prof='Concave unduloid';
    a=0.5*(y_str^2-r^2*sin(dt)^2)/(y_str-r*sin(dt)*sin(dt+th));
    b2=y_str*r*sin(dt)*(y_str*sin(dt+th)-r*sin(dt))/(y_str-r*sin(dt)*sin(dt+th));
    e=sqrt(a^2-b2)/a;
    tau=acos((1/e)*((b2-r^2*sin(dt)^2)/(b2+r^2*sin(dt)^2)))-pi;
    H=-1/abs(a);
    t=[pi-tau:tau/10:tau+pi];
    for i=1:length(t)
        % integral part of the unduloid y-coordinate
        yunduloid=@(t0,e) 1./((1+e*cos(t0)).*sqrt(1-e^2*cos(t0).^2));
        % calculation of the integral
        Y_int=integral(@(t0)yunduloid(t0,e),0,t(i));
        % unduloid y-coordinates
        Yund(i)=(b2/a)*Y_int;
        % unduloid x-coordinates
        Xund(i)=sqrt(b2)*sqrt((1-e*cos(t(i)))/(1+e*cos(t(i))));
    end
    % plotting the unduloid portions (left and right, +/- Xund) centered in the
    % neck center
        plot(Xund+x_nk,Yund,'m-','linewidth',1);
    plot(-Xund+x_nk,Yund,'m-','linewidth',1);
    
     % volume of liquid by integration of Ynod(x)^2 dx - volume of spherical
    % caps
    funvolu=@(tt,e) (1-e*cos(tt))./((1+e*cos(tt)).^2.*sqrt(1-e^2*cos(tt).^2));
    vol_initu=integral(@(tt)funvolu(tt,e),0,tau,'RelTol',0,'AbsTol',1e-12);
    Vol_integral=2*pi*(b2^2/a)*vol_initu-VC1-VC2; % resulting volume
    
    % --> Estimation of maximal error of the profile approximation by nodoid
    y0=[min(Yund):max(Yund)];
    X0=interp1(Yund,Xund,y0,'linear');
    Y1=interp1(yl(hx1)-y_nk,xl(hx1)-x_nk,y0,'linear','extrap');
    Y2=interp1(yr(hx2)-y_nk,xr(hx2)-x_nk,y0,'linear','extrap');
    ERR_max=max(max(abs(1+Y1./y0)),max(abs(1-Y2./y0)));

    
end


%% Volume and forces

% liquid volume estimation by approximate formula for nodoid
xc=dist/2+r*(1-cos(dt)); % triple point position along y  (mm)
% complex formula for volume estimation
Vol_nodoid_c=2*pi*xc*((r*sin(dt)-y_str)^2/5+y_str^2+2/3*y_str*(r*sin(dt)-y_str))-VC1-VC2;

% simplified formula for volume estimation
Vol_nodoid_s=(2*pi*xc/3)*(2*y_str^2+r^2*sin(dt)^2)-VC1-VC2;

%Volume- Pappus-Guldin theorem

mat_X2=[];
mat_X1=[];
area_left=0;
area_right=0;
row=1; 
for i=1:length(hx1)% left bridge profile (out of spheres)
        mat_X1(row)=mean([x_nk,X1(i)]); % mean values of left bridge profile and bridge center
        area_left=area_left+(x_nk-X1(i)); % total left half-profile area
        row=row+1;
end

row=1; 
for i=1:length(hx2)% right bridge profile (out of spheres)
        mat_X2(row)=mean([x_nk,X2(i)]); % mean values of right bridge profile and bridge center
        area_right=area_right+(X2(i)-x_nk); % total right half-profile area 
        row=row+1;
end

centr_left=mean(mat_X1); %left profile centroid (value of x)
centr_right=mean(mat_X2); %right profile centroid (value of x)
r_centr_left=x_nk-centr_left(1); %left radius of rotation for guldin volume
r_centr_right=centr_right(1)-x_nk; %right radius of rotation for guldin volume
V_left=2*pi*r_centr_left*area_left-VC1-VC2; %volume based on left profile
V_right=2*pi*r_centr_right*area_right-VC1-VC2; %volume based on right profile
V_guldin=mean([V_left,V_right]); %mean volume

% Volumes in mm
Vol_integral_m=Vol_integral*discr^3;
Vol_nodoid_c_m=Vol_nodoid_c*discr^3;
Vol_nodoid_s_m=Vol_nodoid_s*discr^3;
V_guldin_m=V_guldin*discr^3;
V_left_m=V_left*discr^3;
V_right_m=V_right*discr^3;

%Capillary force (N)
yst=y_str*discr;         % neck radius in mm
Fcap=pi*gamma*(H*yst^2+2*yst);
deltaP=gamma*H;
FdeltaP=pi*yst^2*deltaP;
FST=2*pi*yst*gamma;


% finalisation of plot
axis([x_nk-400 x_nk+400 y_nk-300 y_nk+300]);
title([file, ' D=', num2str(dist_m) 'mm, y*=' num2str(y_str*discr) 'mm, \delta =' num2str(dt_deg) ', \theta = ' num2str(th_deg) ', ' prof ]);
eval(['print -djpeg90 ' dir 'zoom_' file(1:end-4)]);

% collecte the results for all the cases
res(ifl+1,1:29)={file, dist_m, yst, dt_deg, th_deg, condiam_up_m, condiam_dn_m, H, deltaP, FdeltaP, FST, Fcap, prof, crit_min*discr, crit_max*discr, Vol_integral_m, Vol_nodoid_c_m, Vol_nodoid_s_m, V_guldin_m, ERR_max, dt1/pi*180, dt2/pi*180, dt3/pi*180, dt4/pi*180, th_up/pi*180, th_dn/pi*180, V_left_m, V_right_m, r_m};

%profile_up=pt1(2)-50; % upper bridge profile coordinate y, ~10 pix upper
%profile_dn=pt3(2)+50; % lower bridge profile coordinate y, ~10 pix lower
file
disp([d1_deg,d2_deg,radtodeg(th1),radtodeg(th2), ...
    radtodeg(th3),radtodeg(th4),y_str*discr,dist_m])
%close
end

%xlswrite([dir 'results.xls'],res)
% problems : roznej wiekosci kulki : gdzie srodek mostka? (xn, yn), 
% jakie half-filling angle? (2 rozne)



