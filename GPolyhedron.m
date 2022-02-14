function [V,Vx,Vy,Vz,Vij]=GPolyhedron(G,p,topo,coords)

%% GRAVITY SIGNAL OF A GENERAL POLYHEDRON 
%% USING THE LINE INTEGRAL APPROACH
%
% The function computes analytically the gravitational potential and its up to second 
% order derivatives from a generally shaped homogeneous polyhedral distribution  
% at arbitary points in 3D space, using the line integral approach. The required 
% input data are the constant density of the distribution, Newton's gravitational
% constant, the Cartesian coordinates (X,Y,Z) of the polyhedral vertices with respect   
% to the computation point and the connectivity of the vertices.
%
%
% The function is called using following 4 arguments
%   1. Newton's gravitational constant (G).
%   2. The constant density of the distribution (p).
%   3. The name of the text file containing the topology which defines 
%      the boundary surface of the polyhedron. Each line represents an 
%      independent face (topo).
%   4. The name of the text file containing the known coordinates of all  
%      polyhedral vertices, which enter the topology matrix of file 'topo'.
%      The file contains 3 columns, with the X, Y and Z coordinates of the
%      vertices with respect to the computation point (coords).
%
%----->Example of function use: 
%      [V,Vx,Vy,Vz,Vij]=GPolyhedron(6.67259e-11,2670,'topology.txt','coords.txt');
%      where:
%      G=6.67259e-11[(m^3)*(kg^(-1))*(sec^(-2))],
%      p=2670[Kg*m^(-3)],
%      'topology.txt': the text file containing the connectivity of the vertices,
%      'coords.txt': the text file with the coordinates X[m], Y[m], Z[m].
%
%
% The evaluated parameters are:
%   1. One value for the gravitational potential (V).
%   2. Three values for the fisrt order derivatives of the gravitational 
%      potential (Vx, Vy, Vz).
%   3. A 3x3 matrix of the second order derivatives of the gravitational 
%      potential (Vij).
%
%
% USEFUL NOTES
%
%1. The input text file 'coords' has following structure:
%           X1     Y1     Z1
%           X2     Y2     Z2
%           .      .      .
%           .      .      .
%           .      .      .
%           Xn     Yn     Zn
%
%2. Topology file 'topo' should describe the hull of the distribution, thus
%   it should reffer only to its external boundary surface and not contain any
%   other connections between vertices.
%
%3. There is no limitation on the number of vertices that define each face.
%
%4. The order of the vertices in each line of file 'topo' describing the  
%   individual faces of the distribution, must be given counterclockwise,
%   so that the corresponding plane normal is outward pointing.
%
%5. Each line of file 'topo' must be referring to a different face.
%
%6. The origin of the local reference coordinate system is located at the 
%   computation point P (Xp=0, Yp=0, Zp=0), thus the coordinates of all 
%   vertices of the polyhedral distribution must refer to P.
%
%
% Literature
%
% Tsoulis D (2012) Analytical computation of the full gravity tensor of a 
% homogeneous arbitrarily shaped polyhedral source using line integrals, 
% Geophysics, 77, (2), F1-F11
%
% Tsoulis D and Gavriilidou G (2021) A computational review of the 
% line integral analytical formulation of the polyhedral gravity signal, 
% Geophysical Prospecting, 69, (8-9), 1745-1760
%
%
% (c) Authors: Dimitrios Tsoulis and Georgia Gavriilidou (April 2021)
%

%--------------------------------------------------------------------------
%-------------------- INITIAL DATA SETUP ----------------------------------


%Import the known vertex coordinates from file 'coords' 
fileID=fopen(coords,'r');
formatSpec = '%f';
sizeA = [3 Inf];
xyz = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
xyz=xyz';
nov=length(xyz);   %number of vertices 

x=xyz(:,1);
y=xyz(:,2);
z=xyz(:,3);


%Import and save the connectivity of the vertices from file 'topo'
topology = dlmread(topo);
nop=length(topology);  %Number of planes deriving from the topology file.



%The matrix noe (number of edges) contains the number of vertices defining
%each polyhedral planar face 
noe=zeros(1,nop);
size_topo=size(topology);
for i=1:nop
    loop=0;
    for j=1:size_topo(1,2)
        if topology(i,j)~=0
            loop=loop+1;
        else
            loop=loop+0;
        end
    end
    noe(1,i)=loop;   
end


%Base vectors e1,e2,e3 with origin at the computation point P
e=eye(3);


%Evaluation of segment vectors Gij for each polyhedral planar face
%according to Gij=g1*e1+g2*e2+g3*e3
Gij=zeros(size_topo(1,1),size_topo(1,2),3);
for i=1:nop
    for j=1:noe(1,i)
        if j==noe(1,i)
            Gij(i,j,1)=x(topology(i,1))-x(topology(i,j));
            Gij(i,j,2)=y(topology(i,1))-y(topology(i,j));
            Gij(i,j,3)=z(topology(i,1))-z(topology(i,j));
        else
            Gij(i,j,1)=x(topology(i,j+1))-x(topology(i,j));
            Gij(i,j,2)=y(topology(i,j+1))-y(topology(i,j));
            Gij(i,j,3)=z(topology(i,j+1))-z(topology(i,j));           
        end
    end   
end


%Computation of outward pointing planar unit normal vectors Np for each polyhedral planar face 
%from the equation: Np=Gi1xGi2/|Gi1xGi2|
Np=zeros(nop,3);
Npcosines=zeros(nop,3);
for i=1:nop
    A=[Gij(i,1,1), Gij(i,1,2), Gij(i,1,3)];  %Gi1
    B=[Gij(i,2,1), Gij(i,2,2), Gij(i,2,3)];  %Gi2
    [Ccross(1,1), Ccross(1,2), Ccross(1,3)]=cprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3)); %Gi1xGi2
    Cnorm=vlen(Ccross(1,1),Ccross(1,2),Ccross(1,3));  %|Gi1xGi2|
    Np(i,1)=Ccross(1,1)/Cnorm;
    Np(i,2)=Ccross(1,2)/Cnorm;
    Np(i,3)=Ccross(1,3)/Cnorm;
    Npcosines(i,1)=Np(i,1); 
    Npcosines(i,2)=Np(i,2); 
    Npcosines(i,3)=Np(i,3); 
    A=[0, 0, 0];
    B=[0, 0, 0];
    Ccross=[0, 0, 0];
    Cnorm=0;
end


%Computation of segment unit vectors npq from the equation: nij=GijxNp/|GijxNp|
npq=zeros(size_topo(1,1),size_topo(1,2),3);
npqcosines=zeros(size_topo(1,1),size_topo(1,2),3);
for i=1:nop
    for j=1:noe(i)
        A=[Gij(i,j,1), Gij(i,j,2), Gij(i,j,3)];  %Gij
        B=[Np(i,1), Np(i,2), Np(i,3)];  %Np
        [Ccross(1,1), Ccross(1,2), Ccross(1,3)]=cprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));  %GijxNp 
        Cnorm=vlen(Ccross(1,1),Ccross(1,2),Ccross(1,3));  %|GijxNp| 
        npq(i,j,1)=Ccross(1,1)/Cnorm;
        npq(i,j,2)=Ccross(1,2)/Cnorm;
        npq(i,j,3)=Ccross(1,3)/Cnorm;
        npqcosines(i,j,1)=npq(i,j,1); 
        npqcosines(i,j,2)=npq(i,j,2); 
        npqcosines(i,j,3)=npq(i,j,3); 
        A=[0, 0, 0];
        B=[0, 0, 0];
        Ccross=[0, 0, 0];
        Cnorm=0;
    end
end


%Computation of each vector's Np direction cosines cos(Np,ei) 
%from the equation: (Np)(ei)=|Np||ei|cos(Np,ei)
cosNp=zeros(nop,3);
for i=1:nop
   for j=1:3
       A=[Np(i,1), Np(i,2), Np(i,3)];
       B=[e(j,1), e(j,2), e(j,3)];
       [Cdot]=dprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));
       cosNp(i,j)=Cdot/(vlen(A(1,1),A(1,2),A(1,3)))*(vlen(B(1,1),B(1,2),B(1,3)));     
       A=[0, 0, 0];
       B=[0, 0, 0];
       Cdot=0;
   end
end
%----> cosNp=Npcosines


%Computation of each vector's npq direction cosines cos(npq,ei) 
%from the equation: (npq)(ei)=|npq||ei|cos(npq,ei)
cosnpq=zeros(size_topo(1,1),size_topo(1,2),3);
for i=1:nop
    for j=1:noe(i)
        for m=1:3
            A=[npq(i,j,1), npq(i,j,2), npq(i,j,3)];
            B=[e(m,1), e(m,2), e(m,3)];
            [Cdot]=dprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));
            cosnpq(i,j,m)=Cdot/(vlen(A(1,1),A(1,2),A(1,3)))*(vlen(B(1,1),B(1,2),B(1,3)));    
            A=[0, 0, 0];
            B=[0, 0, 0];
            Cdot=0;
        end
    end
end
%----> cosnpq=npqcosines



%------------------------- P' COORDINATES --------------------------------

%Parameter sp (= 0, -1, +1) from the sign of equation: (Ni)(-Gi1)
sp=zeros(nop,1);
for i=1:nop
    A=[Np(i,1), Np(i,2), Np(i,3)];
    B=[-x(topology(i,1)), -y(topology(i,1)), -z(topology(i,1)) ];
    [Cdot]=dprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));
    if Cdot>0
        sp(i,1)=-1;
    elseif Cdot<0
        sp(i,1)=1;
    elseif Cdot==0
        sp(i,1)=0;       
    end 
    A=[0, 0, 0];
    B=[0, 0, 0];
    Cdot=0;
end


%Calculation of the distances between P and the plane of each face, 
%as well as the intersections of each plane with the three axes x,y,z 
%using function planedist
%[hp,minusDA,minusDB,minusDC]=planedist(x,y,z,top1,top2,top3)
hp=zeros(nop,1);
mDA=zeros(nop,1);
mDB=zeros(nop,1);
mDC=zeros(nop,1);
for i=1:nop
    [hp(i,1), mDA(i,1), mDB(i,1), mDC(i,1)]=planedist(x,y,z,topology(i,1),topology(i,2),topology(i,3));    
end


%Final calculation of the coordinates of each projection P' on every polyhedral plane using 
%the equations: |xp|=|a*hp|, |yp|=|b*hp|, |zp|=|c*hp| where a,b,c are the direction cosines 
%of Np and hp the distance of P from the plane. 
%The coordinate signs derive from a further investigation on the intersections 
%of each plane with the three axes.
xx=zeros(nop,1);
yy=zeros(nop,1);
zz=zeros(nop,1);
for i=1:nop
    if mDA(i,1)>=0 && mDB(i,1)>=0 && mDC(i,1)>=0
        xx(i,1)=abs(Npcosines(i,1)*hp(i,1));
        yy(i,1)=abs(Npcosines(i,2)*hp(i,1));
        zz(i,1)=abs(Npcosines(i,3)*hp(i,1));
    
    elseif mDA(i,1)<0 && mDB(i,1)<0 && mDC(i,1)<0     
        if Npcosines(i,1)>0
            xx(i,1)=-(Npcosines(i,1)*hp(i,1));
        else
            xx(i,1)=(Npcosines(i,1)*hp(i,1));
        end        
        if Npcosines(i,2)>0
            yy(i,1)=-(Npcosines(i,2)*hp(i,1));
        else
            yy(i,1)=(Npcosines(i,2)*hp(i,1));
        end       
        if Npcosines(i,3)>0
            zz(i,1)=-(Npcosines(i,3)*hp(i,1));
        else
            zz(i,1)=(Npcosines(i,3)*hp(i,1));
        end
        
    elseif mDA(i,1)<0 && mDB(i,1)<0    
        if Npcosines(i,1)>0
            xx(i,1)=-(Npcosines(i,1)*hp(i,1));
        else
            xx(i,1)=(Npcosines(i,1)*hp(i,1));
        end        
        if Npcosines(i,2)>0
            yy(i,1)=-(Npcosines(i,2)*hp(i,1));
        else
            yy(i,1)=(Npcosines(i,2)*hp(i,1));
        end
        zz(i,1)=abs(Npcosines(i,3)*hp(i,1));
       
    elseif mDB(i,1)<0 && mDC(i,1)<0     
        xx(i,1)=abs(Npcosines(i,1)*hp(i,1));             
        if Npcosines(i,2)>0
            yy(i,1)=-(Npcosines(i,2)*hp(i,1));
        else
            yy(i,1)=(Npcosines(i,2)*hp(i,1));
        end       
        if Npcosines(i,3)>0
            zz(i,1)=-(Npcosines(i,3)*hp(i,1));
        else
            zz(i,1)=(Npcosines(i,3)*hp(i,1));
        end
    
    elseif mDA(i,1)<0 && mDC(i,1)<0     
        if Npcosines(i,1)>0
            xx(i,1)=-(Npcosines(i,1)*hp(i,1));
        else
            xx(i,1)=(Npcosines(i,1)*hp(i,1));
        end        
        yy(i,1)=abs(Npcosines(i,2)*hp(i,1));      
        if Npcosines(i,3)>0
            zz(i,1)=-(Npcosines(i,3)*hp(i,1));
        else
            zz(i,1)=(Npcosines(i,3)*hp(i,1));
        end
        
    elseif mDA(i,1)<0    
        if Npcosines(i,1)>0
            xx(i,1)=-(Npcosines(i,1)*hp(i,1));
        else
            xx(i,1)=(Npcosines(i,1)*hp(i,1));
        end        
        yy(i,1)=abs(Npcosines(i,2)*hp(i,1));
        zz(i,1)=abs(Npcosines(i,3)*hp(i,1));
        
    elseif mDB(i,1)<0    
        xx(i,1)=abs(Npcosines(i,1)*hp(i,1));
        if Npcosines(i,2)>0
            yy(i,1)=-(Npcosines(i,2)*hp(i,1));
        else
            yy(i,1)=(Npcosines(i,2)*hp(i,1));
        end       
        zz(i,1)=abs(Npcosines(i,3)*hp(i,1));
        
    elseif mDC(i,1)<0     
        xx(i,1)=abs(Npcosines(i,1)*hp(i,1));
        yy(i,1)=abs(Npcosines(i,2)*hp(i,1));
        if Npcosines(i,3)>0
            zz(i,1)=-(Npcosines(i,3)*hp(i,1));
        else
            zz(i,1)=(Npcosines(i,3)*hp(i,1));
        end  
    end
end
%-------->The coordinates of points P' are saved into 
%matrices xx,yy,zz.



%------------------------ P'' COORDINATES ---------------------------------

%Parameter spq (= 0, -1, +1) from the sign of the equation: nij(xp'-xij)
spq=zeros(size_topo(1,1),size_topo(1,2));
for i=1:nop
    for j=1:noe(i)
        A=[npq(i,j,1), npq(i,j,2), npq(i,j,3)];
        B=[(xx(i,1)-x(topology(i,j))), (yy(i,1)-y(topology(i,j))), (zz(i,1)-z(topology(i,j)))];
        [Cdot]=dprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));
        if abs(Cdot)<(1e-10)
            %practically <0
            spq(i,j)=0; 
        elseif Cdot>0
            spq(i,j)=-1;
        elseif Cdot<0
            spq(i,j)=1;
        elseif Cdot==0
            spq(i,j)=0;       
        end 
        A=[0, 0, 0];
        B=[0, 0, 0];
        Cdot=0;
    end
end


%Computation of the coordinates of every point P'' using function Ppq
%[xp,yp,zp]=Ppq(x1,y1,z1,x2,y2,z2,x0,y0,z0)
xxx=zeros(size_topo(1,1),size_topo(1,2));
yyy=zeros(size_topo(1,1),size_topo(1,2));
zzz=zeros(size_topo(1,1),size_topo(1,2));
hpq=zeros(size_topo(1,1),size_topo(1,2));
for i=1:nop
   for j=1:noe(i)
       if j==noe(i)
           if spq(i,j)==0
               hpq(i,j)=0;
               xxx(i,j)=xx(i,1);
               yyy(i,j)=yy(i,1);
               zzz(i,j)=zz(i,1);
           else               
               [xxx(i,j), yyy(i,j), zzz(i,j)]=Ppq(x(topology(i,j)),y(topology(i,j)),z(topology(i,j)),x(topology(i,1)),y(topology(i,1)),z(topology(i,1)),xx(i,1),yy(i,1),zz(i,1));
               A=[(xxx(i,j)-xx(i,1)), (yyy(i,j)-yy(i,1)), (zzz(i,j)-zz(i,1))];
               hpq(i,j)=vlen(A(1,1),A(1,2),A(1,3));  
               A=[0, 0, 0];
           end
       else
           if spq(i,j)==0
               hpq(i,j)=0;
               xxx(i,j)=xx(i,1);
               yyy(i,j)=yy(i,1);
               zzz(i,j)=zz(i,1);
           else
               [xxx(i,j), yyy(i,j), zzz(i,j)]=Ppq(x(topology(i,j)),y(topology(i,j)),z(topology(i,j)),x(topology(i,j+1)),y(topology(i,j+1)),z(topology(i,j+1)),xx(i,1),yy(i,1),zz(i,1));
               A=[(xxx(i,j)-xx(i,1)), (yyy(i,j)-yy(i,1)), (zzz(i,j)-zz(i,1))];
               hpq(i,j)=vlen(A(1,1),A(1,2),A(1,3));
               A=[0, 0, 0];
           end    
       end      
   end
end
%--------->The coordinates of points P'' are saved into
%matrices xxx,yyy,zzz.



%---------------------------- V & Vxi -------------------------------------

V=0;
Vx=0;
Vy=0;
Vz=0;

%The first loop runs through every face
for i=1:nop
    
    %sum1: first common part of the equations for V and Vi
    %sum2: second common part of the equations for V and Vi
    sum1=0;
    sum2=0;
    singarea=0; %variable for controlling whether P' is located inside the corresponding face
    singsegm=0; %variable for controlling whether P' is situated on a segment  
    singvert=0; %variable for controlling whether P' coincides with a vertex
    area=0;
    edge=0;
    vertex=0;
    theta=0;
    
   %The second loop runs through every segment of each face
   for j=1:noe(1,i)
       thisvert=0;
       
       %Segment's Gij length
       A=[Gij(i,j,1), Gij(i,j,2), Gij(i,j,3)]; 
       Gijnorm=0;
       Gijnorm=vlen(A(1,1),A(1,2),A(1,3));
       A=[0, 0, 0];
       %Check if P' is located inside the corresponding face 
       %and calculation of the correction term SING(Ap)=-2*pi*h(p) is due
       if spq(i,j)==1
           singarea=singarea+1;   
       end
       
       %Check if P' lies on the direction of a segment. 
       %This possibility involves two cases of P' situated onto 
       %the corresponding segment or at a vertex.
       if spq(i,j)==0
           r1norm=0;
           r2norm=0;
           r1=[0, 0, 0];
           r2=[0, 0, 0];
           %Computation of vectors r1 and r2.
           if j==noe(i)
               r1x=xx(i,1)-x(topology(i,j));
               r1y=yy(i,1)-y(topology(i,j));
               r1z=zz(i,1)-z(topology(i,j));
               r2x=xx(i,1)-x(topology(i,1));
               r2y=yy(i,1)-y(topology(i,1));
               r2z=zz(i,1)-z(topology(i,1)); 
           else
               r1x=xx(i,1)-x(topology(i,j));
               r1y=yy(i,1)-y(topology(i,j));
               r1z=zz(i,1)-z(topology(i,j));
               r2x=xx(i,1)-x(topology(i,j+1));
               r2y=yy(i,1)-y(topology(i,j+1));
               r2z=zz(i,1)-z(topology(i,j+1));          
           end
           r1=[r1x, r1y, r1z];
           r2=[r2x, r2y, r2z];
           r1norm=vlen(r1(1,1),r1(1,2),r1(1,3));
           r2norm=vlen(r2(1,1),r2(1,2),r2(1,3));
                  
           
           %If the length of r1 or r2 is equal to zero, then P' 
           %coincides with a vertex and thus with P''. 
           %In this case, the correction SING(Ap)=-theta*h(p) is calculated,
           %where theta=arccos{(G2(-G1))/(|G2||-G1|)}.
           if r1norm==0 
               singvert=singvert+1; 
               thisvert=thisvert+1;   
               if j==1
                   A=[-Gij(i,noe(i),1), -Gij(i,noe(i),2), -Gij(i,noe(i),3)]; 
                   B=[Gij(i,1,1), Gij(i,1,2), Gij(i,1,3)]; 
                   G1norm=vlen(A(1,1),A(1,2),A(1,3));
                   G2norm=vlen(B(1,1),B(1,2),B(1,3));
                   [G1dotG2]=dprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));
                   if G1dotG2==0
                       theta=pi/2;
                   else
                       theta=acos(G1dotG2/(G1norm*G2norm));   
                   end
                   A=[0, 0, 0];
                   B=[0, 0, 0];
                   G1norm=0;
                   G2norm=0;
                   G1dotG2=0;
               else
                   A=[-Gij(i,j-1,1), -Gij(i,j-1,2), -Gij(i,j-1,3)]; 
                   B=[Gij(i,j,1), Gij(i,j,2), Gij(i,j,3)]; 
                   G1norm=vlen(A(1,1),A(1,2),A(1,3));
                   G2norm=vlen(B(1,1),B(1,2),B(1,3));
                   [G1dotG2]=dprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));
                   if G1dotG2==0
                       theta=pi/2;
                   else
                       theta=acos(G1dotG2/(G1norm*G2norm));   
                   end 
                   A=[0, 0, 0];
                   B=[0, 0, 0];
                   G1norm=0;
                   G2norm=0;
                   G1dotG2=0;
               end
           end
                                            
           if r2norm==0 
               singvert=singvert+1; 
               thisvert=thisvert+1;
               if j==noe(i)
                   A=[-Gij(i,j,1), -Gij(i,j,2), -Gij(i,j,3)]; 
                   B=[Gij(i,1,1), Gij(i,1,2), Gij(i,1,3)]; 
                   G1norm=vlen(A(1,1),A(1,2),A(1,3));
                   G2norm=vlen(B(1,1),B(1,2),B(1,3));
                   [G1dotG2]=dprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));
                   if G1dotG2==0
                       theta=pi/2;
                   else
                       theta=acos(G1dotG2/(G1norm*G2norm));   
                   end  
                   A=[0, 0, 0];
                   B=[0, 0, 0];
                   G1norm=0;
                   G2norm=0;
                   G1dotG2=0;
               else
                   A=[-Gij(i,j,1), -Gij(i,j,2), -Gij(i,j,3)]; 
                   B=[Gij(i,j+1,1), Gij(i,j+1,2), Gij(i,j+1,3)]; 
                   G1norm=vlen(A(1,1),A(1,2),A(1,3));
                   G2norm=vlen(B(1,1),B(1,2),B(1,3));
                   [G1dotG2]=dprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));
                   if G1dotG2==0
                       theta=pi/2;
                   else
                       theta=acos(G1dotG2/(G1norm*G2norm));   
                   end  
                   A=[0, 0, 0];
                   B=[0, 0, 0];
                   G1norm=0;
                   G2norm=0;
                   G1dotG2=0;
               end     
           end
       
           
           %If the length of both r1 and r2 are less than the segment's length, 
           %then P' is situated on the corresponding segment. In this 
           %case, the correction SING(Ap)=-pi*h(p) is calculated.
           if r1norm<Gijnorm && r2norm<Gijnorm
               singsegm=singsegm+1;
           end 
       end
       
       
       %Definition of s1pq, s2pq, l1pq, l2pq vectors. 
       %First, their lengths are calculated and then with a further 
       %control, the correct sign of each distance is given.
       s1abs=0;
       s2abs=0;
       l1abs=0;
       l2abs=0;
       s1=0;
       s2=0;
       l1=0;
       l2=0;
       if j==noe(i)
           A=[xxx(i,j)-x(topology(i,j),1), yyy(i,j)-y(topology(i,j),1), zzz(i,j)-z(topology(i,j),1)];
           B=[xxx(i,j)-x(topology(i,1),1), yyy(i,j)-y(topology(i,1),1), zzz(i,j)-z(topology(i,1),1)];
           s1abs=vlen(A(1,1),A(1,2),A(1,3));
           s2abs=vlen(B(1,1),B(1,2),B(1,3));
           A=[0, 0, 0];
           B=[0, 0, 0];
           A=[x(topology(i,j),1), y(topology(i,j),1), z(topology(i,j),1)];
           B=[x(topology(i,1),1), y(topology(i,1),1), z(topology(i,1),1)];
           l1abs=vlen(A(1,1),A(1,2),A(1,3));
           l2abs=vlen(B(1,1),B(1,2),B(1,3));  
           A=[0, 0, 0];
           B=[0, 0, 0];
       else
           A=[xxx(i,j)-x(topology(i,j),1), yyy(i,j)-y(topology(i,j),1), zzz(i,j)-z(topology(i,j),1)];
           B=[xxx(i,j)-x(topology(i,j+1),1), yyy(i,j)-y(topology(i,j+1),1), zzz(i,j)-z(topology(i,j+1),1)];
           s1abs=vlen(A(1,1),A(1,2),A(1,3));
           s2abs=vlen(B(1,1),B(1,2),B(1,3));
           A=[0, 0, 0];
           B=[0, 0, 0];
           A=[x(topology(i,j),1), y(topology(i,j),1), z(topology(i,j),1)];
           B=[x(topology(i,j+1),1), y(topology(i,j+1),1), z(topology(i,j+1),1)];
           l1abs=vlen(A(1,1),A(1,2),A(1,3));
           l2abs=vlen(B(1,1),B(1,2),B(1,3));
           A=[0, 0, 0];
           B=[0, 0, 0];
       end
       
       if abs(s1abs-l1abs)<1e-10 && abs(s2abs-l2abs)<1e-10
           %|s1-l1|=0 & |s2-l2|=0 
           
           if s1abs<Gijnorm && s2abs<Gijnorm 
               s1=-s1abs;
               s2=s2abs;
               l1=-l1abs;
               l2=l2abs;
           elseif s1abs<s2abs  
               s1=s1abs;
               s2=s2abs;
               l1=l1abs;
               l2=l2abs;             
           elseif s2abs<s1abs  
               s1=-s1abs;
               s2=-s2abs;
               l1=-l1abs;
               l2=-l2abs;       
           end
          
       else
           if s1abs<Gijnorm && s2abs<Gijnorm
               s1=-s1abs;
               s2=s2abs;
               l1=l1abs;
               l2=l2abs; 
           else
               if s1abs<s2abs 
                   s1=s1abs;
                   s2=s2abs;
                   l1=l1abs;
                   l2=l2abs; 
               elseif s2abs<s1abs 
                   s1=-s1abs;
                   s2=-s2abs;
                   l1=l1abs;
                   l2=l2abs;
               end
           end
       end
      
       
       %Computation of LNpq and ANpq.
       AN=0;
       LN=0;
       if hp(i,1)==0 
           if thisvert>0 
               [LN]=0;
               [AN]=0;
           elseif abs(s1+s2)<1e-10 && abs(l1+l2)<1e-10
               [LN]=0;
               [AN]=0;
           else
               [LN]=LNpq(s2,l2,s1,l1);
               [AN]=0;
           end         
       else 

           if thisvert>0
               [LN]=0;
           else
               [LN]=LNpq(s2,l2,s1,l1);
           end

           if hpq(i,j)==0 
               [AN]=0; 
           else
               [AN]=ANpq(hp(i,1),hpq(i,j),s2,l2,s1,l1);
           end
       end
       
       sum1=sum1+(spq(i,j)*hpq(i,j)*LN);
       sum2=sum2+(spq(i,j)*AN);  
       
   end
   
   
   %The investigations for this plane's segments are finished, so the second loop is concluded. 
   %Before the calculations examine the next polyhedral plane, the exact value of 
   %the correction SING(AP) is computed by checking the variables 
   %singarea, singsegm, singvert. 
   
   if singarea==noe(1,i) 
       area=-(2*pi*hp(i,1));
   else
       area=0;
   end
   if singsegm>0 
       edge=-(pi*hp(i,1));
   else
       edge=0;
   end
   if singvert>0 
       vertex=-(theta*hp(i,1));
   else
       vertex=0;
   end 
   
   
   V=V+(sp(i,1)*hp(i,1)*(sum1+(hp(i,1)*sum2)+area+edge+vertex));
   Vx=Vx+(cosNp(i,1)*(sum1+(hp(i,1)*sum2)+area+edge+vertex));
   Vy=Vy+(cosNp(i,2)*(sum1+(hp(i,1)*sum2)+area+edge+vertex));
   Vz=Vz+(cosNp(i,3)*(sum1+(hp(i,1)*sum2)+area+edge+vertex));
   
end


%Final expression for the potential V and its first derivatives Vi
V=(G*p*V)/2;
Vx=abs(G*p*Vx);
Vy=abs(G*p*Vy);
Vz=abs(G*p*Vz);


%--------------------------- Vxixj ----------------------------------------

Vij=zeros(3);

for k=1:3
    for m=1:3
        for i=1:nop
           
            %Variables sum1 and sum2 are defined as the two parts of 
            %the equation for the second derivatives. 
            %Variables singarea, singsegm, singvert are defined as well 
            %for enabling the investigations about the location of P'. 
            sum1=0;
            sum2=0;
            singarea=0;
            singsegm=0;
            singvert=0;
            area=0;
            edge=0;
            vertex=0;
            theta=0;
            
            for j=1:noe(1,i)
                thisvert=0;
                
                %Vector's Gij length
                A=[Gij(i,j,1), Gij(i,j,2), Gij(i,j,3)]; 
                Gijnorm=0;
                Gijnorm=vlen(A(1,1),A(1,2),A(1,3));
                A=[0, 0, 0];
                
                
                %The case where P' lies inside a face
                if spq(i,j)==1
                    singarea=singarea+1; 
                end
                
                
                %The case where P' lies on the direction of the segment 
                %is investigated. 
                %In this case P' may be located inside the segment or 
                %exactly on one vertex.
                if spq(i,j)==0
                    r1norm=0;
                    r2norm=0;
                    r1=[0, 0, 0];
                    r2=[0, 0, 0];
                    %Computation of vectors r1 and r2
                    if j==noe(1,i)
                        r1x=xx(i,1)-x(topology(i,j));
                        r1y=yy(i,1)-y(topology(i,j));
                        r1z=zz(i,1)-z(topology(i,j));
                        r2x=xx(i,1)-x(topology(i,1));
                        r2y=yy(i,1)-y(topology(i,1));
                        r2z=zz(i,1)-z(topology(i,1));
                    else
                        r1x=xx(i,1)-x(topology(i,j));
                        r1y=yy(i,1)-y(topology(i,j));
                        r1z=zz(i,1)-z(topology(i,j));
                        r2x=xx(i,1)-x(topology(i,j+1));
                        r2y=yy(i,1)-y(topology(i,j+1));
                        r2z=zz(i,1)-z(topology(i,j+1)); 
                    end
                    r1=[r1x, r1y, r1z];
                    r2=[r2x, r2y, r2z];
                    r1norm=vlen(r1(1,1),r1(1,2),r1(1,3));
                    r2norm=vlen(r2(1,1),r2(1,2),r2(1,3));
                    
                    
                    %If the length of r1 or r2 is equal to zero, then P' 
                    %coincides with a vertex and thus with P''. 
                    %In this case, the correction SING(Bp) must be
                    %calculated.
                    if r1norm==0 
                        singvert=singvert+1;
                        thisvert=thisvert+1; 
                        
                        if j==1
                            A=[-Gij(i,noe(1,i),1), -Gij(i,noe(1,i),2), -Gij(i,noe(1,i),3)]; 
                            B=[Gij(i,1,1), Gij(i,1,2), Gij(i,1,3)]; 
                            G1norm=vlen(A(1,1),A(1,2),A(1,3));
                            G2norm=vlen(B(1,1),B(1,2),B(1,3));
                            [G1dotG2]=dprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));
                            if G1dotG2==0
                                theta=pi/2;
                            else
                                theta=acos(G1dotG2/(G1norm*G2norm));
                            end
                            A=[0, 0, 0];
                            B=[0, 0, 0];
                            G1norm=0;
                            G2norm=0;
                            G1dotG2=0;
                        else
                            A=[-Gij(i,j-1,1), -Gij(i,j-1,2), -Gij(i,j-1,3)]; 
                            B=[Gij(i,j,1), Gij(i,j,2), Gij(i,j,3)]; 
                            G1norm=vlen(A(1,1),A(1,2),A(1,3));
                            G2norm=vlen(B(1,1),B(1,2),B(1,3));
                            [G1dotG2]=dprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));
                            if G1dotG2==0
                                theta=pi/2;
                            else
                                theta=acos(G1dotG2/(G1norm*G2norm));
                            end  
                            A=[0, 0, 0];
                            B=[0, 0, 0];
                            G1norm=0;
                            G2norm=0;
                            G1dotG2=0;
                        end
                    end
                    
                    if r2norm==0 
                        singvert=singvert+1;
                        thisvert=thisvert+1; 
                        
                        if j==noe(1,i)
                            A=[-Gij(i,j,1), -Gij(i,j,2), -Gij(i,j,3)]; 
                            B=[Gij(i,1,1), Gij(i,1,2), Gij(i,1,3)]; 
                            G1norm=vlen(A(1,1),A(1,2),A(1,3));
                            G2norm=vlen(B(1,1),B(1,2),B(1,3));
                            [G1dotG2]=dprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));
                            if G1dotG2==0
                                theta=pi/2;
                            else
                                theta=acos(G1dotG2/(G1norm*G2norm));
                            end
                            A=[0, 0, 0];
                            B=[0, 0, 0];
                            G1norm=0;
                            G2norm=0;
                            G1dotG2=0;
                        else
                            A=[-Gij(i,j,1), -Gij(i,j,2), -Gij(i,j,3)]; 
                            B=[Gij(i,j+1,1), Gij(i,j+1,2), Gij(i,j+1,3)]; 
                            G1norm=vlen(A(1,1),A(1,2),A(1,3));
                            G2norm=vlen(B(1,1),B(1,2),B(1,3));
                            [G1dotG2]=dprod(A(1,1),A(1,2),A(1,3),B(1,1),B(1,2),B(1,3));
                            if G1dotG2==0
                                theta=pi/2;
                            else
                                theta=acos(G1dotG2/(G1norm*G2norm));
                            end  
                            A=[0, 0, 0];
                            B=[0, 0, 0];
                            G1norm=0;
                            G2norm=0;
                            G1dotG2=0;
                        end
                    end   
                   
                    
                    %If the length of both r1 and r2 are less than the segment's length, 
                    %then P' is situated on the corresponding segment and the 
                    %correction SING(Bp) is calculated.
                    if r1norm<Gijnorm && r2norm<Gijnorm
                        singsegm=singsegm+1;
                    end
                    
                end
                
                                
                %The vectors s1pq, s2pq, l1pq, l2pq are defined. 
                %First their lengths are calculated and then 
                %the correct sign of each distance is given.  
                s1abs=0;
                s2abs=0;
                l1abs=0;
                l2abs=0;
                s1=0;
                s2=0;
                l1=0;
                l2=0;
                if j==noe(1,i)
                    A=[xxx(i,j)-x(topology(i,j)), yyy(i,j)-y(topology(i,j)), zzz(i,j)-z(topology(i,j))];
                    B=[xxx(i,j)-x(topology(i,1)), yyy(i,j)-y(topology(i,1)), zzz(i,j)-z(topology(i,1))];
                    s1abs=vlen(A(1,1),A(1,2),A(1,3));
                    s2abs=vlen(B(1,1),B(1,2),B(1,3));
                    A=[0, 0, 0];
                    B=[0, 0, 0];
                    A=[x(topology(i,j)), y(topology(i,j)), z(topology(i,j))];
                    B=[x(topology(i,1)), y(topology(i,1)), z(topology(i,1))];
                    l1abs=vlen(A(1,1),A(1,2),A(1,3));
                    l2abs=vlen(B(1,1),B(1,2),B(1,3));    
                    A=[0, 0, 0];
                    B=[0, 0, 0];
                else
                    A=[xxx(i,j)-x(topology(i,j)), yyy(i,j)-y(topology(i,j)), zzz(i,j)-z(topology(i,j))];
                    B=[xxx(i,j)-x(topology(i,j+1)), yyy(i,j)-y(topology(i,j+1)), zzz(i,j)-z(topology(i,j+1))];
                    s1abs=vlen(A(1,1),A(1,2),A(1,3));
                    s2abs=vlen(B(1,1),B(1,2),B(1,3));
                    A=[0, 0, 0];
                    B=[0, 0, 0];
                    A=[x(topology(i,j)), y(topology(i,j)), z(topology(i,j))];
                    B=[x(topology(i,j+1)), y(topology(i,j+1)), z(topology(i,j+1))];
                    l1abs=vlen(A(1,1),A(1,2),A(1,3));
                    l2abs=vlen(B(1,1),B(1,2),B(1,3));
                    A=[0, 0, 0];
                    B=[0, 0, 0];
                end
                
                if abs(s1abs-l1abs)<1e-10 && abs(s2abs-l2abs)<1e-10
                    %|s1-l1|=0 & |s2-l2|=0 
                    if s1abs<Gijnorm && s2abs<Gijnorm 
                        s1=-s1abs;
                        s2=s2abs;
                        l1=-l1abs;
                        l2=l2abs;
                    elseif s1abs<s2abs 
                        s1=s1abs;
                        s2=s2abs;
                        l1=l1abs;
                        l2=l2abs;             
                    elseif s2abs<s1abs  
                        s1=-s1abs;
                        s2=-s2abs;
                        l1=-l1abs;
                        l2=-l2abs;                        
                    end
                
                else
                    if s1abs<Gijnorm && s2abs<Gijnorm
                        s1=-s1abs;
                        s2=s2abs;
                        l1=l1abs;
                        l2=l2abs; 
                    else
                        if s1abs<s2abs 
                            s1=s1abs;
                            s2=s2abs;
                            l1=l1abs;
                            l2=l2abs; 
                        elseif s2abs<s1abs 
                            s1=-s1abs;
                            s2=-s2abs;
                            l1=l1abs;
                            l2=l2abs;
                        end
                    end
                end
        
                %Computation of LNpq and ANpq.
                AN=0;
                LN=0;
                if hp(i,1)==0 
                   
                    if thisvert>0 
                        [LN]=0;
                        [AN]=0;
                    elseif abs(s1+s2)<1e-10 && abs(l1+l2)<1e-10
                        [LN]=0;
                        [AN]=0;
                    else
                        [LN]=LNpq(s2,l2,s1,l1);
                        [AN]=0;
                    end  
                    
                else 
                    
                    if thisvert>0
                        [LN]=0;
                    else
                        [LN]=LNpq(s2,l2,s1,l1);
                    end
                    
                    if hpq(i,j)==0 
                        [AN]=0; 
                    else
                        [AN]=ANpq(hp(i,1),hpq(i,j),s2,l2,s1,l1);
                    end
                end
                
                sum1=sum1+(cosnpq(i,j,m)*LN);
                sum2=sum2+(spq(i,j)*AN);
                
            end
            
            
            %Calculation of correction SING(Bp) by checking the values 
            %of variables singarea, singsegm, singvert.
            if singarea==noe(1,i) 
                area=-(2*pi);
            else
                area=0;
            end
            if singsegm>0 
                edge=-(pi);
            else
                edge=0;
            end
            if singvert>0
                vertex=-(theta); 
            else
                vertex=0;
            end
           
            
            %Second order derivatives 
            Vij(k,m)=Vij(k,m)+((cosNp(i,k)*(sum1+(sp(i,1)*cosNp(i,m))*(sum2+area+vertex+edge))));
            
        end        
    end   
end

%Final expression 
for i=1:3
    for j=1:3
        Vij(i,j)=G*p*Vij(i,j);
    end
end


%Print the results into a text file
fileID = fopen('results.txt','w');
fprintf(fileID,'V = %.15d \r\n',V);
fprintf(fileID,'Vx = %.15d \r\n',Vx);
fprintf(fileID,'Vy = %.15d \r\n',Vy);
fprintf(fileID,'Vz = %.15d \r\n',Vz);
fprintf(fileID,'Vxx = %.15d \r\n',Vij(1,1));
fprintf(fileID,'Vyy = %.15d \r\n',Vij(2,2));
fprintf(fileID,'Vzz = %.15d \r\n',Vij(3,3));
fprintf(fileID,'Vxy = %.15d \r\n',Vij(1,2));
fprintf(fileID,'Vxz = %.15d \r\n',Vij(1,3));
fprintf(fileID,'Vyz = %.15d \r\n',Vij(2,3));
fclose(fileID);

%Print the results on the command window
fprintf('V = %.15d \n',V)
fprintf('Vx = %.15d \n',Vx)
fprintf('Vy = %.15d \n',Vy)
fprintf('Vz = %.15d \n',Vz)
fprintf('Vxx = %.15d \n',Vij(1,1))
fprintf('Vyy = %.15d \n',Vij(2,2))
fprintf('Vzz = %.15d \n',Vij(3,3))
fprintf('Vxy = %.15d \n',Vij(1,2))
fprintf('Vxz = %.15d \n',Vij(1,3))
fprintf('Vyz = %.15d \n',Vij(2,3))


%--------------------------------------------------------------------------
%------------------------- FUNCTIONS --------------------------------------


% A. Function for computing LNpq
%    Input data: the distances s2, s1, l2, l1
%    Output data: one value for quantity LNpq
%    Function form: [a]=LNpq(s2,l2,s1,l1)

function [a]=LNpq(s2,l2,s1,l1)
a=log((s2+l2)/(s1+l1));
end


% B. Function for computing ANpq 
%    Input data: the distances s2, s1, l2, l1 as well as the distances hp
%    and hpq
%    Output data: one value for quantity ANpq
%    Function form: [a]=ANpq(hp,hpq,s2,l2,s1,l1)

function [a]=ANpq(hp,hpq,s2,l2,s1,l1)
a=(atan((hp*s2)/(hpq*l2)))-(atan((hp*s1)/(hpq*l1)));
end


% C. Function for computing the distance of a point from a plane.
%    Input data: the coordinates of the considered point and 
%    of three random points belonging to the plane.
%    Output data: the distance of the point from the plane and the intersections 
%    of the plane with the three axes.
%    Function form: [hp,minusDA,minusDB,minusDC]=planedist(x,y,z,top1,top2,top3)
%    --->Here the variables top1, top2, top3 are providing the number of the 
%        vertex given from the coordinate file, while x,y,z are the coordinates of 
%        the considered point. 

function [hp,minusDA,minusDB,minusDC]=planedist(x,y,z,top1,top2,top3)
a=((y(top2)-y(top1))*(z(top3)-z(top1)))-((y(top3)-y(top1))*(z(top2)-z(top1)));
b=((z(top2)-z(top1))*(x(top3)-x(top1)))-((z(top3)-z(top1))*(x(top2)-x(top1)));
c=((x(top2)-x(top1))*(y(top3)-y(top1)))-((x(top3)-x(top1))*(y(top2)-y(top1)));
d=-(a*x(top1))-(b*y(top1))-(c*z(top1));
hp=abs(d/sqrt((a^2)+(b^2)+(c^2)));
%Controlling special cases
if z(top1)==z(top2) && z(top2)==z(top3)
    minusDA=0;
    minusDB=0;
    minusDC=-d/c;
elseif y(top1)==y(top2) && y(top2)==y(top3)
    minusDA=0;
    minusDB=-d/b;
    minusDC=0;
elseif x(top1)==x(top2) && x(top2)==x(top3)
    minusDA=-d/a;
    minusDB=0;
    minusDC=0;
else
    if a == 0
        minusDA=0;
    else
        minusDA=-d/a;
    end
    if b == 0
        minusDB=0;
    else
        minusDB=-d/b;
    end
    if c == 0
        minusDC=0;
    else
        minusDC=-d/c;
    end
end
end


% D. Function for computing the coordinates of point P''.
%    Input data: the coordinates of the two vertices constructing the considered 
%    segment and the coordinates of the corresponding point P'.
%    Output data: the three coordinates of P''.
%    Function form: [xp,yp,zp]=Ppq(x1,y1,z1,x2,y2,z2,x0,y0,z0)
%    --->Here the coordinates x0,y0,z0 are referring to P' and 
%        x1,y1,z1,x2,y2,z2 to the two vertices defining the segment.

function [xp,yp,zp]=Ppq(x1,y1,z1,x2,y2,z2,x0,y0,z0)

%Vector G
a1=x2-x1;
b1=y2-y1;
c1=z2-z1;
A1=[a1, b1, c1];  %A1=G

%Vector y
k1=x1-x0;
k2=y1-y0;
k3=z1-z0;
A0=[k1, k2, k3];  %A0=y

%First equation: d1=X'G where X' is the position vector of point P'
d1=a1*x0+b1*y0+c1*z0;

%Second equation: d2=X'(Gxy) 
[A2(1,1), A2(1,2), A2(1,3)]=cprod(A1(1,1),A1(1,2),A1(1,3),A0(1,1),A0(1,2),A0(1,3)); %A2=Gxy
a2=A2(1,1);
b2=A2(1,2);
c2=A2(1,3);
d2=a2*x0+b2*y0+c2*z0;

%Third equation: d3=X1((Gxy)xG)
[A3(1,1), A3(1,2), A3(1,3)]=cprod(A2(1,1),A2(1,2),A2(1,3),A1(1,1),A1(1,2),A1(1,3)); %A3=(Gxy)xG
a3=A3(1,1);
b3=A3(1,2);
c3=A3(1,3);
d3=a3*x1+b3*y1+c3*z1;

%Solving the three equation system using determinants
det=(a1*b2*c3)-(a1*b3*c2)-(b1*a2*c3)+(b1*a3*c2)+(c1*a2*b3)-(c1*a3*b2);
dx=(d1*b2*c3)-(d1*b3*c2)-(b1*d2*c3)+(b1*d3*c2)+(c1*d2*b3)-(c1*d3*b2);
dy=(a1*d2*c3)-(a1*d3*c2)-(d1*a2*c3)+(d1*a3*c2)+(c1*a2*d3)-(c1*a3*d2);
dz=(a1*b2*d3)-(a1*b3*d2)-(b1*a2*d3)+(b1*a3*d2)+(d1*a2*b3)-(d1*a3*b2);
xp=dx/det;
yp=dy/det;
zp=dz/det;

end


% E. Function for computing the Euclidean norm of a vector.
%    Input data: the components x1, y1, z1 of a vector.
%    Output data: the calculated norm. 
%    Function form: [v]=vlen(x1,y1,z1)

function [v]=vlen(x1,y1,z1)
v=sqrt((x1*x1)+(y1*y1)+(z1*z1));
end

% F. Function for the computation of the cross product of two vectors.
%    Input data: the components a1,a2,a3 and b1,b2,b3 of the two vectors.
%    Output data: the cross product vector components. 
%    Function form: [c1, c2, c3]=cprod(a1,a2,a3,b1,b2,b3)

function [c1, c2, c3]=cprod(a1,a2,a3,b1,b2,b3)
c1=(a2*b3)-(a3*b2);
c2=(a3*b1)-(a1*b3);
c3=(a1*b2)-(a2*b1);
end


% G. Function for the computation of the dot product of two vectors.
%    Input data: the components a1,a2,a3 and b1,b2,b3 of the two vectors.
%    Output data: the dot product of the two considered vectors.
%    Function form: [c]=dprod(a1,a2,a3,b1,b2,b3)

function [c]=dprod(a1,a2,a3,b1,b2,b3)
c=(a1*b1)+(a2*b2)+(a3*b3);
end



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
end

