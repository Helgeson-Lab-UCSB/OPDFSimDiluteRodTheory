function [Y,a2,a4,Ytot,phivect,thetavect] = FokkerPlanckDiluteRodSolver(nangle,k,Dr,ARG,Yinitcond,tstart,store,tstop,dt)
%%Outputs orientation distribution in Y (also called N or f - see  
%Leal, Hinch 1976 or Bird et al DPL Vol 2)
%Adapted from R. Bay Thesis 1991, Crank Nicholson implemented according
%to Ferec et al 2008
%   
%   Patrick Corona (corona@ucsb.edu)
%
%   OUTPUTS
%   Y: Final Orientation Distribution Function
%   a2: Second moment tensor of final orientation distribution function
%       order of moments is 11, 12, 13, 22, 23, 33
%   a4: Fourth moment tensor of final orientation distribution function
%       order of moments is 1111, 2222, 3333, 1122, 1133, 2233, 1222, 2333, 
%       1333, 1112, 1113, 2223, 1233, 1223, 1123
%   Ytot: integral of the orientaiton distribution function, should be
%   1 (the normalization value of the OPDF), provides a proxy
%   for if the OPDF is finely discritized enough
%   phivect: phi values for Y
%   thetavect: theta values for Y
%
%   INPUTS
%   nangle: number of angles in the theta and phi directions, implemented 
%       as a square grid currently
%   k: velocity gradient tensor as a 3x3 matrix
%   Dr: rotational diffusivity
%   ARG: effective aspect ratio of the rods
%   Yinitcond: initial condition for OPDF (e.g., uniform is Y=1/(4*pi)
%   everywhere)
%   tstart: absolute start time of simulation
%   store: number of timesteps to store information, not currently
%       implemented to store, keep large
%   tstop: stop time of simulation
%   dt: timestep size



%Solver Parameters
Cr=1;%Courant Number, keep as 1
alpha=0.5;%Solver Tuning Parameter (0.5 for Crank Nicholson) - keep as 0.5

%matrix organization to setup 1D system of linear equations
nphi=nangle;
ntheta=nangle;
x=zeros((nphi)*(ntheta),1);
Ax=zeros((nphi)*(ntheta),1);
Ay=Ax;
Aval=Ax;
xinitial=x;
b=x;
r=zeros(3,1);
t=r;
p=r;
Fe=zeros(ntheta,nphi);
Fw=Fe;
Fn=Fe;
Fs=Fe;
De=Fe;
Dw=Fe;
Dn=Fe;
Ds=Fe;
Y=Fe;
psi=Fe;
am=Fe;
bmat=Fe;
ap=zeros(ntheta,1);
phivect=ap;
thetavect=ap;
Ytot=0;
a2=zeros(6,1);
a4=zeros(15,1);



dtheta=pi/ntheta;
dphi=pi/nphi;





%Setup Fluxes and Initial Condition
for Jp=1:nphi %create phi and theta vectors
    JJ=Jp-0.5;
    phi=(Jp-0.5)*dphi;
    phivect(Jp)=phi*180/pi;
    for Jt=1:ntheta
        theta=(Jt-0.5)*dtheta;
        thetavect(Jt)=theta*180/pi;
        II=Jt-0.5;
        
        r(1)=sin(theta)*cos(phi);
        r(2)=sin(theta)*sin(phi);
        r(3)=cos(theta);
        
        t(1)=cos(theta)*cos(phi);
        t(2)=cos(theta)*sin(phi);
        t(3)=-sin(theta);
        
        p(1)=-sin(phi);
        p(2)=cos(phi);
        p(3)=0;
        
        %Time rate of change of theta and phi
        tdot = (ARG+1)*0.5*(...
            (k(1,1)*r(1)+k(1,2)*r(2)+k(1,3)*r(3))*t(1)+...
            (k(2,1)*r(1)+k(2,2)*r(2)+k(2,3)*r(3))*t(2)+...
            (k(3,1)*r(1)+k(3,2)*r(2)+k(3,3)*r(3))*t(3))+...
            (ARG-1)*0.5*(...
            (k(1,1)*r(1)+k(2,1)*r(2)+k(3,1)*r(3))*t(1)+...
            (k(1,2)*r(1)+k(2,2)*r(2)+k(3,2)*r(3))*t(2)+...
            (k(1,3)*r(1)+k(2,3)*r(2)+k(3,3)*r(3))*t(3));
        
        pdot = (ARG+1)*0.5*(...
            (k(1,1)*r(1)+k(1,2)*r(2)+k(1,3)*r(3))*p(1)+...
            (k(2,1)*r(1)+k(2,2)*r(2)+k(2,3)*r(3))*p(2)+...
            (k(3,1)*r(1)+k(3,2)*r(2)+k(3,3)*r(3))*p(3))+...
            (ARG-1)*0.5*(...
            (k(1,1)*r(1)+k(2,1)*r(2)+k(3,1)*r(3))*p(1)+...
            (k(1,2)*r(1)+k(2,2)*r(2)+k(3,2)*r(3))*p(2)+...
            (k(1,3)*r(1)+k(2,3)*r(2)+k(3,3)*r(3))*p(3));
        
        %Calculate Timestep
        if tdot<dtheta/dt
            tdot=dtheta/dt;
        end
        if abs(pdot)<dphi/dt
            pdot=dphi/dt;
        end
        
        dt=min([abs(Cr*dtheta/tdot),abs(Cr*dphi/pdot),dt]);%Timestep
        steps=tstop/dt;%calculate number of steps required
        
        %East Flux, theta+0.5
        
        r(1)=sin((II+0.5)*dtheta)*cos(phi);
        r(2)=sin((II+0.5)*dtheta)*sin(phi);
        r(3)=cos((II+0.5)*dtheta);
        
        t(1)=cos((II+0.5)*dtheta)*cos(phi);
        t(2)=cos((II+0.5)*dtheta)*sin(phi);
        t(3)=-sin((II+0.5)*dtheta);
        
        Fe(Jt,Jp)=-sin((II+0.5)*dtheta)*dphi*(ARG+1)*0.25*(...
            (k(1,1)*r(1)+k(1,2)*r(2)+k(1,3)*r(3))*t(1)+...
            (k(2,1)*r(1)+k(2,2)*r(2)+k(2,3)*r(3))*t(2)+...
            (k(3,1)*r(1)+k(3,2)*r(2)+k(3,3)*r(3))*t(3))+...
            -sin((II+0.5)*dtheta)*dphi*(ARG-1)*0.25*(...
            (k(1,1)*r(1)+k(2,1)*r(2)+k(3,1)*r(3))*t(1)+...
            (k(1,2)*r(1)+k(2,2)*r(2)+k(3,2)*r(3))*t(2)+...
            (k(1,3)*r(1)+k(2,3)*r(2)+k(3,3)*r(3))*t(3));
        
        De(Jt,Jp)=Dr*dphi*sin((II+0.5)*dtheta)/dtheta;
        
        if Jt==ntheta
            Fe(Jt,Jp)=0;
            De(Jt,Jp)=0;
        end
        
        %West Flux, theta-0.5
        
        r(1)=sin((II-0.5)*dtheta)*cos(phi);
        r(2)=sin((II-0.5)*dtheta)*sin(phi);
        r(3)=cos((II-0.5)*dtheta);
        
        t(1)=cos((II-0.5)*dtheta)*cos(phi);
        t(2)=cos((II-0.5)*dtheta)*sin(phi);
        t(3)=-sin((II-0.5)*dtheta);
        
        Fw(Jt,Jp)=sin((II-0.5)*dtheta)*dphi*(ARG+1)*0.25*(...
            (k(1,1)*r(1)+k(1,2)*r(2)+k(1,3)*r(3))*t(1)+...
            (k(2,1)*r(1)+k(2,2)*r(2)+k(2,3)*r(3))*t(2)+...
            (k(3,1)*r(1)+k(3,2)*r(2)+k(3,3)*r(3))*t(3))+...
            sin((II-0.5)*dtheta)*dphi*(ARG-1)*0.25*(...
            (k(1,1)*r(1)+k(2,1)*r(2)+k(3,1)*r(3))*t(1)+...
            (k(1,2)*r(1)+k(2,2)*r(2)+k(3,2)*r(3))*t(2)+...
            (k(1,3)*r(1)+k(2,3)*r(2)+k(3,3)*r(3))*t(3));
        Dw(Jt,Jp)=Dr*dphi*sin((II-0.5)*dtheta)/dtheta;
        
        if Jt==1
            Fw(Jt,Jp)=0;
            Dw(Jt,Jp)=0;
        end
        
        %North Flux, phi+0.5
        
        r(1)=sin(theta)*cos((JJ+0.5)*dphi);
        r(2)=sin(theta)*sin((JJ+0.5)*dphi);
        r(3)=cos(theta);
        
        p(1)=-sin((JJ+0.5)*dphi);
        p(2)=cos((JJ+0.5)*dphi);
        p(3)=0;
        
        Fn(Jt,Jp)=-dtheta*(ARG+1)*0.25*(...
            (k(1,1)*r(1)+k(1,2)*r(2)+k(1,3)*r(3))*p(1)+...
            (k(2,1)*r(1)+k(2,2)*r(2)+k(2,3)*r(3))*p(2)+...
            (k(3,1)*r(1)+k(3,2)*r(2)+k(3,3)*r(3))*p(3))+...
            -dtheta*(ARG-1)*0.25*(...
            (k(1,1)*r(1)+k(2,1)*r(2)+k(3,1)*r(3))*p(1)+...
            (k(1,2)*r(1)+k(2,2)*r(2)+k(3,2)*r(3))*p(2)+...
            (k(1,3)*r(1)+k(2,3)*r(2)+k(3,3)*r(3))*p(3));
        Dn(Jt,Jp)=Dr*dtheta/sin(theta)/dphi;
        
        %South Flux, phi-0.5
        
        r(1)=sin(theta)*cos((JJ-0.5)*dphi);
        r(2)=sin(theta)*sin((JJ-0.5)*dphi);
        %r(3)=cos(theta);
        
        p(1)=-sin((JJ-0.5)*dphi);
        p(2)=cos((JJ-0.5)*dphi);
        %p(3)=0;
        
        Fs(Jt,Jp)=dtheta*(ARG+1)*0.25*(...
            (k(1,1)*r(1)+k(1,2)*r(2)+k(1,3)*r(3))*p(1)+...
            (k(2,1)*r(1)+k(2,2)*r(2)+k(2,3)*r(3))*p(2)+...
            (k(3,1)*r(1)+k(3,2)*r(2)+k(3,3)*r(3))*p(3))+...
            dtheta*(ARG-1)*0.25*(...
            (k(1,1)*r(1)+k(2,1)*r(2)+k(3,1)*r(3))*p(1)+...
            (k(1,2)*r(1)+k(2,2)*r(2)+k(3,2)*r(3))*p(2)+...
            (k(1,3)*r(1)+k(2,3)*r(2)+k(3,3)*r(3))*p(3));
        Ds(Jt,Jp)=Dr*dtheta/sin(theta)/dphi;
        
        %Initial Distribution Function
        
        r(1)=sin(theta)*cos(phi);
        r(2)=sin(theta)*sin(phi);
        r(3)=cos(theta);
        
        
        
        if Jp==1 %only need to evaluate once per phi
            ap(Jt)=sin(theta)*dtheta*dphi/dt;
        end
        
        am(Jt,Jp)=Fs(Jt,Jp)-Ds(Jt,Jp)+Fn(Jt,Jp)-Dn(Jt,Jp)+Fw(Jt,Jp)-Dw(Jt,Jp)+...
            Fe(Jt,Jp)-De(Jt,Jp);
        
        %A diagonal creation
        Ax(Jt+ntheta*(Jp-1))=Jt+ntheta*(Jp-1);
        Ay(Jt+ntheta*(Jp-1))=Jt+ntheta*(Jp-1);
        Aval(Jt+ntheta*(Jp-1))=ap(Jt);
        
        
    end
end

Y=Yinitcond; %Set initial condition

dt
steps

%shearstress=zeros(ceil(steps),2);

%Create A Matrix
%Center is 0x, east is 1x, west is 2x, north is 3x, south is 4x
%Case 1: non-boundaries
for j=2:nphi-1
    for i=2:ntheta-1
        bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(i+1,j)+...%East is theta+1
            (Fw(i,j)+Dw(i,j))*Y(i-1,j)+...%West is theta-1
            (Fn(i,j)+Dn(i,j))*Y(i,j+1)+...%North is phi+1
            (Fs(i,j)+Ds(i,j))*Y(i,j-1))+...%South is phi-1
            (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
        
        Ax(i+ntheta*(j-1))=i+ntheta*(j-1);%Center
        Ay(i+ntheta*(j-1))=i+ntheta*(j-1);
        Aval(i+ntheta*(j-1))=ap(i)-alpha*am(i,j);
        
        Ax(1*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%East
        Ay(1*nphi*ntheta+(i+ntheta*(j-1)))=i+1+ntheta*(j-1);
        Aval(1*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fe(i,j)+De(i,j));
        
        Ax(2*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%West
        Ay(2*nphi*ntheta+(i+ntheta*(j-1)))=i-1+ntheta*(j-1);
        Aval(2*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fw(i,j)+Dw(i,j));
        
        Ax(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%North
        Ay(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j+1-1);
        Aval(3*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fn(i,j)+Dn(i,j));
        
        Ax(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%South
        Ay(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1-1);
        Aval(4*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fs(i,j)+Ds(i,j));
    end
end
%Case 2:Upper phi boundary
j=nphi;
for i=2:ntheta-1
    bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(i+1,j)+...%East is theta+1
        (Fw(i,j)+Dw(i,j))*Y(i-1,j)+...%West is theta-1
        (Fn(i,j)+Dn(i,j))*Y(ntheta-i+1,1)+...%North is phi+1, reflective boundary conditions
        (Fs(i,j)+Ds(i,j))*Y(i,j-1))+...%South is phi-1
        (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
    
    Ax(i+ntheta*(j-1))=i+ntheta*(j-1);%Center
    Ay(i+ntheta*(j-1))=i+ntheta*(j-1);
    Aval(i+ntheta*(j-1))=ap(i)-alpha*am(i,j);
    
    Ax(1*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%East
    Ay(1*nphi*ntheta+(i+ntheta*(j-1)))=i+1+ntheta*(j-1);
    Aval(1*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fe(i,j)+De(i,j));
    
    Ax(2*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%West
    Ay(2*nphi*ntheta+(i+ntheta*(j-1)))=i-1+ntheta*(j-1);
    Aval(2*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fw(i,j)+Dw(i,j));
    
    Ax(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%North
    Ay(3*nphi*ntheta+(i+ntheta*(j-1)))=ntheta-i+1+ntheta*(1-1);
    Aval(3*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fn(i,j)+Dn(i,j));
    
    Ax(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%South
    Ay(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1-1);
    Aval(4*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fs(i,j)+Ds(i,j));
end
%Case 3:Lower phi boundary
j=1;
for i=2:ntheta-1
    bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(i+1,j)+...%East is theta+1
        (Fw(i,j)+Dw(i,j))*Y(i-1,j)+...%West is theta-1
        (Fn(i,j)+Dn(i,j))*Y(i,j+1)+...%North is phi+1
        (Fs(i,j)+Ds(i,j))*Y(ntheta-i+1,nphi))+...%South is phi-1, reflective boundary condition
        (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
    
    Ax(i+ntheta*(j-1))=i+ntheta*(j-1);%Center
    Ay(i+ntheta*(j-1))=i+ntheta*(j-1);
    Aval(i+ntheta*(j-1))=ap(i)-alpha*am(i,j);
    
    Ax(1*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%East
    Ay(1*nphi*ntheta+(i+ntheta*(j-1)))=i+1+ntheta*(j-1);
    Aval(1*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fe(i,j)+De(i,j));
    
    Ax(2*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%West
    Ay(2*nphi*ntheta+(i+ntheta*(j-1)))=i-1+ntheta*(j-1);
    Aval(2*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fw(i,j)+Dw(i,j));
    
    Ax(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%North
    Ay(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j+1-1);
    Aval(3*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fn(i,j)+Dn(i,j));
    
    Ax(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%South
    Ay(4*nphi*ntheta+(i+ntheta*(j-1)))=ntheta-i+1+ntheta*(nphi-1);
    Aval(4*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fs(i,j)+Ds(i,j));
end
%Case 4:Lower theta boundary
i=1;
for j=2:nphi-1
    bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(i+1,j)+...%East is theta+1
        (Fw(i,j)+Dw(i,j))*Y(ntheta,j)+...%West is theta-1, periodic boundary condition
        (Fn(i,j)+Dn(i,j))*Y(i,j+1)+...%North is phi+1
        (Fs(i,j)+Ds(i,j))*Y(i,j-1))+...%South is phi-1
        (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
    
    Ax(i+ntheta*(j-1))=i+ntheta*(j-1);%Center
    Ay(i+ntheta*(j-1))=i+ntheta*(j-1);
    Aval(i+ntheta*(j-1))=ap(i)-alpha*am(i,j);
    
    Ax(1*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%East
    Ay(1*nphi*ntheta+(i+ntheta*(j-1)))=i+1+ntheta*(j-1);
    Aval(1*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fe(i,j)+De(i,j));
    
    Ax(2*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%West
    Ay(2*nphi*ntheta+(i+ntheta*(j-1)))=ntheta+ntheta*(j-1);
    Aval(2*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fw(i,j)+Dw(i,j));
    
    Ax(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%North
    Ay(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j+1-1);
    Aval(3*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fn(i,j)+Dn(i,j));
    
    Ax(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%South
    Ay(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1-1);
    Aval(4*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fs(i,j)+Ds(i,j));
end
%Case 5:Upper theta boundary
i=ntheta;
for j=2:nphi-1
    bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(1,j)+...%East is theta+1, periodic boundary condition
        (Fw(i,j)+Dw(i,j))*Y(i-1,j)+...%West is theta-1
        (Fn(i,j)+Dn(i,j))*Y(i,j+1)+...%North is phi+1
        (Fs(i,j)+Ds(i,j))*Y(i,j-1))+...%South is phi-1
        (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
    
    Ax(i+ntheta*(j-1))=i+ntheta*(j-1);%Center
    Ay(i+ntheta*(j-1))=i+ntheta*(j-1);
    Aval(i+ntheta*(j-1))=ap(i)-alpha*am(i,j);
    
    Ax(1*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%East
    Ay(1*nphi*ntheta+(i+ntheta*(j-1)))=1+ntheta*(j-1);
    Aval(1*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fe(i,j)+De(i,j));
    
    Ax(2*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%West
    Ay(2*nphi*ntheta+(i+ntheta*(j-1)))=i-1+ntheta*(j-1);
    Aval(2*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fw(i,j)+Dw(i,j));
    
    Ax(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%North
    Ay(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j+1-1);
    Aval(3*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fn(i,j)+Dn(i,j));
    
    Ax(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%South
    Ay(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1-1);
    Aval(4*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fs(i,j)+Ds(i,j));
end
%Case 6:Lower Left Corner
j=1;
i=1;
bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(i+1,j)+...%East is theta+1
    (Fw(i,j)+Dw(i,j))*Y(ntheta,j)+...%West is theta-1, periodic boundary condition
    (Fn(i,j)+Dn(i,j))*Y(i,j+1)+...%North is phi+1
    (Fs(i,j)+Ds(i,j))*Y(ntheta-i+1,nphi))+...%South is phi-1, reflective boundary condition
    (ap(i)+(1-alpha)*am(i,j))*Y(i,j);

Ax(i+ntheta*(j-1))=i+ntheta*(j-1);%Center
Ay(i+ntheta*(j-1))=i+ntheta*(j-1);
Aval(i+ntheta*(j-1))=ap(i)-alpha*am(i,j);

Ax(1*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%East
Ay(1*nphi*ntheta+(i+ntheta*(j-1)))=i+1+ntheta*(j-1);
Aval(1*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fe(i,j)+De(i,j));

Ax(2*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%West
Ay(2*nphi*ntheta+(i+ntheta*(j-1)))=ntheta+ntheta*(j-1);
Aval(2*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fw(i,j)+Dw(i,j));

Ax(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%North
Ay(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j+1-1);
Aval(3*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fn(i,j)+Dn(i,j));

Ax(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%South
Ay(4*nphi*ntheta+(i+ntheta*(j-1)))=ntheta-i+1+ntheta*(nphi-1);
Aval(4*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fs(i,j)+Ds(i,j));
%Case 7:Lower Right Corner
j=1;
i=ntheta;
bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(1,j)+...%East is theta+1, periodic boundary condition
    (Fw(i,j)+Dw(i,j))*Y(i-1,j)+...%West is theta-1
    (Fn(i,j)+Dn(i,j))*Y(i,j+1)+...%North is phi+1
    (Fs(i,j)+Ds(i,j))*Y(ntheta-i+1,nphi))+...%South is phi-1, reflective boundary condition
    (ap(i)+(1-alpha)*am(i,j))*Y(i,j);

Ax(i+ntheta*(j-1))=i+ntheta*(j-1);%Center
Ay(i+ntheta*(j-1))=i+ntheta*(j-1);
Aval(i+ntheta*(j-1))=ap(i)-alpha*am(i,j);

Ax(1*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%East
Ay(1*nphi*ntheta+(i+ntheta*(j-1)))=1+ntheta*(j-1);
Aval(1*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fe(i,j)+De(i,j));

Ax(2*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%West
Ay(2*nphi*ntheta+(i+ntheta*(j-1)))=i-1+ntheta*(j-1);
Aval(2*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fw(i,j)+Dw(i,j));

Ax(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%North
Ay(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j+1-1);
Aval(3*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fn(i,j)+Dn(i,j));

Ax(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%South
Ay(4*nphi*ntheta+(i+ntheta*(j-1)))=ntheta-i+1+ntheta*(nphi-1);
Aval(4*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fs(i,j)+Ds(i,j));
%Case 8:Upper Left Corner
j=nphi;
i=1;
bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(i+1,j)+...%East is theta+1
    (Fw(i,j)+Dw(i,j))*Y(ntheta,j)+...%West is theta-1, periodic boundary condition
    (Fn(i,j)+Dn(i,j))*Y(ntheta-i+1,1)+...%North is phi+1, reflective boundary condition
    (Fs(i,j)+Ds(i,j))*Y(i,j-1))+...%South is phi-1
    (ap(i)+(1-alpha)*am(i,j))*Y(i,j);

Ax(i+ntheta*(j-1))=i+ntheta*(j-1);%Center
Ay(i+ntheta*(j-1))=i+ntheta*(j-1);
Aval(i+ntheta*(j-1))=ap(i)-alpha*am(i,j);

Ax(1*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%East
Ay(1*nphi*ntheta+(i+ntheta*(j-1)))=i+1+ntheta*(j-1);
Aval(1*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fe(i,j)+De(i,j));

Ax(2*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%West
Ay(2*nphi*ntheta+(i+ntheta*(j-1)))=ntheta+ntheta*(j-1);
Aval(2*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fw(i,j)+Dw(i,j));

Ax(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%North
Ay(3*nphi*ntheta+(i+ntheta*(j-1)))=ntheta-i+1+ntheta*(1-1);
Aval(3*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fn(i,j)+Dn(i,j));

Ax(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%South
Ay(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1-1);
Aval(4*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fs(i,j)+Ds(i,j));
%Case 9:Upper Right Corner
j=nphi;
i=ntheta;
bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(1,j)+...%East is theta+1, periodic boundary condition
    (Fw(i,j)+Dw(i,j))*Y(i-1,j)+...%West is theta-1
    (Fn(i,j)+Dn(i,j))*Y(ntheta-i+1,1)+...%North is phi+1, reflective boundary condition
    (Fs(i,j)+Ds(i,j))*Y(i,j-1))+...%South is phi-1
    (ap(i)+(1-alpha)*am(i,j))*Y(i,j);

Ax(i+ntheta*(j-1))=i+ntheta*(j-1);%Center
Ay(i+ntheta*(j-1))=i+ntheta*(j-1);
Aval(i+ntheta*(j-1))=ap(i)-alpha*am(i,j);

Ax(1*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%East
Ay(1*nphi*ntheta+(i+ntheta*(j-1)))=1+ntheta*(j-1);
Aval(1*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fe(i,j)+De(i,j));

Ax(2*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%West
Ay(2*nphi*ntheta+(i+ntheta*(j-1)))=i-1+ntheta*(j-1);
Aval(2*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fw(i,j)+Dw(i,j));

Ax(3*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%North
Ay(3*nphi*ntheta+(i+ntheta*(j-1)))=ntheta-i+1+ntheta*(1-1);
Aval(3*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fn(i,j)+Dn(i,j));

Ax(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1);%South
Ay(4*nphi*ntheta+(i+ntheta*(j-1)))=i+ntheta*(j-1-1);
Aval(4*nphi*ntheta+(i+ntheta*(j-1)))=-alpha*(Fs(i,j)+Ds(i,j));

A=sparse(Ax,Ay,Aval);%generate sparse matrix A
clear Aval Ax Ay


%Calculate Initial Moments of Distribution
for Jp=1:nphi
    JJ=Jp-0.5;
    phi=(Jp-0.5)*dphi;
    for Jt=1:ntheta
        theta=(Jt-0.5)*dtheta;
        II=Jt-0.5;
        r(1)=sin(theta)*cos(phi);
        r(2)=sin(theta)*sin(phi);
        r(3)=cos(theta);
        
        Ytot=Ytot+2*Y(Jt,Jp)*dtheta*dphi*sin(theta);
        
        
        a2(1)=a2(1)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1);
        a2(2)=a2(2)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(2);
        a2(3)=a2(3)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(3);
        a2(4)=a2(4)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(2);
        a2(5)=a2(5)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(3);
        a2(6)=a2(6)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(3)*r(3);
        
        a4(1)=a4(1)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(1)*r(1);
        a4(2)=a4(2)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(2)*r(2)*r(2);
        a4(3)=a4(3)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(3)*r(3)*r(3)*r(3);
        a4(4)=a4(4)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(2)*r(2);
        a4(5)=a4(5)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(3)*r(3);
        a4(6)=a4(6)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(2)*r(3)*r(3);
        a4(7)=a4(7)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(2)*r(2)*r(2);
        a4(8)=a4(8)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(3)*r(3)*r(3);
        a4(9)=a4(9)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(3)*r(3)*r(3);
        a4(10)=a4(10)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(1)*r(2);
        a4(11)=a4(11)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(1)*r(3);
        a4(12)=a4(12)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(2)*r(2)*r(3);
        a4(13)=a4(13)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(2)*r(3)*r(3);
        a4(14)=a4(14)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(2)*r(2)*r(3);
        a4(15)=a4(15)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(2)*r(3);
    end
end

Ytotinit=Ytot;
time=tstart;
count=0;
stresscount=0;

while time<tstop
    count=count+1;
    time=time+dt;
    %Build b vector
    %Case 1: non-boundaries
    for j=2:nphi-1
        for i=2:ntheta-1
            bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(i+1,j)+...%East is theta+1
                (Fw(i,j)+Dw(i,j))*Y(i-1,j)+...%West is theta-1
                (Fn(i,j)+Dn(i,j))*Y(i,j+1)+...%North is phi+1
                (Fs(i,j)+Ds(i,j))*Y(i,j-1))+...%South is phi-1
                (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
        end
    end
    %Case 2:Upper phi boundary
    j=nphi;
    for i=2:ntheta-1
        bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(i+1,j)+...%East is theta+1
            (Fw(i,j)+Dw(i,j))*Y(i-1,j)+...%West is theta-1
            (Fn(i,j)+Dn(i,j))*Y(ntheta-i+1,1)+...%North is phi+1, reflective boundary conditions
            (Fs(i,j)+Ds(i,j))*Y(i,j-1))+...%South is phi-1
            (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
    end
    %Case 3:Lower phi boundary
    j=1;
    for i=2:ntheta-1
        bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(i+1,j)+...%East is theta+1
            (Fw(i,j)+Dw(i,j))*Y(i-1,j)+...%West is theta-1
            (Fn(i,j)+Dn(i,j))*Y(i,j+1)+...%North is phi+1
            (Fs(i,j)+Ds(i,j))*Y(ntheta-i+1,nphi))+...%South is phi-1, reflective boundary condition
            (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
    end
    %Case 4:Lower theta boundary
    i=1;
    for j=2:nphi-1
        bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(i+1,j)+...%East is theta+1
            (Fw(i,j)+Dw(i,j))*Y(ntheta,j)+...%West is theta-1, periodic boundary condition
            (Fn(i,j)+Dn(i,j))*Y(i,j+1)+...%North is phi+1
            (Fs(i,j)+Ds(i,j))*Y(i,j-1))+...%South is phi-1
            (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
    end
    %Case 5:Upper theta boundary
    i=ntheta;
    for j=2:nphi-1
        bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(1,j)+...%East is theta+1, periodic boundary condition
            (Fw(i,j)+Dw(i,j))*Y(i-1,j)+...%West is theta-1
            (Fn(i,j)+Dn(i,j))*Y(i,j+1)+...%North is phi+1
            (Fs(i,j)+Ds(i,j))*Y(i,j-1))+...%South is phi-1
            (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
    end
    %Case 6:Lower Left Corner
    j=1;
    i=1;
    bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(i+1,j)+...%East is theta+1
        (Fw(i,j)+Dw(i,j))*Y(ntheta,j)+...%West is theta-1, periodic boundary condition
        (Fn(i,j)+Dn(i,j))*Y(i,j+1)+...%North is phi+1
        (Fs(i,j)+Ds(i,j))*Y(ntheta-i+1,nphi))+...%South is phi-1, reflective boundary condition
        (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
    %Case 7:Lower Right Corner
    j=1;
    i=ntheta;
    bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(1,j)+...%East is theta+1, periodic boundary condition
        (Fw(i,j)+Dw(i,j))*Y(i-1,j)+...%West is theta-1
        (Fn(i,j)+Dn(i,j))*Y(i,j+1)+...%North is phi+1
        (Fs(i,j)+Ds(i,j))*Y(ntheta-i+1,nphi))+...%South is phi-1, reflective boundary condition
        (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
    %Case 8:Upper Left Corner
    j=nphi;
    i=1;
    bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(i+1,j)+...%East is theta+1
        (Fw(i,j)+Dw(i,j))*Y(ntheta,j)+...%West is theta-1, periodic boundary condition
        (Fn(i,j)+Dn(i,j))*Y(ntheta-i+1,1)+...%North is phi+1, reflective boundary condition
        (Fs(i,j)+Ds(i,j))*Y(i,j-1))+...%South is phi-1
        (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
    %Case 9:Upper Right Corner
    j=nphi;
    i=ntheta;
    bmat(i,j)=(1-alpha)*((Fe(i,j)+De(i,j))*Y(1,j)+...%East is theta+1, periodic boundary condition
        (Fw(i,j)+Dw(i,j))*Y(i-1,j)+...%West is theta-1
        (Fn(i,j)+Dn(i,j))*Y(ntheta-i+1,1)+...%North is phi+1, reflective boundary condition
        (Fs(i,j)+Ds(i,j))*Y(i,j-1))+...%South is phi-1
        (ap(i)+(1-alpha)*am(i,j))*Y(i,j);
    
    b=reshape(bmat,nphi*ntheta,1);
    x=A\b;%Solve linear system of equations
    
    Ynew=reshape(x,nphi,ntheta);
    Y=Ynew;
    if any(Y<0)==1
        disp('Error: distribution has negative values, try increasing nangle')
        time
        count
        break
    end
    
    if mod(count,store)==0
        stresscount=stresscount+1;
        Ytot=0;
        a2=zeros(6,1);
        a4=zeros(15,1);
        for Jp=1:nphi
            JJ=Jp-0.5;
            phi=(Jp-0.5)*dphi;
            for Jt=1:ntheta
                theta=(Jt-0.5)*dtheta;
                II=Jt-0.5;
                r(1)=sin(theta)*cos(phi);
                r(2)=sin(theta)*sin(phi);
                r(3)=cos(theta);
                
                Ytot=Ytot+2*Y(Jt,Jp)*dtheta*dphi*sin(theta);
                
                a2(1)=a2(1)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1);
                a2(2)=a2(2)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(2);
                a2(3)=a2(3)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(3);
                a2(4)=a2(4)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(2);
                a2(5)=a2(5)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(3);
                a2(6)=a2(6)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(3)*r(3);
                
                a4(1)=a4(1)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(1)*r(1);
                a4(2)=a4(2)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(2)*r(2)*r(2);
                a4(3)=a4(3)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(3)*r(3)*r(3)*r(3);
                a4(4)=a4(4)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(2)*r(2);
                a4(5)=a4(5)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(3)*r(3);
                a4(6)=a4(6)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(2)*r(3)*r(3);
                a4(7)=a4(7)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(2)*r(2)*r(2);
                a4(8)=a4(8)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(3)*r(3)*r(3);
                a4(9)=a4(9)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(3)*r(3)*r(3);
                a4(10)=a4(10)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(1)*r(2);
                a4(11)=a4(11)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(1)*r(3);
                a4(12)=a4(12)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(2)*r(2)*r(3);
                a4(13)=a4(13)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(2)*r(3)*r(3);
                a4(14)=a4(14)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(2)*r(2)*r(3);
                a4(15)=a4(15)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(2)*r(3);
                
                psi(Jt,Jp)=sin(theta)*Y(Jt,Jp);
            end
        end
        shearstress(stresscount,1)=time;
        
        Aconst=AR^2/(2*(log(2*AR)-3/2));
        Bconst=(6*log(2*AR)-11)/AR^2;
        Cconst=2;
        Fconst=3*AR^2/(log(2*AR)-1/2);
        
        shearstress(stresscount,2)=4*Aconst*G*a4(4)+2*Bconst*G*(a2(1)+a2(4))+Cconst*G+2*Fconst*Dr*a2(2);
        %time
    end
end

Ytot=0;
a2=zeros(6,1);
a4=zeros(15,1);
for Jp=1:nphi
    JJ=Jp-0.5;
    phi=(Jp-0.5)*dphi;
    for Jt=1:ntheta
        theta=(Jt-0.5)*dtheta;
        II=Jt-0.5;
        r(1)=sin(theta)*cos(phi);
        r(2)=sin(theta)*sin(phi);
        r(3)=cos(theta);
        
        Ytot=Ytot+2*Y(Jt,Jp)*dtheta*dphi*sin(theta);
        
        a2(1)=a2(1)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1);
        a2(2)=a2(2)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(2);
        a2(3)=a2(3)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(3);
        a2(4)=a2(4)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(2);
        a2(5)=a2(5)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(3);
        a2(6)=a2(6)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(3)*r(3);
        
        a4(1)=a4(1)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(1)*r(1);
        a4(2)=a4(2)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(2)*r(2)*r(2);
        a4(3)=a4(3)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(3)*r(3)*r(3)*r(3);
        a4(4)=a4(4)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(2)*r(2);
        a4(5)=a4(5)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(3)*r(3);
        a4(6)=a4(6)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(2)*r(3)*r(3);
        a4(7)=a4(7)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(2)*r(2)*r(2);
        a4(8)=a4(8)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(3)*r(3)*r(3);
        a4(9)=a4(9)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(3)*r(3)*r(3);
        a4(10)=a4(10)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(1)*r(2);
        a4(11)=a4(11)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(1)*r(3);
        a4(12)=a4(12)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(2)*r(2)*r(2)*r(3);
        a4(13)=a4(13)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(2)*r(3)*r(3);
        a4(14)=a4(14)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(2)*r(2)*r(3);
        a4(15)=a4(15)+2*Y(Jt,Jp)*dtheta*dphi*sin(theta)*r(1)*r(1)*r(2)*r(3);
        
        
        psi(Jt,Jp)=2*dtheta*dphi*sin(theta)*Y(Jt,Jp);
    end
end
end