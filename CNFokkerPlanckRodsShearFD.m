%PTC last updated 4/22/2020

%Example use of function FokkerPlanckDiluteRodSolver
%Simple shear flow, Per=100/17, AR=78


clear all

%flow and discretization parameters
nangle=20;%angle discretization
Dr=17; %same units as deformation rate
G=100;%Deformation Rate
Per=G/Dr;
AR=78;%Particle Effective Aspect Ratio
ARG=(AR.^2-1)/(AR.^2+1);%Shape Factor
tstart=0;
store=10000;%store results per this many timesteps
strain=max(100,0.1*pi*AR);%strain units
tstop=strain./G; %simulation time
%tstop=max(1/(Dr*AR^2),pi.*AR./G); %for large Per
dt=1000;%Only used for no flow calculations, otherwise keep large
%dt=0.01;%uncomment when in no flow

%Set initial condition (can be from another simulation ;])
Yinitcond=ones(nangle)./(4*pi); %uniform (equilibrium) initial condition

%Velocity Gradient Tensor

% k = [k11 k12 k13]
%     [k21 k22 k23]
%     [k31 k32 k33]

k11=0;
k12=G;
k13=0;
k21=0;
k22=0;
k23=0;
k31=0;
k32=0;
k33=0;

k=[k11 k12 k13; k21 k22 k23; k31 k32 k33];

[Y,a2,a4,Ytot,phivect,thetavect] = FokkerPlanckDiluteRodSolver(nangle,k,Dr,ARG,Yinitcond,tstart,store,tstop,dt);

%if you want to calculate the shear stress from slender body theory

% Aconst=AR^2/(2*(log(2*AR)-3/2));
% Bconst=(6*log(2*AR)-11)/AR^2;
% Cconst=2;
% Fconst=3*AR^2/(log(2*AR)-1/2);

%shearstressfinal=4*Aconst*G*a4(4)+2*Bconst*G*(a2(1)+a2(4))+Cconst*G+2*Fconst*Dr*a2(2);


%plot OPDF
dangle=pi/nangle;
anglecount=1;
for Jp=1:nangle
    JJ=Jp-0.5;
    phi=(Jp-0.5)*dangle;%relationship between Jp and phi
    for Jt=1:nangle
        theta=(Jt-0.5)*dangle;%relationship between Jt and theta
        II=Jt-0.5;
        
        r1(anglecount,1)=sin(theta)*cos(phi+pi);
        r2(anglecount,1)=sin(theta)*sin(phi+pi);
        r3(anglecount,1)=cos(theta);
        Yvector(anglecount,1)=Y(Jt,Jp);
        
        anglecount=anglecount+1;
        r1(anglecount,1)=sin(theta)*cos(phi);
        r2(anglecount,1)=sin(theta)*sin(phi);
        r3(anglecount,1)=cos(theta);
        Yvector(anglecount,1)=Y(Jt,Jp);
        
        anglecount=anglecount+1;
    end
end

Y=Y./Ytot;
Yvector=Yvector./Ytot;

figure

simp = convhulln( [r1,r2,r3] );
trisurf(simp,r1,r2,r3,Yvector,'EdgeColor','none','FaceColor','interp');
colormap jet
xlim([-1,1]);
xlabel('u')
ylim([-1,1]);
ylabel('\nabla u')
zlim([-1,1]);
zlabel('\omega')
axis square
set(gca,'FontSize',16)
colorbar
