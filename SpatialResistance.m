%This function reproduces most results from Krieger et al. 2020, "Population structure across scales facilitates
%coexistence and spatial heterogeneity of antibiotic-resistant infections"

%Input: adjmat - an NxN unweighted undirected adjacency matrix on N demes
%       IC - the initial condition, an Nx2 matrix giving the initial condition (first column - sensitive 
%            infections, second column - resistant infections), which in the
%            paper was [1/R0, 0.001] in every deme
%       paramvec - a 5-tuple [kappa,beta,g,epsilon,c]
%       treatvec - a Nx1 vector with a 0 for untreated deme, 1 for treated
%       threshold - the number below which, when demes on average have
%                   derivatives less than this value, we will consider equilibrium 
%                   to have been reached
%
%Output: popvec - an Nx2 matrix giving the susceptible-resistant strain
%                 frequencies in each of the N demes
%        tt - the number of days (from 0) for equilibrium to have been
%             reached


function [popvec,tt]=SpatialResistance(adjmat,IC,paramvec,treatvec,threshold)
%interpret parameters
kappa=paramvec(1);
beta=paramvec(2);
g=paramvec(3);
effic=paramvec(4);
cost=paramvec(5);
tstep = 40; %The number of days to simulate before checking if equilibrium has been reached

%initialize
N=sum(size(adjmat))/2; 
popvec=IC; %should be 2 columns (susceptible and resistant) and N rows

tt=0; %Current "day" in simulation
hitswitch=0; %Will become 1 when equilibrium is reached
yold=0; %Population structure in last step

while hitswitch==0
    y0(1:N)=popvec(:,1);   
    y0(N+1:2*N)=popvec(:,2);
    opts=odeset('Events',@(t,y) EVENTFUN(t,y));
    [t,y,te,ye,ie] = ode45(@(t,y) ODEFUN(t,y),[0,tstep],y0,opts); %Solve ODEs from paper
    yend=y(end,:);
    popvec(:,1)=yend(1:N)';
    popvec(:,2)=yend(N+1:2*N)';
    tt=tt+tstep;
    
    if abs(sum(yend-yold))/sum(yend)<threshold/length(adjmat) %original value used: 0.0001/length(adjmat)
        hitswitch=1;
    else
        yold=yend;
    end


end

function dydt=ODEFUN(t,y) %The ODEs from the paper
      dydt=zeros(2*N,1); %a column vector
      for ii=1:N %susceptible strain
          dydt(ii)=-g*y(ii)+kappa*y(ii)*(1-y(ii)-y(ii+N))*(1-effic*treatvec(ii));
          for deme=1:N
              if deme~=ii
                  newterm=beta*adjmat(deme,ii)*(1-effic*treatvec(deme))*y(deme)*(1-y(ii)-y(ii+N));
                  dydt(ii)=dydt(ii)+newterm;
              end
          end
      end
      
      for jj=N+1:2*N %resistant strain
          dydt(jj)=-g*y(jj)+kappa*y(jj)*(1-y(jj)-y(jj-N))*(1-cost);
          for deme=1:N
              if deme~=jj
                  newterm=beta*adjmat(deme,jj-N)*(1-cost)*y(deme+N)*(1-y(jj)-y(jj-N));
                  dydt(jj)=dydt(jj)+newterm;
              end
          end
      end
    end

    function [position,isterminal,direction] = EVENTFUN(t,y) % to break the code if it runs away
position=t-1e6;
isterminal = 1;
direction=-1;
    end

end %of all things