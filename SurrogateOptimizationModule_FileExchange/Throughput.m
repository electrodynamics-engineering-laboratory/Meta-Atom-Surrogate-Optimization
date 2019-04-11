function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov] = Throughput(x, runlength, seed, other)

% x is a vector containing the service rates (1-by-n)
% runlength is the number of hours of simulated time to simulate
% seed is the index of the substreams to use (integer >= 1)
% other is a vector of buffer spaces before stage i (1-by-(n-1))
% Returns Mean response time, no var or gradient estimates.

%   ***************************************
%   *** Code written by German Gutierrez***
%   ***         gg92@cornell.edu        ***
%   ***************************************
% 
% other=randperm(11);
% runlength=10;
% seed=2;
% x=randperm(12);

%Returns throughput as defined in problem statement, Mean and C.I. Half
%Width of the time between the last 50 job releases.
FnGrad = NaN;
FnGradCov = NaN;
constraint = NaN;
ConstraintCov = NaN;
ConstraintGrad = NaN;
ConstraintGradCov = NaN;

if (sum(x < 0) > 0) || (runlength <= 0) || (seed <= 0) || (round(seed) ~= seed),
    fprintf('x (row vector with %d components)\nx components should be between 0 and 1\nrunlength should be positive and real,\nseed should be a positive integer\n', nAmbulances*2);
    fn = NaN;
    FnVar = NaN;
else % main simulation
  
r=x';             %rates of service times
b=other;               %buffer spaces before stage i
njobs=2050;             %Number of jobs to be released in simulation
nstages=length(x);             %Number of stages
ExitTimes=zeros(nstages,njobs+1);     %matrix with the time at which each job 
                                      %leaves stage i (row i).
%sTime=zeros(nstages,njobs);
                                      
ServiceT=zeros(nstages,njobs+1);

%NOTE: I will consider buffer space as the number of pieces waiting to be
%worked on. i.e. the current job does not count as part of the buffer space

% Generate new streams for call arrivals, call 
%ServiceStream = RandStream.create('mrg32k3a');

% Set the substream to the "seed"
%ServiceStream.Substream = seed;

% Generate random data
%OldStream = RandStream.setDefaultStream(ServiceStream); % Temporarily store old stream
rand('seed',seed);

%Generate matrix of service times. 
sTimeOne=exprnd(1/r(1),1, njobs);

% for k=2:nstages
%     sTimeTwo=exprnd(1/r(k),njobs, nstages);
% end

% Restore old random number stream
%RandStream.setDefaultStream(OldStream);


%Cycle through each job to determine the time at which it leaves every stage
for i= 1:njobs
    t=sTimeOne(1,i);
    ServiceT(1,i)=t;
    %First stage: if time at which the i-b(1)th job left stage 2 is less
    %than the time at which this job is ready to leave stage one (there is
    %buffer space available) job leaves stage one as soon as it is ready.
    if(ExitTimes(2,max(1,i-b(1))) <= ExitTimes(1,i)+t)
        ExitTimes(1,i+1)=ExitTimes(1,i)+t;
    else
        %No buffer space available, job leaves as soon as buffer space is
        %available.
        ExitTimes(1,i+1)=ExitTimes(2,max(1,i-b(1)));
    end
    
    for j=2:nstages-1
        t=exprnd(1/r(j));
        ServiceT(j,i)=t;
        %if (time at which previous order left > time this order left
        %previous stage)Then there must be a queue
        if(ExitTimes(j,i)>ExitTimes(j-1,i+1))
            %if there is buffer space available
            if(ExitTimes(j+1,max(1,i-b(j))) <= ExitTimes(j,i)+t)
                %leaves as soon as ready to leave
                ExitTimes(j,i+1)= ExitTimes(j,i)+t;
            else
                %No buffer space available, job leaves as soon as buffer space is
                %available.
                ExitTimes(j,i+1)=ExitTimes(j+1,i-b(j));
            end
        else
            %There is no queue, job starts to be worked on as soon as it
            %leaves previous stage.
            if(ExitTimes(j+1,max(1,i-b(j)))<= ExitTimes(j-1,i+1)+ t)
                %If there is buffer space available, job leaves as soon as
                %ready
                ExitTimes(j,i+1) = ExitTimes(j-1,i+1)+ t;
            else
                %No buffer space available, job leaves as soon as buffer space is
                %available.
                ExitTimes(j,i+1) = ExitTimes(j+1,i-b(j));
            end
        end
    end
        
        %last stage: All jobs leave as soon as ready
        
        %if there is not a queue (i.e. previous job finished before current job
        %leaves previous stage)
        if(ExitTimes(nstages,i)<=ExitTimes(nstages-1,i+1))
            t=exprnd(1/r(nstages));
            ServiceT(nstages,i)=t;
            %leaves as soon as ready, entered to be worked on as soon as
            %finished in previous stage
            ExitTimes(nstages,i+1)=ExitTimes(nstages-1,i+1)+t;
        else
            %Entered to be worked on after previous job was finished,
            %leaves as soon as ready.
            ExitTimes(nstages,i+1)=ExitTimes(nstages,i)+t;
        end   
end
end

%Throughput= 50/(time last job released - time at which the 50th to last job was released)
fn= 50/(ExitTimes(1,njobs+1)-ExitTimes(1,njobs+1-50));
end
