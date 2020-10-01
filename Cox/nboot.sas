%macro nboot(p_target=0.04, p_limit=, conlevp=, maxiter=50, nboot=, side=lower);


* Makro zur Planung der Anzahl der Bootstrap-Resamples;

* besonders bei RELIMPLR und RELIMPCR;



* side can be:  lower or upper;

%let limit=&side;

%if &nboot= %then %do;
data t_up;
z_clp=-probit((1-&conlevp)/2);
t_p_target=probit(1-&p_target/2);
t_p_limit=probit(1-&p_limit/2);
if &p_limit>&p_target then z_clp=-z_clp;
ratio=t_p_limit**2/t_p_target**2;

if &p_target<&p_limit then p=(1-&conlevp)/2;
else p=1-(1-&conlevp)/2;

n=1500;

optlag=cinv(p,n-1)/(n-1)-ratio;
*optlag=0.5*(z_clp+sqrt(2*n-3))**2/(n-1)-ratio;
nlag=n;
n=50;
*opt=0.5*(z_clp+sqrt(2*n-3))**2/(n-1)-ratio;
opt=cinv(p,n-1)/(n-1)-ratio;

iter=1;
put iter= opt= n=;
do while(iter<=&maxiter and (abs(floor(nlag)-floor(n))>0));
 iter=iter+1;
 nsave=n;
 n=nlag-(n-nlag)*(optlag)/(opt-optlag);
 if n<3 then n=3;
 nlag=nsave;
 optlag=opt;
* opt=0.5*(z_clp+sqrt(2*n-3))**2/(n-1)-ratio;
 opt=cinv(p,n-1)/(n-1)-ratio;
 put iter= opt= n=;
end;
eff=(&p_target*(1-&p_target)*z_clp**2/(&p_target-&p_limit)**2+1)/(ceil(n));
p_target=&p_target;
p_limit=&p_limit;
conlevp=&conlevp;
keep iter n eff opt p_target p_limit conlevp;
label iter="# Iterations" n="Bootstrap resamples needed" 
 eff="Efficiency of normal/nonparametric"
 p_target="Target p-value" p_limit="Limit of C.I." conlevp="Confidence level" opt="Target function";
output t_up;
run;

%end;
%else %if &conlevp= %then %do;

data t_up;
clp=0.99;
z_clp=-probit((1-clp)/2);
t_p_target=probit(1-&p_target/2);
t_p_limit=probit(1-&p_limit/2);
if &p_limit>&p_target then z_clp=-z_clp;
ratio=t_p_limit**2/t_p_target**2;

if &p_target<&p_limit then p=(1-clp)/2;
else p=1-(1-clp)/2;

n=&nboot;

optlag=cinv(p,n-1)/(n-1)-ratio;
*optlag=0.5*(z_clp+sqrt(2*n-3))**2/(n-1)-ratio;
clplag=clp;
clp=0.5;
*opt=0.5*(z_clp+sqrt(2*n-3))**2/(n-1)-ratio;
z_clp=-probit((1-clp)/2);
if p_limit>&p_target then z_clp=-z_clp;
if &p_target<&p_limit then p=(1-clp)/2;
else p=1-(1-clp)/2;

opt=cinv(p,n-1)/(n-1)-ratio;

iter=1;
put iter= opt= clp=;
do while(iter<=&maxiter and (abs(clplag-clp)>0.005));
 iter=iter+1;
 clpsave=clp;
 clp=clplag-(clp-clplag)*(optlag)/(opt-optlag);
 if clp<0.001 then clp=0.001;
 if clp>0.9999 then clp=0.9999;
 clplag=clpsave;
 optlag=opt;
* opt=0.5*(z_clp+sqrt(2*n-3))**2/(n-1)-ratio;
 if &p_target<&p_limit then p=(1-clp)/2;
 else p=1-(1-clp)/2;
 opt=cinv(p,n-1)/(n-1)-ratio;
 put iter= opt= clp=;
end;
z_clp=-probit((1-clp)/2);
if &p_limit>&p_target then z_clp=-z_clp;
eff=(&p_target*(1-&p_target)*z_clp**2/(&p_target-&p_limit)**2+1)/(ceil(n));
p_target=&p_target;
p_limit=&p_limit;
keep iter n eff opt p_target p_limit conlevp;
label iter="# Iterations" n="Bootstrap resamples needed" 
 eff="Efficiency of normal/nonparametric"
 p_target="Target p-value" p_limit="Limit of C.I." conlevp="Confidence level" opt="Target function";
conlevp=clp;
output t_up;

run;



%end;
%else %if &p_limit= %then %do;

data t_up;
%if %upcase(&limit)=UPPER %then %do;
 p_limit=min(&p_target*2, 0.99);
%end;
%else %do;
 p_limit=max(&p_target/2, 0.001);
%end;
clp=&conlevp;
z_clp=-probit((1-clp)/2);
t_p_target=probit(1-&p_target/2);
t_p_limit=probit(1-p_limit/2);
if p_limit>&p_target then z_clp=-z_clp;
ratio=t_p_limit**2/t_p_target**2;

if &p_target<p_limit then p=(1-clp)/2;
else p=1-(1-clp)/2;

n=&nboot;

optlag=cinv(p,n-1)/(n-1)-ratio;
*optlag=0.5*(z_clp+sqrt(2*n-3))**2/(n-1)-ratio;
plimlag=p_limit;
%if %upcase(&limit)=UPPER %then %do;
 p_limit=min(&p_target*1.5, 0.99);
%end;
%else %do;
 p_limit=max(&p_target/1.5, 0.001);
%end;

*opt=0.5*(z_clp+sqrt(2*n-3))**2/(n-1)-ratio;
z_clp=-probit((1-clp)/2);
if p_limit>&p_target then z_clp=-z_clp;
if &p_target<p_limit then p=(1-clp)/2;
else p=1-(1-clp)/2;
t_p_limit=probit(1-p_limit/2);
ratio=t_p_limit**2/t_p_target**2;
opt=cinv(p,n-1)/(n-1)-ratio;

iter=1;
put iter= opt= p_limit=;
do while(iter<=&maxiter and (abs(plimlag-p_limit)>0.0001));
 iter=iter+1;
 plimsave=p_limit;
 p_limit=plimlag-(p_limit-plimlag)*(optlag)/(opt-optlag);
 %if %upcase(&limit)=UPPER %then %do;
  if p_limit<&p_target then p_limit=&p_target+(&p_target-p_limit);
 %end;
 if p_limit<0.0001 then clp=0.0001;
 if p_limit>0.9999 then p_limit=0.9999;
 plimlag=plimsave;
 optlag=opt;
* opt=0.5*(z_clp+sqrt(2*n-3))**2/(n-1)-ratio;
 if &p_target<p_limit then p=(1-clp)/2;
 else p=1-(1-clp)/2;
 t_p_limit=probit(1-p_limit/2);
 ratio=t_p_limit**2/t_p_target**2;
 opt=cinv(p,n-1)/(n-1)-ratio;
 put iter= opt= p_limit=;
end;
z_clp=-probit((1-clp)/2);
if p_limit>&p_target then z_clp=-z_clp;
eff=(&p_target*(1-&p_target)*z_clp**2/(&p_target-p_limit)**2+1)/(ceil(n));
p_target=&p_target;
conlevp=clp;
output t_up;

keep iter n eff opt p_target p_limit conlevp;

label iter="# Iterations" n="Bootstrap resamples needed" 
 eff="Efficiency of normal/nonparametric"
 p_target="Target p-value" p_limit="Limit of C.I." conlevp="Confidence level" opt="Target function";



run;



%end;






proc print noobs label data=t_up;
format conlevp percent8.2;
var p_target p_limit conlevp n eff iter opt;
run;




%mend;

*example: %nboot(p_target=0.05, conlevp=0.99, p_limit=0.04);
