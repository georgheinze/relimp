%macro relimplr(data=_last_, y=, varlist=, cand=, class=, nboot=350, shrink=0, firth=0, seed=395748, 
          correctp=0, alpha=0.05, clp=0.99, groups=, gnames=, printclp=1, histogram=0, noint=0,
           maxit=50,  epsilon=0.0001, maxstep=5, maxhs=5);

* Version 2001-11;
* Build   200111081700;
* balanced bootstrap;
* angle transformation;
* C.L. for p-values;
* groups of prognostic factors;

%macro srcore(srcdata=_srcwork, srcy=&y, srck=1, rby=_rby, srcout=_srcout, pevname=pev);

* core program to compute PEV for a data set in by-groups       ;
* covariate names must be _src1 to _src&srck                    ;
* by variable must be numeric integer, starting with 1          ;
* and no missing values, please                                 ;

%if &srck=0 %then %do;
 data &srcout;
 _sample_=0;
 &pevname=0;
 output;
 run;
%end;
%else %do;

options nonotes nosource;
*options nomprint nomacrogen;
*filename SASCBTBL 's:\sasdll\surev.def';

data _ttt_;
x=time();
output;
run;

%if &firth=0 %then %do;
 ods listing close;
 ods output GlobalTests=_glob;

 proc logistic descending data=&srcdata outest=_srmy;
 model &srcy=%do srj=1 %to &srck; _src&srj %end; ;
 output out=_srcs pred=_pred_ xbeta=_xbeta_;
 by &rby;
 run;
 ods listing;

 %if &shrink=1 %then %do;
  data _glob;
  set _glob;
  if Test="Likelihood Ratio";
  _shrink_=(ChiSq-DF)/ChiSq;
  run;


  data _srmy;
  set _srmy;
  keep &rby %if &noint=0 %then %do; Intercept %end;;
  run;

  data _srcs;
  merge _srcs _glob _srmy;
  by &rby;
  _xbeta_=(_xbeta_-Intercept)*_shrink_+Intercept;
  _pred_=1/(1+exp(-_xbeta_));
  run;
 %end;


 proc corr data=_srcs outp=&srcout noprint;
 var &srcy _pred_;
 by &rby;
 run;


 data &srcout;
 set &srcout;
 &pevname=&y * &y;
 if _name_="_pred_";
 keep &pevname &rby;
 run;
%end;
%else %do;

%if &noint=0 %then %do;
 data &srcdata;
 set &srcdata;
 intercep=1;
 run;
%end;


%let srcnvar=%eval(&srck+1-&noint);

proc iml;

*******************************************************************;
start pl_comp;
 do i=1 to n;
  pi[i]=1/(1+exp(-x[i,]*beta));
*  wvec[i]=pi[i]*(1-pi[i])*(1+hvec[i]);
  wvec[i]=pi[i]*(1-pi[i]);
 end;
* w=diag(wvec);

* Fisher=repeat(0,&nvar,&nvar);
* Fisher=t(x)*w*x;
 loglike=0;
 do i=1 to n;
  pi[i]=1/(1+exp(-x[i,]*beta));
  wvec[i]=pi[i]*(1-pi[i]);
  rtwvec[i]=sqrt(wvec[i]);
 end;
 rtwx=repeat(0,n,&srcnvar);
 do i=1 to n;
  rtwx[i,]=x[i,]*rtwvec[i];
 end;
 Fisher=t(rtwx)*rtwx;
 help1=rtwx*inv(Fisher);
 do i=1 to n;
  hvec[i]=help1[i,]*t(rtwx[i,]);
 end;

 do i=1 to n;
*  if (y[i]=1 & pi[i]<&critpi) |(y[i]=0 &(1-pi[i]<&critpi) then do;
*   put "WARNING: pi = 0 or 1 (see output)";
*   pii=pi[i];
*   print i pii beta;
*  end;
  if y[i]=1 then loglike=loglike+log(pi[i]);
  else loglike=loglike+log(1-pi[i]);
 end;
 detfish=det(Fisher);
* if detfish<&critpi then do;
*  put "WARNING: det(Fisher)=0 (see output)";
*  print detfish beta;
* end;
 penlike=loglike+0.5*log(detfish);
 g=t(x)*(y+0.5*hvec-pi-diag(hvec)*pi);
 v=-Fisher;
finish pl_comp;





**********************************************************************;
start newraph;

 * uses:    x, y, ew, os, n, &nvar, &noint, and the subroutine pl_comp;
 * output:  beta, loglike, penlike, vm, U, it;

 dimx=ew[+];
 xwot=repeat(0,n,dimx);
 map=repeat(0,dimx,1);   * to map x to xwot ;

 nvar=&srcnvar;
 jj=0;
 do j=1 to nvar;
  if ew[j]=1 then do;
   jj=jj+1;
   map[jj]=j;
  end;
 end;


 do j=1 to dimx;
  xwot[,j]=x[,map[j]];
 end;


 *print iby start stop n x y;

 beta=os;
 if &noint=0 & ew[1]=1 then do;
  eta=x*beta;
  eta_bar=eta[+]/n;
  y_bar=y[+]/n;
  beta[1]=log(y_bar/(1-y_bar))-eta_bar;
*  beta[1]=-eta_bar;
 end;


 pi=repeat(0.5,n,1);
 wvec=repeat(0.25,n,1);
 rtwvec=repeat(0.5,n,1);
 rtpi=repeat(0.5,n,1);
 rtwxwot=repeat(0,n,ncol(xwot));
 di=repeat(1,&srcnvar,1);
 hvec=repeat(0,n,1);
 diff=1;
 it=0;

 run pl_comp;

* print it beta penlike;

 do while(it<&maxit & diff>abs(&epsilon));
  it=it+1;
  oldbeta=beta;
  do i=1 to n;
   rtwxwot[i,]=xwot[i,]*rtwvec[i];
  end;
  vmwot=inv(t(rtwxwot)*rtwxwot);
  vm=repeat(0,&srcnvar,&srcnvar);
  do j=1 to dimx;
   do jj=1 to dimx;
    vm[map[j],map[jj]]=vmwot[j,jj];
   end;
  end;
  fd=t(x)*(y+0.5*hvec-pi-diag(hvec)*pi);
*  beta=beta+vm*fd;
  do j=1 to &srcnvar;
   incr=0;
   do jj=1 to &srcnvar;
    if ew[jj]=1 then incr=incr+vm[j,jj]*fd[jj];
   end;
   if abs(incr)>&maxstep then incr=sign(incr)*&maxstep;
   if ew[j]=1 then beta[j]=beta[j]+incr;
  end;

  hs=0;
  oldlike=penlike;
  run pl_comp;
  do while(hs<&maxhs & penlike<oldlike);
   hs=hs+1;
   beta=(beta+oldbeta)/2;
   run pl_comp;
  end;

  U=t(g);
  di=beta-oldbeta;
  diff=0;
  do j=1 to &srcnvar;
   diff=diff+abs(di[j]);
  end;

 end;

 run pl_comp;

finish newraph;

**************************************************************************;


use &srcdata;
 read all var("&rby") into by;


%if &noint=0 %then %do;
 read all var("intercep"
  %do srj=1 %to &srck;
   ||"_src&srj"
  %end;
  ) into x_all;
%end;
%else %do;
 read all var("_src1"
  %if &srck>=2 %then %do;
   %do srj=2 %to &srck;
    ||"_src&srj"
   %end;
  %end;
  ) into x_all;
%end;
read all var("&y") into y_all;
close &srcdata;



n_all=nrow(x_all);
maxrby=by[n_all,1];                 
pev=repeat(0,maxrby+1,2);   ********** +1 is because the first dataset is _sample_=0!!! ;

stop=0;
h_index=0;


do iby=0 to maxrby;

 start=stop+1;
 stop=start;

 cond=1;
 do while (by[stop,1]=iby & cond=1);
  stop=stop+1;
  if stop > nrow(by) then do;
   stop=stop-1;
   cond=0;
  end;
 end;

 if cond=1 then stop=stop-1;

 n=stop-start+1;

 x=repeat(0,n,&srck+1-&noint);
 y=repeat(0,n);

 do i=1 to n;
  do k=1 to &srcnvar;
   x[i,k]=x_all[start-1+i,k];
  end;
  y[i]=y_all[start-1+i];
 end;

 ew=repeat(1,&srcnvar,1);
 os=repeat(0,&srcnvar,1);

 %if &shrink=1 %then %do;
  %if &noint=0 %then %do;
   do j=2 to &srcnvar;
    ew[j]=0;
   end;
  %end;
  %if &noint=1 %then %do;
   do j=1 to &srcnvar;
    ew[j]=0;
   end;
  %end;

  run newraph;

  pl_save=penlike;
 %end;

  ew=repeat(1,&srcnvar,1);
  os=repeat(0,&srcnvar,1);

  run newraph;

 %if &shrink=1 %then %do;
  ChiSq=2*(penlike-pl_save);
  shrink=(ChiSq-&srck)/ChiSq;
 %end;
 %else %do;
  shrink=1;
 %end;

 it_save=it;
 
 * computation of R-Square;
 sumpy=0;
 sump=0;
 sumy=0;
 sumpp=0;

 do j=2-&noint to &srcnvar;
  beta[j]=beta[j]*shrink;
 end;
 do i=1 to n;
  xbeta=x[i,]*beta;
  pred=1/(1+exp(-xbeta));
  *print i pred;
  sumpy=sumpy+pred*y[i];
  sump=sump+pred;
  sumy=sumy+y[i];
  sumpp=sumpp+pred*pred;
 end;

 covpy=sumpy/n-sump/n*sumy/n;
 varp=sumpp/n-(sump/n)**2;
 vary=sumy/n-(sumy/n)**2;
 pev[iby+1,1]=iby;
 pev[iby+1,2]=(covpy**2)/varp/vary;
end;

*print covpy varp vary pev;

create &srcout from pev[colname={"&rby","&pevname"}];
append from pev;
close &srcout;

quit;

%end;


options notes source;
*options mprint macrogen;
%end;

data _ttt_;
set _ttt_;
time=time()-x;
time=floor(time*100+0.5)/100;
put "NOTE: Execution of srcore used " time "seconds.";
run;


%mend;

%let clpval=%sysevalf(&clp*100);

%let nvar=0;
%if &varlist ne %then %do;
 %do %while(%scan(&varlist,&nvar+1)~=);
  %let nvar=%eval(&nvar+1);
  %let var&nvar=%scan(&varlist,&nvar);
 %end;
%end;

%let ngroup=0;
%if &groups ne %then %do;
 %do k=1 %to &nvar;
  %let group&k=%scan(&groups,&k);
  %if &&group&k gt &ngroup %then %do;
   %let ngroup=%eval(&&group&k);
  %end;
 %end;
 %do k=1 %to &ngroup;
  %if &gnames ne %then %do;
   %let gname&k=%scan(&gnames,&k);
  %end;
  %else %do;
   %let gname&k=&k;
  %end;
 %end;
%end;

data _xxx;
put "NOTE: &ngroup groups found.";
time=time();
output;
run;

%let ncand=0;
%if &cand ne %then %do;
 %do %while(%scan(&cand,&ncand+1)~=);
  %let ncand=%eval(&ncand+1);
  %let cand&ncand=%scan(&cand,&ncand);
 %end;
%end;

%if &ncand ne 0 %then %do;
 %let mode=CANDIDATE;
 %let nobj=&ncand;
 %do k=1 %to &nobj;
  %let obj&k=&&cand&k;
 %end;
%end;
%else %if &ngroup ne 0 %then %do;
 %let mode=GROUP;
 %let nobj=&ngroup;
 %do k=1 %to &ngroup;
  %let obj&k=&&gname&k;
 %end;
%end;
%else %if &varlist ne  %then %do;
 %let mode=TOTAL;
 %let nobj=&nvar;
  %do k=1 %to &nobj;
  %let obj&k=&&var&k;
 %end;
%end;

%if &nobj=1 %then %do;
 %let nboot=1;
 %let printclp=0;
 %let histogram=0;
%end;

data _worksr;
set &data;
if &y ne .;
%if &varlist ne %then %do;
 %do j=1 %to &nvar;
  if &&var&j ne .;
 %end;
%end;
%if &cand ne %then %do;
 %do j=1 %to &ncand;
  if &&cand&j ne .;
 %end;
%end;
run;

%let nclass=0;

%if &class ne %then %do; 

 %do %while(%scan(&class,&nclass+1)~=);
  %let nclass=%eval(&nclass+1);
  %let class&nclass=%scan(&class,&nclass);
 %end; 


 proc freq data=_worksr;
 %do j=1 %to &nclass;
  tables &&class&j / out=_clt&j noprint;
 %end;
 run;

 proc iml;
 use _worksr;
 read all var("&class1" 
  %if &nclass>1 %then %do;
   %do j=2 %to &nclass;
    || "&&class&j"
   %end;
  %end; ) into x;
 close _worksr;
 newvar=0;
 n=nrow(x);
 df=repeat(0,1,&nclass);
 %do j=1 %to &nclass;
  use _clt&j;
  read all var("&&class&j") into factor&j;
  close _clt&j;
  df[&j]=nrow(factor&j)-1;
  newvar=newvar+df[&j];
 %end;

 newx=repeat(0,n,newvar); 

 do i=1 to n;
  index=0;
  %do j=1 %to &nclass;
   do ii=1 to df[&j];
    index=index+1;
    newx[i,index]=(x[i,&j]=factor&j[ii+1]);
   end;
  %end;
 end;

 create _df from df[colname={%do j=1 %to &nclass; "df&j" %end;}];
 append from df;
 close _df; 

 create _dummies from newx;
 append from newx;
 close _dummies;
 

 quit; 

 data _df;
 set _df;
  %do j=1 %to &nclass;
   call symput("df&j", df&j);
  %end;
 run; 


 data _dummies;
 set _dummies;
  %let index=0;
  %do j=1 %to &nclass;
   %do jj=1 %to &&df&j;
    %let index=%eval(&index+1);
    rename col&index=_dum&j._&jj;
   %end;
  %end;
 run;

 data _worksr;
 merge _worksr _dummies;
 run;

%end;  *class;
 
 *** mapping of factors to dummies ***;

 %if &varlist ne %then %do;
  %do j=1 %to &nvar;
   %let dumv&j=0;
   %do k=1 %to &nclass;
    %if &&var&j=&&class&k %then %do;
	 %let dumv&j=&k;
	%end;
   %end;
  %end;
 %end;
 %if &cand ne %then %do;
  %do j=1 %to &ncand;
   %let dumc&j=0;
   %do k=1 %to &nclass;
    %if &&cand&j=&&class&k %then %do;
	 %let dumc&j=&k;
	%end;
   %end;
  %end;
 %end;


  ********* erzeuge grossen bootstrap-Datensatz **********;
  ********* balanced bootstrap                  **********;

%if &nboot=1 %then %do;
  data _bootsr;
  set _worksr;
   _sample_=0;
  run;
 %end;
 %else %do;
  proc iml; 
  use _worksr; 
  read all var("&y"
   %if &varlist ne %then %do;
    %do k=1 %to &nvar;
	 %if &&dumv&k=0 %then %do;
      ||"&&var&k"
	 %end;
	 %else %do;
	  %let help=&&dumv&k;
	  %do kk=1 %to &&df&help;
	   || "_dum&help._&kk"
	  %end;
	 %end;
    %end; 
   %end;
   %if &cand ne %then %do;
    %do k=1 %to &ncand;
	 %if &&dumc&k = 0 %then %do;
      ||"&&cand&k"
	 %end;
	 %else %do;
	  %let help=&&dumc&k;
	  %do kk=1 %to &&df&help;
	   || "_dum&help._&kk"
	  %end;
     %end; 
    %end; 
   %end;
   ) into data; 
  close _worksr;
  n=nrow(data);
  bootdata=repeat(0,n*&nboot,ncol(data)+1);
  longrow=repeat(0,n*&nboot,1);
  index=1;
  do iboot=1 to &nboot;
   do i=1 to n;
    longrow[index]=i;
	index=index+1;
   end;
  end;
  m=n*&nboot;
  do iboot=1 to &nboot;
   do i=1 to n;
    draw=floor(ranuni(&seed)*m)+1;
    take=longrow[draw];
    longrow[draw]=longrow[m];
	m=m-1;
    bootdata[(iboot-1)*n+i,1]=iboot;
    do j=1 to ncol(data);
     bootdata[(iboot-1)*n+i,j+1]=data[take,j];
    end;
   end;
  end;
  create _bootsr from bootdata[colname={"_sample_","&y"
   %if &varlist ne %then %do;
    %do k=1 %to &nvar;
	 %if &&dumv&k=0 %then %do;
      ,"&&var&k"
	 %end;
	 %else %do;
	  %let help=&&dumv&k;
	  %do kk=1 %to &&df&help;
	   ,"_dum&help._&kk"
	  %end;
	 %end;
    %end; 
   %end;
   %if &cand ne %then %do;
    %do k=1 %to &ncand;
	 %if &&dumc&k = 0 %then %do;
      ,"&&cand&k"
	 %end;
	 %else %do;
	  %let help=&&dumc&k;
	  %do kk=1 %to &&df&help;
	   ,"_dum&help._&kk"
	  %end;
     %end; 
    %end; 
   %end;
   }];
  append from bootdata;
  close _bootsr;
  quit;
  data _worksr;
  set _worksr;
  _sample_=0;
  run;
  data _bootsr;
  set _worksr _bootsr;
  run;
 %end;
  
  
  *************** sorge vor fuer srcore ****************;

  PROC SORT DATA=_bootsr;
  BY _sample_ &y;
  run;

 %let varlist2=; 
 %if &varlist ne %then %do;
  %do j=1 %to &nvar;
   %let jump=0;
   %if &nclass gt 0 %then %do;
    %do jj=1 %to &nclass;
     %if &&class&jj=&&var&j %then %do;
	  %let jump=1;
      %do jjj=1 %to &&df&jj;
       %let varlist2=&varlist2 _dum&jj._&jjj;
      %end;
     %end;
	%end;
   %end;
   %if &jump=0 %then %do;
 	%let varlist2=&varlist2 &&var&j;
   %end;
  %end;
 %end;
 %if &cand ne %then %do;
  %do j=1 %to &ncand;
   %let jump=0;
   %if &nclass gt 0 %then %do;
    %do jj=1 %to &nclass;
     %if &&class&jj=&&cand&j %then %do;
	  %let jump=1;
      %do jjj=1 %to &&df&jj;
       %let varlist2=&varlist2 _dum&jj._&jjj;
      %end;
     %end;
	%end;
   %end;
   %if &jump=0 %then %do;
 	%let varlist2=&varlist2 &&cand&j;
   %end;
  %end;
 %end;

  PROC STANDARD DATA=_bootsr OUT=_bootsr1 MEAN=0;
  VAR &varlist2;
  by _sample_;
  run;

  proc means data=_bootsr1 noprint;
  var &y;
  output out=_ndat n=ndat;
  by _sample_;
  run;
  


%if &nboot=1 %then %do;
 data _ndat;
 set _ndat;
  call symput("n",ndat);
 run;
%end;


%if &cand ne %then %do;

**********************************************************************;
****************** candidate mode ************************************;
**********************************************************************;

  ******************* marginale PEVs *****************************;
  %do j=1 %to &ncand;
   data _take_;
   set _bootsr1;
    %let jump=0;
    %do jj=1 %to &nclass;
     %if &&class&jj=&&cand&j %then %do;
	  %let jump=&&df&jj;
      %do jjj=1 %to &&df&jj;
       _src&jjj=_dum&jj._&jjj;
      %end;
     %end;
	%end;
	%if &jump=0 %then %do;
	 _src1=&&cand&j;
	 %let jump=1;
	%end;
   run;

   %srcore(srcdata=_take_, srck=&jump, rby=_sample_, srcout=_pevmc&j, pevname=_pevmc&j);  

  %end;
  %if &varlist ne %then %do;
   %do j=1 %to &nvar;
    data _take_;
    set _bootsr1;
	if _sample_=0;
     %let jump=0;
     %do jj=1 %to &nclass;
      %if &&class&jj=&&var&j %then %do; 
 	   %let jump=&&df&jj;
       %do jjj=1 %to &&df&jj;
        _src&jjj=_dum&jj._&jjj;
       %end;
      %end;
  	 %end;
	 %if &jump=0 %then %do;
	  _src1=&&var&j;
	  %let jump=1;
	 %end;
    run;

    %srcore(srcdata=_take_, srck=&jump, rby=_sample_, srcout=_pevmv&j, pevname=_pevmv&j);  
   %end;
  %end;



  ******************* partielle PEVs *****************************;
  ******************* Leermodell (nur Variablen aus varlist) *****;
  
  %if &varlist ne %then %do;
   data _take_;
   set _bootsr1;
   if _sample_=0;
    %let srckj=0;
    %do j=1 %to &nvar;
     %let jump=0;
     %do jj=1 %to &nclass;
      %if &&class&jj=&&var&j %then %do;
	   %let jump=1;
       %do jjj=1 %to &&df&jj;
	    %let srckj=%eval(&srckj+1);
        _src&srckj=_dum&jj._&jjj;
       %end;
      %end;
	 %end;
	 %if &jump=0 %then %do;
	   %let srckj=%eval(&srckj+1);
 	  _src&srckj=&&var&j;
     %end;
    %end;
   run;


   %srcore(srcdata=_take_, srck=&srckj, rby=_sample_, srcout=_pevfull, pevname=_pevfull);  

  ************************ Partielle Modelle (mit jeweils einem Kandidaten dazu) ****;
   %do k=1 %to &ncand;
    data _take_;
    set _bootsr1;
     %let srckj=0;    
     %do j=1 %to &nvar;
      %let jump=0;
      %do jj=1 %to &nclass;
       %if &&class&jj=&&var&j %then %do; 
	    %let jump=1;
        %do jjj=1 %to &&df&jj;
	     %let srckj=%eval(&srckj+1);
         _src&srckj=_dum&jj._&jjj;
        %end; *do jjj ;
       %end; *if classjj ;

	  %end; * do jj;
	  %if &jump=0 %then %do;
	   %let srckj=%eval(&srckj+1);
 	   _src&srckj=&&var&j;
	  %end; *if jump;
	 %end; *do j;

     %let jump=0;
     %do jj=1 %to &nclass;
      %if &&class&jj=&&cand&k %then %do; 
       %let jump=1;
       %do jjj=1 %to &&df&jj;
        %let srckj=%eval(&srckj+1);
        _src&srckj=_dum&jj._&jjj;
       %end; *do jjj ;
      %end; *if classjj ;
	 %end; * do jj;
	 %if &jump=0 %then %do;
	  %let srckj=%eval(&srckj+1);
 	  _src&srckj=&&cand&k;
     %end; *if jump;
    run;
    %srcore(srcdata=_take_, srck=&srckj, rby=_sample_, srcout=_pevp&k, pevname=_pevp&k);  
   %end;
  %end;



********************                           *****************************;
  data _pev_;
  merge %if &varlist ne %then %do; _pevfull  %end;
   %do j=1 %to &ncand;
    _pevmc&j %if &varlist ne %then %do; _pevp&j %end;
   %end; 
   %if &varlist ne %then %do;
    %do j=1 %to &nvar;
     _pevmv&j
	%end;
   %end; ;
   by _sample_;
   %do j=1 %to &ncand;
    %if &varlist ne %then %do;
	 _asfull=arsin(sqrt(_pevfull));
     if _sample_=0 then do;
      _pevp&j=max(_pevp&j-_pevfull,0);
	 end;
	 _asp&j=arsin(sqrt(_pevp&j));
	%end;
    _asmc&j=arsin(sqrt(_pevmc&j));
   %end;
   %if &varlist ne %then %do;
    %do j=1 %to &nvar;
	 _asmv&j=arsin(sqrt(_pevmv&j));
	%end;
   %end;
   %do j=1 %to &ncand-1;
    %do k=&j+1 %to &ncand;
     %if &varlist ne %then %do;
      dp&j._&k=_asp&j-_asp&k;
	 %end;
	 dm&j._&k=_asmc&j-_asmc&k;
	%end;
   %end;
  run;

  proc univariate noprint data=_pev_;
  where _sample_>0;
  var %if &varlist ne %then %do;
/*       _asfull   */
	  %end;
      %do j=1 %to &ncand; 
       %if &varlist ne %then %do;
        _asp&j 
	   %end;
       _asmc&j 
      %end;
	  %if &varlist ne %then %do;
	   %do j=1 %to &nvar;
	    _asmv&j
	   %end;
	  %end;
   %do j=1 %to &ncand-1; 
    %do k=&j+1 %to &ncand;
     %if &varlist ne %then %do;
      dp&j._&k 
     %end;
     dm&j._&k
	%end;
   %end; ;
  output out=_relimp_ 

   std=    
   %if &varlist ne %then %do;
/*    sasfull   */
   %end;
   %do j=1 %to &ncand;
    %if &varlist ne %then %do;
	 sasp&j
	%end;
	sasmc&j
   %end;
   %if &varlist ne %then %do;
    %do j=1 %to &nvar;
     sasmv&j 
	%end;
   %end;
   %do j=1 %to &ncand-1; 
    %do k=&j+1 %to &ncand;
	 %if &varlist ne %then %do;
 	  sdp&j._&k 
     %end;
     sdm&j._&k
	%end;
   %end; 
  ;
  run;

 
  data _srcbbm_;
  set _pev_;
  if _sample_=0;
   %if &varlist ne %then %do;
    masfull=_asfull;
   %end;
   %do j=1 %to &ncand;
    %if &varlist ne %then %do;
	 masp&j=_asp&j;
	%end;
	masmc&j=_asmc&j;
   %end;
   %if &varlist ne %then %do;
    %do j=1 %to &nvar;
     masmv&j =_asmv&j;
	%end;
   %end;
   %do j=1 %to &ncand-1; 
    %do k=&j+1 %to &ncand;
	 %if &varlist ne %then %do;
 	  mdp&j._&k=dp&j._&k;
     %end;
     mdm&j._&k=dm&j._&k;
	%end;
   %end; 
   run;

  data _relimp_;
  merge _srcbbm_ _relimp_;
  file print;
   length name $ 8;
   %if &nboot>1 %then %do;
    _chia=0.5*(probit(1-(1-&clp)/2)+sqrt(2*&nboot-3))**2;
	_chi1ma=0.5*(-probit(1-(1-&clp)/2)+sqrt(2*&nboot-3))**2;
    %do j=1 %to &ncand-1;
     %do k=&j+1 %to &ncand;
	  sdm&j._&k.up=sqrt((&nboot-1)*sdm&j._&k**2/_chia);
	  sdm&j._&k.lo=sqrt((&nboot-1)*sdm&j._&k**2/_chi1ma);
      tm&j._&k=abs(mdm&j._&k)/(sdm&j._&k);
      tm&j._&k.lo=abs(mdm&j._&k)/(sdm&j._&k.up);
      tm&j._&k.up=abs(mdm&j._&k)/(sdm&j._&k.lo);
      %if &varlist ne %then %do;
	   sdp&j._&k.up=sqrt((&nboot-1)*sdp&j._&k**2/_chia);
	   sdp&j._&k.lo=sqrt((&nboot-1)*sdp&j._&k**2/_chi1ma);
       tp&j._&k=abs(mdp&j._&k)/(sdp&j._&k);
       tp&j._&k.lo=abs(mdp&j._&k)/(sdp&j._&k.up);
       tp&j._&k.up=abs(mdp&j._&k)/(sdp&j._&k.lo);
      %end;

      %if %eval(&ncand*&correctp)>2 %then %do;
       pm&j._&k=floor((1-probmc('RANGE',tm&j._&k,.,.,&ncand))*10000+0.5)/10000;
	   pm&j._&k.lo=floor((1-probmc('RANGE',tm&j._&k.lo,.,.,&ncand))*10000+0.5)/10000;
	   pm&j._&k.up=floor((1-probmc('RANGE',tm&j._&k,up,.,.,&ncand))*10000+0.5)/10000;
       %if &varlist ne %then %do;
           pp&j._&k=floor((1-probmc('RANGE',tp&j._&k,.,.,&ncand))*10000+0.5)/10000;
           pp&j._&k.lo=floor((1-probmc('RANGE',tp&j._&k.lo,.,.,&ncand))*10000+0.5)/10000;
           pp&j._&k.up=floor((1-probmc('RANGE',tp&j._&k.up,.,.,&ncand))*10000+0.5)/10000;
       %end;
      %end;
      %else %do;
       pm&j._&k=floor((1-probnorm(tm&j._&k))*2*10000+0.5)/10000;
       pm&j._&k.lo=floor((1-probnorm(tm&j._&k.lo))*2*10000+0.5)/10000;
       pm&j._&k.up=floor((1-probnorm(tm&j._&k.up))*2*10000+0.5)/10000;
       %if &varlist ne %then %do;
        pp&j._&k=floor((1-probnorm(tp&j._&k))*2*10000+0.5)/10000;
        pp&j._&k.lo=floor((1-probnorm(tp&j._&k.lo))*2*10000+0.5)/10000;
        pp&j._&k.up=floor((1-probnorm(tp&j._&k.up))*2*10000+0.5)/10000;
       %end;
      %end;
        label %if &varlist ne %then %do;
             pp&j._&k="Partial P-value for &&cand&j &&cand&k"
                %end;
            pm&j._&k="Marginal P-value for &&cand&j &&cand&k";
     %end;
    %end; ;
    hsd=probmc('RANGE',.,0.95,.,3);
   %end;
   %if &varlist ne %then %do;
    mpfull=floor(10000*sin(masfull)**2+0.5)/100;
	label mpfull="PEV for model without candidates";
	%do j=1 %to &nvar;
	 mmv&j=floor(10000*sin(masmv&j)**2+0.5)/100;
	 label mmv&j="Marginal PEV for &&var&j";
	%end;
   %end;
   %do j=1 %to &ncand;
    %if &varlist ne %then %do;
     mpc&j=floor(10000*sin(masp&j)**2+0.5)/100;
	%end;
    mmc&j=floor(10000*sin(masmc&j)**2+0.5)/100;
	label    
     %if &varlist ne %then %do;
      mpc&j="Partial PEV for &&cand&j"
	 %end;
     mmc&j="Marginal PEV for &&cand&j";
   %end;
   put "Proportion of Explained Variation (PEV)";
   put "=======================================";
   put " ";
   %if &firth=1 %then %do;
    put "Parameter estimates based on penalized maximum likelihood.";
	put " ";
   %end;
   %if &shrink=1 %then %do;
    put "Parameter estimates have been shrinked.";
	put " ";
   %end;
   %if &varlist ne %then %do;
    put "PEV         Marginal   Partial";
    put "========    ========   =======";
    put " ";
   %end;
   %else %do;
    put "PEV         Marginal";
	put "========    ========";
	put " ";
   %end;
   %if &varlist ne %then %do;
    %do j=1 %to &nvar;
     name="&&var&j"; 
 	 perc="% ";
     put name 1-12 mmv&j F7.2 perc;
    %end;
    put "-------------------------------";
	put "Model       " mpfull F7.2 perc;
	put " ";
	put "Candidates:";
	put "===========";
	put " ";
   %end;
   %do j=1 %to &ncand;
    name="&&cand&j";
	perc="% ";
	put name 1-12 mmc&j f7.2 perc @@;
    %if &varlist ne %then %do;
     put mpc&j f8.2 perc @@;
	%end;
	put;
   %end;
   %if &nboot>1 %then %do;
    put " ";
	put " ";
    put "All comparisons based on &nboot bootstrap replicates.";
   
    put " ";
    put " ";
    put "Comparison (" @@;
    %if %eval(&ncand*&correctp)>2 %then %do;
     put "Tukey-HSD-corrected " @@;
    %end;
    put "p-values): marginal PEV";
    %if %eval(&ncand*&correctp)>2 %then %do;
     put "====================" @@;
    %end;
    put "===================================" ;
    put " ";
    %do k=1 %to &ncand; 
     %let stcol=%eval(&k*9+1);
     %let encol=%eval(&stcol+8);  
 	 name="&&cand&k";
     put name &stcol-&encol @@;
    %end; 
    put;
    put " ";
    %do k=1 %to &ncand;
     name="&&cand&k";
     put name 1-8 @@;
 	 %let stcol=10;
     %do j=1 %to &ncand;
      %let encol=%eval(&stcol+8);
	  %if &j le &k %then %do; 
	   name=" .";
       put name &stcol-&encol @@;
	  %end;
	  %else %do;
         put pm&k._&j pvalue6. @@;
         if pm&k._&j.lo<&alpha and pm&k._&j.up>&alpha then do; 
          put "*  " @@;
         end;
         else do;
          put "   " @@;
		 end;
	  %end;
	  %let stcol=%eval(&stcol+9);
	 %end;
	 put;
    %end;
    %if &varlist ne %then %do;
     put " ";
     put " ";
     put "Comparison (" @@;
     %if %eval(&ncand*&correctp)>2 %then %do;
      put "Tukey-HSD-corrected " @@;
     %end;
     put "p-values): partial PEV";
     %if %eval(&ncand*&correctp)>2 %then %do;
      put "====================" @@;
     %end;
     put "==================================" ;
     put " ";
     %do k=1 %to &ncand; 
      name="&&cand&k";
      %let stcol=%eval(&k*9+1);
      %let encol=%eval(&stcol+8); 
      put name &stcol-&encol @@;
     %end; 
     put;
     put " ";
     %do k=1 %to &ncand;
      name="&&cand&k";
      put name 1-8 @@; 

 	  %let stcol=10;
      %do j=1 %to &ncand;
       %let encol=%eval(&stcol+8);
	   %if &j le &k %then %do; 
	    name=" .";
        put name &stcol-&encol @@;
	   %end;
	   %else %do;
         put pp&k._&j pvalue6. @@;
         if pp&k._&j.lo<&alpha and pp&k._&j.up>&alpha then do; 
          put "*  " @@;
         end;
         else do;
          put "   " @@;
		 end;
	   %end;
	   %let stcol=%eval(&stcol+9);
	  %end;
	  put;
	 %end;
    %end;
    put "   ";
    put "*: &clpval.%-C.L. for p-value contains the critical value &alpha..";
   %end;
  run;

%end;

****************************************************************************;
********************** total mode ******************************************;
********************** no groups  ******************************************;
****************************************************************************;

%if &ngroup=0 %then %do;
%if &varlist ne %then %do;
 %if &cand = %then %do;

 *****************************  partielle und marginale PEVs, alle mit allen vergleichen ***;
  ******************* marginale PEVs *****************************;
  %do j=1 %to &nvar;
   data _take_;
   set _bootsr1;
    %let jump=0;
    %do jj=1 %to &nclass;
     %if &&class&jj=&&var&j %then %do;
	  %let jump=&&df&jj;
      %do jjj=1 %to &&df&jj;
       _src&jjj=_dum&jj._&jjj;
      %end;
     %end;
	%end;
	%if &jump=0 %then %do;
	 _src1=&&var&j;
	 %let jump=1;
	%end;
   run;

   %srcore(srcdata=_take_, srck=&jump, rby=_sample_, srcout=_pevm&j, pevname=_pevm&j);  

  %end;
  ******************* partielle PEVs *****************************;
  ******************* Vollmodell      *****************************;
   data _take_;
   set _bootsr1;
    %let srckj=0;
    %do j=1 %to &nvar;
	 %let jump=0;
     %do jj=1 %to &nclass;
      %if &&class&jj=&&var&j %then %do;
	   %let jump=1;
       %do jjj=1 %to &&df&jj;
	    %let srckj=%eval(&srckj+1);
        _src&srckj=_dum&jj._&jjj;
       %end;
      %end;
	 %end;
	 %if &jump=0 %then %do;
	   %let srckj=%eval(&srckj+1);
 	  _src&srckj=&&var&j;
     %end;
    %end;

   run;

   %srcore(srcdata=_take_, srck=&srckj, rby=_sample_, srcout=_pevfull, pevname=_pevfull);  

  ************************ Partielle Modelle *************************************;
  %do k=1 %to &nvar;
   data _take_;
   set _bootsr1;
    %let srckj=0;    
    %do j=1 %to &nvar;
	 %if &j ne &k %then %do;
      %let jump=0;
      %do jj=1 %to &nclass;
       %if &&class&jj=&&var&j %then %do;
	    %let jump=1;
        %do jjj=1 %to &&df&jj;
	     %let srckj=%eval(&srckj+1);
         _src&srckj=_dum&jj._&jjj;
        %end;
       %end;
	  %end;
	  %if &jump=0 %then %do;
	   %let srckj=%eval(&srckj+1);
 	   _src&srckj=&&var&j;
	  %end;
	 %end;
	%end;
   run;

   %srcore(srcdata=_take_, srck=&srckj, rby=_sample_, srcout=_pevp&k, pevname=_pevp&k);  

  %end;
  
  data _pev_;
  merge _pevfull 
   %do j=1 %to &nvar;
    _pevm&j _pevp&j
   %end; ;
   by _sample_;
   %do j=1 %to &nvar;
*    _pevp&j=max(_pevfull-_pevp&j,0);
	_asp&j=arsin(sqrt(_pevp&j));
	_asm&j=arsin(sqrt(_pevm&j));
   %end;
   %do j=1 %to &nvar-1;
    %do k=&j+1 %to &nvar;
	 dp&j._&k=_asp&k-_asp&j;
	 dm&j._&k=_asm&j-_asm&k;
	%end;
   %end;
   %do j=1 %to &nvar;
    _pevp&j=max(_pevfull-_pevp&j,0);
	_asp&j=arsin(sqrt(_pevp&j));
	_asm&j=arsin(sqrt(_pevm&j));
   %end;
  run;

  proc univariate noprint data=_pev_;
  where _sample_>0;
  var %do j=1 %to &nvar; _asp&j _asm&j %end;
   %do j=1 %to &nvar-1; 
    %do k=&j+1 %to &nvar;
	 dp&j._&k dm&j._&k
	%end;
   %end; ;
  output out=_relimp_ 
   std=   %do j=1 %to &nvar; sasp&j sasm&j %end;
   %do j=1 %to &nvar-1; 
    %do k=&j+1 %to &nvar;
	 sdp&j._&k sdm&j._&k
	%end;
   %end; ;
  run;

  data _srcbbm_;
  set _pev_;
  if _sample_=0;
   %do j=1 %to &nvar; 
    masp&j=_asp&j;
    masm&j=_asm&j;
   %end;
   %do j=1 %to &nvar-1; 
    %do k=&j+1 %to &nvar;
	 mdp&j._&k=dp&j._&k;
     mdm&j._&k=dm&j._&k;
	%end;
   %end; 
  run;


  data _relimp_;
  merge _srcbbm_ _relimp_;
  file print;
   length name $ 8;
   %if &nboot>1 %then %do;
    _chia=0.5*(probit(1-(1-&clp)/2)+sqrt(2*&nboot-3))**2;
	_chi1ma=0.5*(-probit(1-(1-&clp)/2)+sqrt(2*&nboot-3))**2;
    %do j=1 %to &nvar-1;
     %do k=&j+1 %to &nvar;
	  sdm&j._&k.up=sqrt((&nboot-1)*sdm&j._&k**2/_chia);
	  sdm&j._&k.lo=sqrt((&nboot-1)*sdm&j._&k**2/_chi1ma);
	  sdp&j._&k.up=sqrt((&nboot-1)*sdp&j._&k**2/_chia);
	  sdp&j._&k.lo=sqrt((&nboot-1)*sdp&j._&k**2/_chi1ma);

      tm&j._&k=abs(mdm&j._&k)/(sdm&j._&k*(&nboot-1)/(&nboot));
      tp&j._&k=abs(mdp&j._&k)/(sdp&j._&k*(&nboot-1)/(&nboot));
      tm&j._&k.lo=abs(mdm&j._&k)/(sdm&j._&k.lo*(&nboot-1)/(&nboot));
      tp&j._&k.lo=abs(mdp&j._&k)/(sdp&j._&k.lo*(&nboot-1)/(&nboot));
      tm&j._&k.up=abs(mdm&j._&k)/(sdm&j._&k.up*(&nboot-1)/(&nboot));
      tp&j._&k.up=abs(mdp&j._&k)/(sdp&j._&k.up*(&nboot-1)/(&nboot));
      %if %eval(&nvar*&correctp)>3 %then %do;
       pm&j._&k=floor((1-probmc('RANGE',tm&j._&k,.,.,&nvar))*10000+0.5)/10000;
       pp&j._&k=floor((1-probmc('RANGE',tp&j._&k,.,.,&nvar))*10000+0.5)/10000;
       pm&j._&k.up=floor((1-probmc('RANGE',tm&j._&k.lo,.,.,&nvar))*10000+0.5)/10000;
       pp&j._&k.up=floor((1-probmc('RANGE',tp&j._&k.lo,.,.,&nvar))*10000+0.5)/10000;
       pm&j._&k.lo=floor((1-probmc('RANGE',tm&j._&k.up,.,.,&nvar))*10000+0.5)/10000;
       pp&j._&k.lo=floor((1-probmc('RANGE',tp&j._&k.up,.,.,&nvar))*10000+0.5)/10000;
      %end;
      %else %do;
       pm&j._&k=floor((1-probnorm(tm&j._&k))*2*10000+0.5)/10000;
       pp&j._&k=floor((1-probnorm(tp&j._&k))*2*10000+0.5)/10000;
       pm&j._&k.up=floor((1-probnorm(tm&j._&k.lo))*2*10000+0.5)/10000;
       pp&j._&k.up=floor((1-probnorm(tp&j._&k.lo))*2*10000+0.5)/10000;
       pm&j._&k.lo=floor((1-probnorm(tm&j._&k.up))*2*10000+0.5)/10000;
       pp&j._&k.lo=floor((1-probnorm(tp&j._&k.up))*2*10000+0.5)/10000;
      %end;
      label pp&j._&k="Partial P-value for &&var&j &&var&k"
            pm&j._&k="Marginal P-value for &&var&j &&var&k";
     %end;
    %end; ;
    hsd=probmc('RANGE',.,0.95,.,3);
   %end;
   %do j=1 %to &nvar;
    mp&j=floor(10000*sin(masp&j)**2+0.5)/100;
    mm&j=floor(10000*sin(masm&j)**2+0.5)/100;
	label mp&j="Partial PEV for &&var&j"
          mm&j="Marginal PEV for &&var&j";
   %end;
   put "Proportion of Explained Variation (PEV)";
   put "=======================================";
   put " ";
   %if &firth=1 %then %do;
    put "Parameter estimates based on penalized maximum likelihood.";
	put " ";
   %end;
   %if &shrink=1 %then %do;
    put "Parameter estimates have been shrinked.";
	put " ";
   %end;
   put "PEV         Marginal   Partial";
   put "========    ========   =======";
   put " ";
   %do j=1 %to &nvar;
    name="&&var&j";
	perc="% ";
    put name 1-12 mm&j F7.2 perc mp&j F8.2 perc;
   %end;
    put "-------------------------------";
	_pevfull=_pevfull*100;
	put "Model       " _pevfull F7.2 perc;
    put " ";
    put " ";
   %if &nboot>1 %then %do;
    put "All comparisons based on &nboot bootstrap replicates.";
	put " ";
    put "Comparison (" @@;
    %if %eval(&nvar*&correctp)>3 %then %do;
     put "Tukey-HSD-corrected " @@;
    %end;
    put "p-values): marginal PEV";
    %if %eval(&nvar*&correctp)>3 %then %do;
     put "====================" @@;
    %end;
    put "===================================" ;
    put " ";
    %do k=1 %to &nvar;
     %let stcol=%eval(&k*9+1);
     %let encol=%eval(&stcol+8);
       name="&&var&k";
     put name &stcol-&encol @@;
    %end;
    put;
    put " ";
    %do k=1 %to &nvar;
     name="&&var&k";
     put name 1-8 @@;
       %let stcol=10;
     %do j=1 %to &nvar;
      %let encol=%eval(&stcol+8);
        %if &j le &k %then %do;
         name=" .";
       put name &stcol-&encol @@;
        %end;
        %else %do;
         put pm&k._&j pvalue6. @@;
         if pm&k._&j.lo<&alpha and pm&k._&j.up>&alpha then do; 
          put "*  " @@;
         end;
         else do;
          put "   " @@;
		 end;
        %end;
        %let stcol=%eval(&stcol+9);
       %end;
       put;
    %end;
    put " ";
    put " ";
    put "Comparison (" @@;
    %if %eval(&nvar*&correctp)>3 %then %do;
     put "Tukey-HSD-corrected " @@;
    %end;
    put "p-values): partial PEV";

    %if %eval(&nvar*&correctp)>3 %then %do;
     put "====================" @@;
    %end;
    put "==================================" ;
    put " ";
    %do k=1 %to &nvar;
     name="&&var&k";
     %let stcol=%eval(&k*9+1);
     %let encol=%eval(&stcol+8);
     put name &stcol-&encol @@;
    %end;
    put;
    put " ";
    %do k=1 %to &nvar;
     name="&&var&k";
     put name 1-8 @@;

       %let stcol=10;
     %do j=1 %to &nvar;
      %let encol=%eval(&stcol+8);
        %if &j le &k %then %do;
         name=" .";
       put name &stcol-&encol @@;
        %end;
        %else %do;
         put pp&k._&j pvalue6. @@;
         if pp&k._&j.lo<&alpha and pp&k._&j.up>&alpha then do; 
          put "*  " @@;
         end;
         else do;
          put "   " @@;
		 end;
        %end;
        %let stcol=%eval(&stcol+9);
       %end;
       put;
    %end;
	put "   ";
	put "*: &clpval.%-C.L. for p-value contains the critical value &alpha..";
   %end;
   run;

 %end;
%end;
%end;

****************************************************************************;
****************************************************************************;
********************** total mode ******************************************;
********************** groups     ******************************************;
****************************************************************************;

%if &ngroup>0 %then %do;
 %if &varlist ne %then %do;
  %if &cand = %then %do;

  *****************************  partielle und marginale PEVs, alle Gruppen mit allen vergleichen ***;
  ******************* marginale PEVs *****************************;
  %do g=1 %to &ngroup;
   data _take_;
   set _bootsr1;
    %let jump=0;
    %do j=1 %to &nvar;
	 %if &&group&j=&g %then %do;
	  %let classok=0;
      %do jj=1 %to &nclass;
       %if &&class&jj=&&var&j %then %do; 
	    %let classok=1;
        %do jjj=1 %to &&df&jj;
		 %let jump=%eval(&jump+1);
         _src&jump=_dum&jj._&jjj;
        %end; *jjj;
       %end; *if;
	  %end; *jj;
	  %if &classok=0 %then %do;
	   %let jump=%eval(&jump+1);
	   _src&jump=&&var&j;
	  %end; *if classok;
	 %end; *if groupj=g;
	%end; *j;
   run;

   %srcore(srcdata=_take_, srck=&jump, rby=_sample_, srcout=_pevm&g, pevname=_pevm&g);  

  %end;
  ******************* partielle PEVs *****************************;
  ******************* Vollmodell      *****************************;
   data _take_;
   set _bootsr1;
    %let srckj=0;
    %do j=1 %to &nvar;
	 %let jump=0;
     %do jj=1 %to &nclass;
      %if &&class&jj=&&var&j %then %do;
	   %let jump=1;
       %do jjj=1 %to &&df&jj;
	    %let srckj=%eval(&srckj+1);
        _src&srckj=_dum&jj._&jjj;
       %end;
      %end;
	 %end;
	 %if &jump=0 %then %do;
	   %let srckj=%eval(&srckj+1);
 	  _src&srckj=&&var&j;
     %end;
    %end;

   run;

   %srcore(srcdata=_take_, srck=&srckj, rby=_sample_, srcout=_pevfull, pevname=_pevfull);  

  ************************ Partielle Modelle *************************************;
  %do k=1 %to &ngroup;
   data _take_;
   set _bootsr1;
    %let srckj=0;    
    %do j=1 %to &nvar;
	 %if &&group&j ne &k %then %do;
	  %let classok=0;
      %do jj=1 %to &nclass;
       %if &&class&jj=&&var&j %then %do;
	    %let classok=1;
        %do jjj=1 %to &&df&jj;
	     %let srckj=%eval(&srckj+1);
         _src&srckj=_dum&jj._&jjj;
        %end;
       %end;
	  %end;
	  %if &classok=0 %then %do;
      %let srckj=%eval(&srckj+1);
       _src&srckj=&&var&j;
	  %end;
	 %end;
    %end;
   run;

   %srcore(srcdata=_take_, srck=&srckj, rby=_sample_, srcout=_pevp&k, pevname=_pevp&k);  

  %end;
  
  data _pev_;
  merge _pevfull 
   %do j=1 %to &ngroup;
    _pevm&j _pevp&j
   %end; ;
   by _sample_;
   %do j=1 %to &ngroup;
*    _pevp&j=max(_pevfull-_pevp&j,0);
	_asp&j=arsin(sqrt(_pevp&j));
	_asm&j=arsin(sqrt(_pevm&j));
   %end;
   %do j=1 %to &ngroup-1;
    %do k=&j+1 %to &ngroup;
	 dp&j._&k=_asp&k-_asp&j;
	 dm&j._&k=_asm&j-_asm&k;
	%end;
   %end;
   %do j=1 %to &ngroup;
    _pevp&j=max(_pevfull-_pevp&j,0);
	_asp&j=arsin(sqrt(_pevp&j));
	_asm&j=arsin(sqrt(_pevm&j));
   %end;
  run;

  proc univariate noprint data=_pev_;
  where _sample_>0;
  var %do j=1 %to &ngroup; _asp&j _asm&j %end;
   %do j=1 %to &ngroup-1; 
    %do k=&j+1 %to &ngroup;
	 dp&j._&k dm&j._&k
	%end;
   %end; ;
  output out=_relimp_ 
   std=   %do j=1 %to &ngroup; sasp&j sasm&j %end;
   %do j=1 %to &ngroup-1; 
    %do k=&j+1 %to &ngroup;
	 sdp&j._&k sdm&j._&k
	%end;
   %end; ;
  run;

  data _srcbbm_;
  set _pev_;
  if _sample_=0;
   %do j=1 %to &ngroup; 
    masp&j=_asp&j;
    masm&j=_asm&j;
   %end;
   %do j=1 %to &ngroup-1; 
    %do k=&j+1 %to &ngroup;
	 mdp&j._&k=dp&j._&k;
     mdm&j._&k=dm&j._&k;
	%end;
   %end; 
  run;


  data _relimp_;
  merge _srcbbm_ _relimp_;
  file print;
   length name $ 8 vars $ %eval(9*&nvar);
   %if &nboot>1 %then %do;
    _chia=0.5*(probit(1-(1-&clp)/2)+sqrt(2*&nboot-3))**2;
	_chi1ma=0.5*(-probit(1-(1-&clp)/2)+sqrt(2*&nboot-3))**2;
    %do j=1 %to &ngroup-1;
     %do k=&j+1 %to &ngroup;
	  sdm&j._&k.up=sqrt((&nboot-1)*sdm&j._&k**2/_chia);
	  sdm&j._&k.lo=sqrt((&nboot-1)*sdm&j._&k**2/_chi1ma);
	  sdp&j._&k.up=sqrt((&nboot-1)*sdp&j._&k**2/_chia);
	  sdp&j._&k.lo=sqrt((&nboot-1)*sdp&j._&k**2/_chi1ma);

      tm&j._&k=abs(mdm&j._&k)/(sdm&j._&k*(&nboot-1)/(&nboot));
      tp&j._&k=abs(mdp&j._&k)/(sdp&j._&k*(&nboot-1)/(&nboot));
      tm&j._&k.lo=abs(mdm&j._&k)/(sdm&j._&k.lo*(&nboot-1)/(&nboot));
      tp&j._&k.lo=abs(mdp&j._&k)/(sdp&j._&k.lo*(&nboot-1)/(&nboot));
      tm&j._&k.up=abs(mdm&j._&k)/(sdm&j._&k.up*(&nboot-1)/(&nboot));
      tp&j._&k.up=abs(mdp&j._&k)/(sdp&j._&k.up*(&nboot-1)/(&nboot));
      %if %eval(&ngroup*&correctp)>2 %then %do;
       pm&j._&k=floor((1-probmc('RANGE',tm&j._&k,.,.,&ngroup))*10000+0.5)/10000;
       pp&j._&k=floor((1-probmc('RANGE',tp&j._&k,.,.,&ngroup))*10000+0.5)/10000;
       pm&j._&k.up=floor((1-probmc('RANGE',tm&j._&k.lo,.,.,&ngroup))*10000+0.5)/10000;
       pp&j._&k.up=floor((1-probmc('RANGE',tp&j._&k.lo,.,.,&ngroup))*10000+0.5)/10000;
       pm&j._&k.lo=floor((1-probmc('RANGE',tm&j._&k.up,.,.,&ngroup))*10000+0.5)/10000;
       pp&j._&k.lo=floor((1-probmc('RANGE',tp&j._&k.up,.,.,&ngroup))*10000+0.5)/10000;
      %end;
      %else %do;
       pm&j._&k=floor((1-probnorm(tm&j._&k))*2*10000+0.5)/10000;
       pp&j._&k=floor((1-probnorm(tp&j._&k))*2*10000+0.5)/10000;
       pm&j._&k.up=floor((1-probnorm(tm&j._&k.lo))*2*10000+0.5)/10000;
       pp&j._&k.up=floor((1-probnorm(tp&j._&k.lo))*2*10000+0.5)/10000;
       pm&j._&k.lo=floor((1-probnorm(tm&j._&k.up))*2*10000+0.5)/10000;
       pp&j._&k.lo=floor((1-probnorm(tp&j._&k.up))*2*10000+0.5)/10000;
      %end;
      label pp&j._&k="Partial P-value for &&gname&j &&gname&k"
            pm&j._&k="Marginal P-value for &&gname&j &&gname&k";
     %end;
    %end; ;
    hsd=probmc('RANGE',.,0.95,.,3);
   %end;
   %do j=1 %to &ngroup;
    mp&j=floor(10000*sin(masp&j)**2+0.5)/100;
    mm&j=floor(10000*sin(masm&j)**2+0.5)/100;
	label mp&j="Partial PEV for &&gname&j"
          mm&j="Marginal PEV for &&gname&j";
   %end;
   put "Proportion of Explained Variation (PEV)";
   put "=======================================";
   put " ";
   %if &firth=1 %then %do;
    put "Parameter estimates based on penalized maximum likelihood.";
	put " ";
   %end;
   %if &shrink=1 %then %do;
    put "Parameter estimates have been shrinked.";
	put " ";
   %end;
   put "PEV         Marginal   Partial    Variables";
   put "========    ========   =======    =========";
   put " ";
   %do j=1 %to &ngroup;
    %let vars=;
    %do k=1 %to &nvar; 
     %if &&group&k=&j %then %do; 
      %let vars=&vars &&var&k ; 
     %end; 
    %end;
    vars="&vars";
    name="&&gname&j";
	perc="% ";
    put name 1-12 mm&j F7.2 perc mp&j F8.2 perc "   " vars;
   %end;
    put "-------------------------------";
	_pevfull=_pevfull*100;
	put "Model       " _pevfull F7.2 perc;
    put " ";
    put " ";
   %if &nboot>1 %then %do;
    put "All comparisons based on &nboot bootstrap replicates.";
	put " ";
    put "Comparison (" @@;
    %if %eval(&ngroup*&correctp)>2 %then %do;
     put "Tukey-HSD-corrected " @@;
    %end;
    put "p-values): marginal PEV";
    %if %eval(&ngroup*&correctp)>2 %then %do;
     put "====================" @@;
    %end;
    put "===================================" ;
    put " ";
    %do k=1 %to &ngroup;
     %let stcol=%eval(&k*9+1);
     %let encol=%eval(&stcol+8);
       name="&&gname&k";
     put name &stcol-&encol @@;
    %end;
    put;
    put " ";
    %do k=1 %to &ngroup;
     name="&&gname&k";
     put name 1-8 @@;
       %let stcol=10;
     %do j=1 %to &ngroup;
      %let encol=%eval(&stcol+8);
        %if &j le &k %then %do;
         name=" .";
       put name &stcol-&encol @@;
        %end;
        %else %do;
         put pm&k._&j pvalue6. @@;
         if pm&k._&j.lo<&alpha and pm&k._&j.up>&alpha then do; 
          put "*  " @@;
         end;
         else do;
          put "   " @@;
		 end;
        %end;
        %let stcol=%eval(&stcol+9);
       %end;
       put;
    %end;
    put " ";
    put " ";
    put "Comparison (" @@;
    %if %eval(&ngroup*&correctp)>2 %then %do;
     put "Tukey-HSD-corrected " @@;
    %end;
    put "p-values): partial PEV";

    %if %eval(&ngroup*&correctp)>2 %then %do;
     put "====================" @@;
    %end;
    put "==================================" ;
    put " ";
    %do k=1 %to &ngroup;
     name="&&gname&k";
     %let stcol=%eval(&k*9+1);
     %let encol=%eval(&stcol+8);
     put name &stcol-&encol @@;
    %end;
    put;
    put " ";
    %do k=1 %to &ngroup;
     name="&&gname&k";
     put name 1-8 @@;

       %let stcol=10;
     %do j=1 %to &ngroup;
      %let encol=%eval(&stcol+8);
        %if &j le &k %then %do;
         name=" .";
       put name &stcol-&encol @@;
        %end;
        %else %do;
         put pp&k._&j pvalue6. @@;
         if pp&k._&j.lo<&alpha and pp&k._&j.up>&alpha then do; 
          put "*  " @@;
         end;
         else do;
          put "   " @@;
		 end;
        %end;
        %let stcol=%eval(&stcol+9);
       %end;
       put;
    %end;
	put "   ";
	put "*: &clpval.%-C.L. for p-value contains the critical value &alpha..";
   %end;
   run;

 %end;
%end;
%end;

%if &nboot*&printclp>1 %then %do;
 data _relimp_;
 set _relimp_;
 file print;
 length na1 $ 8 na2 $ 8;
 put "&clpval.%-confidence intervals for P-values";
 put       "=====================================";
 put " ";
 put "                        Marginal PEV           Partial PEV";
 put "                      =================     =================";
 put " ";
 %if &ngroup>0 %then %let bis=&ngroup;
 %else %if &ncand>1 %then %let bis=&ncand;
 %else %if &varlist ne %then %let bis=&nvar;
 %do j=1 %to %eval(&bis-1);
  %do k=%eval(1+&j) %to &bis;
   %if &ngroup>0 %then %do; 
    na1="&&gname&j";
	na2="&&gname&k";
   %end;
   %else %if &ncand>1 %then %do;
    na1="&&cand&j";
	na2="&&cand&k";
   %end;
   %else %if &varlist ne %then %do;
    na1="&&var&j";
	na2="&&var&k";
   %end;
   put na1 1-8 na2 10-17 "     [" pm&&j._&k.lo pvalue6. " , " pm&&j._&k.up pvalue6. "]     [" pp&&j._&k.lo pvalue6. " , " pp&&j._&k.up pvalue6. "]";
  %end;
 %end;
 run;
%end;


%if &histogram=1 %then %do;

 proc capability  noprint data=_pev_;
 var
 %if &ngroup>0 %then %let bis=&ngroup;
 %else %if &ncand>1 %then %let bis=&ncand;
 %else %if &varlist ne %then %let bis=&nvar;
 %do j=1 %to %eval(&bis-1);
  %do k=%eval(1+&j) %to &bis;
   dm&j._&k dp&j._&k
  %end;
 %end;
 ;
 histogram
 %if &ngroup>0 %then %let bis=&ngroup;
 %else %if &ncand>1 %then %let bis=&ncand;
 %else %if &varlist ne %then %let bis=&nvar;
 %do j=1 %to %eval(&bis-1);
  %do k=%eval(1+&j) %to &bis;
   dm&j._&k dp&j._&k
  %end;
 %end;
 ;
 run;

%end;


data _xxx;
set _xxx;
time=time()-time;
time=floor(time*100+0.5)/100;
minutes=floor(time/60);
seconds=time-minutes*60;
call symput('minute',compress(minutes));
call symput('second',compress(seconds));
run;

options nonotes nosource;
proc datasets;
delete _bootsr _bootsr1 _ndat  _pevfull  _srcbbm_ _srcs 
        _take_ _ttt_ _worksr _xxx _glob _srmy
	   %if &cand ne %then %do; 
        %do c=1 %to &ncand; 
         _pevmc&c _pevp&c
        %end; 
	    %if &varlist ne %then %do;
		 %do v=1 %to &nvar;
		  _pevmv&v
		 %end;
		%end;
	   %end;
	   %else %do;
	    %do v=1 %to &nvar;
	     _pevm&v _pevp&v
    	%end;
	   %end;
	   %if &class ne %then %do;
	    %do c=1 %to &nclass;
		 _clt&c
		%end;
		_df _dummies
	   %end;
	   ;
run;
quit;

options notes source;


%put NOTE: Execution of macro used &minute minutes &second seconds.;
%mend;
