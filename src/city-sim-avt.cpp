#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <time.h>
#include <string>
using namespace std;


#define pow2(n) (1 << (n))

#define bit_map(x, n) (((x) & (1 << (n))) != 0)

#define bx(x)        ((x) - (((x)>>1) & 0x77777777)  \
			   - (((x)>>2) & 0x33333333)  \
	 		   - (((x)>>3) & 0x11111111))
#define count_bits(x)  (((bx(x)+(bx(x)>>4)) & 0x0F0F0F0F) % 255)

#define hdist(x, y) (count_bits((x)^(y)))


////////////////////// functions

void alloc_mem();

void alloc_memm();

void make_graph();

void make_ma();

void make_mobility();

void make_to();

void make_moves_D(int cycle);

void make_path(int a);

void compute_score();

void report();

double rang(int& seed);


////////////////////// end functions




////////////////////// variables

const int L=10;
const int N=L*L;

int M;
const int  K=4;
const int  Ncycle=300;
const int  T0=200;

const double Dt=1.0;
const double lambda=0.5;

int no,nd;
int seed;

double g=0.0;
double mu=3.0;
double alpha=0.0;

double TO,TOD;
double TD_min,TD_max;

double eta,etaT,etaF,Sv,Sm,Sq,Sd,Pt,Qt,St;
double ND,FD,R;

vector<int> X,Y;

vector<int>  ka;
vector<vector<int> > va,vl;

vector<int>  Oi,Di,ma;
vector<vector<int> > fin,fout,ro;

vector<int> Mi,Li;

vector<double> td_min,td_max;
vector<double> tio,tid,tiw,tif;

vector<vector<double> > fo,fD;
vector<vector<double> > to,te,tD,tab,TD;


////////////////////// end variables






////////////////////// main

int main()
{

int nav;  
double e1,e2,et1,et2,ef1,ef2,d1,d2,m1,m2,s1,s2,v1,v2,p1,p2,q1,q2,n1,n2;
double eta1,eta2,eT1,eT2,eF1,eF2,Sv1,Sv2,Sq1,Sq2,Sm1,Sm2,Sd1,Sd2,Pt1,Pt2,Qt1,Qt2,St1,St2;


alloc_mem();

seed=time(NULL);
  
make_graph();

ofstream output("avC.dat");

g=2.0;
alpha=0.0;
mu=3.0;
TO=8;
no=TO/Dt;

for(int la=0;la<21;la++){
alpha=la*0.02;
//for(int lg=1;lg<21;lg++){
//g=lg*0.2;


eta1=0;
eta2=0;
eT1=0;
eT2=0;
eF1=0;
eF2=0;
Sv1=0;
Sv2=0;
Sq1=0;
Sq2=0;
Sm1=0;
Sm2=0;
Sd1=0;
Sd2=0;
Pt1=0;
Pt2=0;
Qt1=0;
Qt2=0;
St1=0;
St2=0;

int Nsample=100;
for(int sample=0;sample<Nsample;sample++){

make_ma();

alloc_memm();

make_mobility();
  
make_to();

for(int a=0;a<N;a++){
  for(int la=0;la<ka[a];la++){
    tD[a][la]=to[a][la];
  }
}

eta=0;
etaT=0;
etaF=0;
Sv=0;
Sd=0;
Sm=0;
Sq=0;
Pt=0;
Qt=0;
St=0;
nav=0;
for(int cycle=0;cycle<Ncycle;cycle++){		/// loop over cycles
  
  
  ////////////////////////////// first stage
  
  /// initials
  
  FD=0;
  ND=0;
  TOD=0;
  TD_max=0;
  TD_min=1e+10;
  for(int a=0;a<N;a++){
      td_max[a]=0;
      td_min[a]=1e+10;
      for(int la=0;la<ka[a];la++){
        fD[a][la]=0;
        TD[a][la]=0;
        tab[a][la]=tD[a][la];
      }
      for(int b=0;b<N;b++)te[a][b]=0;
  }
  
  /// movments to destinations

  
  for(int a=0;a<N;a++)make_path(a);

  
  make_moves_D(cycle);		  

  
  /// time averaging
  
  for(int a=0; a<N; a++){
    for(int la=0; la<ka[a]; la++){
      
      if(fD[a][la]>0){
	tD[a][la]=(1-lambda)*tab[a][la]+lambda*TD[a][la]/fD[a][la];
      }else{
	tD[a][la]=to[a][la];
      }
      
    }
  }
  

  ////////////////////////////// end of cycle

  if(cycle>T0){
    nav+=1;
    compute_score();
  }
}
eta=eta/nav;
etaT=etaT/nav;
etaF=etaF/nav;
Sv=Sv/nav;
Sd=Sd/nav;
Sm=Sm/nav;
Sq=Sq/nav;
Pt=Pt/nav;
Qt=Qt/nav;
St=St/nav;


eta1+=eta;
eta2+=eta*eta;
eT1+=etaT;
eT2+=etaT*etaT;
eF1+=etaF;
eF2+=etaF*etaF;
Sv1+=Sv;
Sv2+=Sv*Sv;
Sd1+=Sd;
Sd2+=Sd*Sd;
Sm1+=Sm;
Sm2+=Sm*Sm;
Sq1+=Sq;
Sq2+=Sq*Sq;
Pt1+=Pt;
Pt2+=Pt*Pt;
Qt1+=Qt;
Qt2+=Qt*Qt;
St1+=St;
St2+=St*St;

}

e1=eta1/Nsample;
e2=eta2/Nsample;
e2=sqrt((e2-e1*e1)/(Nsample-1));

et1=eT1/Nsample;
et2=eT2/Nsample;
et2=sqrt((et2-et1*et1)/(Nsample-1));

ef1=eF1/Nsample;
ef2=eF2/Nsample;
ef2=sqrt((ef2-ef1*ef1)/(Nsample-1));

d1=Sd1/Nsample;
d2=Sd2/Nsample;
d2=sqrt((d2-d1*d1)/(Nsample-1));

m1=Sm1/Nsample;
m2=Sm2/Nsample;
m2=sqrt((m2-m1*m1)/(Nsample-1));

s1=Sq1/Nsample;
s2=Sq2/Nsample;
s2=sqrt((s2-s1*s1)/(Nsample-1));

v1=Sv1/Nsample;
v2=Sv2/Nsample;
v2=sqrt((v2-v1*v1)/(Nsample-1));

p1=Pt1/Nsample;
p2=Pt2/Nsample;
p2=sqrt((p2-p1*p1)/(Nsample-1));

q1=Qt1/Nsample;
q2=Qt2/Nsample;
q2=sqrt((q2-q1*q1)/(Nsample-1));

n1=St1/Nsample;
n2=St2/Nsample;
n2=sqrt((n2-n1*n1)/(Nsample-1));


cout<<alpha<<"  "<<e1<<"  "<<e2<<"  "<<et1<<"  "<<et2<<" "<<ef1<<"  "<<ef2<<"  "<<v1<<"  "<<v2<<"  "<<s1<<"  "<<s2<<" "<<m1<<"  "<<m2<<" "<<d1<<"  "<<d2<<"  "<<p1<<"  "<<p2<<"  "<<q1<<"  "<<q2<<"  "<<n1<<"  "<<n2<<endl;
output<<alpha<<"  "<<e1<<"  "<<e2<<"  "<<et1<<"  "<<et2<<" "<<ef1<<"  "<<ef2<<"  "<<v1<<"  "<<v2<<"  "<<s1<<"  "<<s2<<" "<<m1<<"  "<<m2<<" "<<d1<<"  "<<d2<<" "<<p1<<"  "<<p2<<"  "<<q1<<"  "<<q2<<"  "<<n1<<"  "<<n2<<endl;

}



 return 0;
}
////////////////////// end  main



/////////////////// alloc_mem

void alloc_mem()
{
  
  
ma.resize(N);

X.resize(N);
Y.resize(N);

ka.resize(N);
va.resize(N);
vl.resize(N);
for(int a=0;a<N;a++){
va[a].resize(K);
vl[a].resize(K);
}

fo.resize(N);
fD.resize(N);
to.resize(N);
te.resize(N);
tD.resize(N);
TD.resize(N);
tab.resize(N);
for(int a=0;a<N;a++){
fo[a].resize(K);
fD[a].resize(K);
to[a].resize(K);
te[a].resize(N);
tD[a].resize(K);
TD[a].resize(K);
tab[a].resize(K);
}

td_min.resize(N);
td_max.resize(N);

ro.resize(N);
fin.resize(N);
fout.resize(N);
for(int a=0; a<N; a++){
ro[a].resize(K);
fin[a].resize(K);
fout[a].resize(K);
}


}
/////////////////// end  alloc_mem




/////////////////// alloc_memm

void alloc_memm()
{
 
  
Mi.resize(M);
Li.resize(M);
  
Oi.resize(M);
Di.resize(M);

tio.resize(M);
tid.resize(M);
tiw.resize(M);
tif.resize(M);


}
/////////////////// end  alloc_memm



////////////////////// make_graph

void make_graph()
{

int a,b,la,lb;
int xx,yy;

for(int a=0;a<N;a++)ka[a]=0;
  

for(int x=0;x<L;x++){
  for(int y=0;y<L;y++){
    
    a=y*L+x;
    X[a]=x;
    Y[a]=y;
    
    ///
    
    xx=x+1;
    yy=y;
    if(xx<L){
    b=yy*L+xx;
          
    la=ka[a];
    lb=ka[b];
    va[a][la]=b;
    va[b][lb]=a;
    vl[a][la]=lb;
    vl[b][lb]=la;
    ka[a]=la+1;
    ka[b]=lb+1;
          
    to[a][la]=1.0;
    to[b][lb]=1.0;
    fo[a][la]=1.0;
    fo[b][lb]=1.0;
    }
    
    ///

    xx=x;
    yy=y+1;
    if(yy<L){
    b=yy*L+xx;
          
    la=ka[a];
    lb=ka[b];
    va[a][la]=b;
    va[b][lb]=a;
    vl[a][la]=lb;
    vl[b][lb]=la;
    ka[a]=la+1;
    ka[b]=lb+1;
          
    to[a][la]=1.0;
    to[b][lb]=1.0;
    fo[a][la]=1.0;
    fo[b][lb]=1.0;
    }
    
    ///
    
    
  }
}


}
////////////////////// end make_graph




////////////////////// make_ma

void make_ma()
{

int a,b;
int x,y,z;
int C=1;
int l0=1;
int m0=1000;
short int check;
double p,summ;


M=0;
summ=N*C;
for(int a=0;a<N;a++)ma[a]=0;

a=(L/2)*L+L/2;
for(int lx=-l0;lx<=l0;lx++){
for(int ly=-l0;ly<=l0;ly++){
    x=X[a]+lx;
    y=Y[a]+ly;
    b=y*L+x;
    ma[b]=1;
    summ+=1;
    M+=1;
}
}


while(M/double(N)<m0){

for(int s=0;s<N;s++){
  a=rang(seed)*N;
  
  check=0;
  for(int lx=-l0;lx<=l0;lx++){
  for(int ly=-l0;ly<=l0;ly++){
    x=X[a]+lx;
    y=Y[a]+ly;
    if((x>=0)&&(x<L)){
    if((y>=0)&&(y<L)){
      b=y*L+x;
      if(ma[b]>0)check=1;
    }
    }
  }
  }
    
  p=(ma[a]+C)/summ;
  if((check==1)&&(rang(seed)<p)){  
    ma[a]+=1;
    summ+=1;
    M+=1;
  }
}

 
}



}
////////////////////// end make_ma




////////////////////// make_mobility

void make_mobility()
{

int i,c,lb,lc,lm,x,y,z;
short int check;
double Za,Zab;
double Dab,Dcb;
double qab;

vector<vector<int> > Ia;
Ia.resize(N);
for(int a=0;a<N;a++){
Ia[a].resize(ma[a]);
}


vector<double> Qab;
Qab.resize(N);


i=0;
for(int a=0;a<N;a++){
  for(int la=0;la<ma[a];la++){
    Oi[i]=a;
    Di[i]=a;
    Ia[a][la]=i;
    i+=1;
  }
  for(int la=0;la<ka[a];la++)fo[a][la]=double(M)/(N*K);
}

///

for(int a=0;a<N;a++){ 

  Za=0;
  Qab[a]=0;
  for(int b=0;b<N;b++){ 
    if(b!=a){
      
      x=(X[a]-X[b]);
      y=(Y[a]-Y[b]);
      Dab=sqrt(x*x+y*y);
      
      Zab=1;
      for(int c=0;c<N;c++){
        x=(X[c]-X[b]);
        y=(Y[c]-Y[b]);
        Dcb=sqrt(x*x+y*y);
	if((c!=a)&&(c!=b)&&(Dcb<=Dab))Zab+=ma[c];
      }
      
      qab=ma[b]/Zab;
      Qab[b]=ma[a]*qab;
      Za+=qab;
    }
  }
  
  ///
  
  lm=0;
  for(int b=0;b<N;b++){
    z=Qab[b]/Za;
    for(int la=lm;la<lm+z;la++){
      i=Ia[a][la];
      Di[i]=b;
    }
    lm+=z;
  }
  
  for(int la=lm;la<ma[a];la++){
    i=Ia[a][la];
    check=0;
    while(check==0){
      c=rang(seed)*N;
      if(c!=a){
	Di[i]=c;
	check=1;
      }
    }
  }

  ///
  
}


///



}
////////////////////// end make_mobility




////////////////////// make_to

void make_to()
{

double r;

for(int i=0;i<M;i++){
  r=0;
  for(int l=0;l<12;l++)r+=rang(seed);
  tio[i]=TO*r/12.0;
  tio[i]=TO*rang(seed);
}



/////

/*
int ai;  
double mmax;
  
mmax=0;
for(int i=0;i<M;i++){
  ai=Oi[i];
  if(ma[ai]>mmax)mmax=ma[ai];
}
for(int i=0;i<M;i++){
  ai=Oi[i];
  tio[i]=TO*rang(seed)*ma[ai]/mmax;
}
*/



}
////////////////////// end make_to




/////////////////// make_moves_D
void make_moves_D(int cycle)
{

  
int nar,nac;
int ai,bi,ci,di,la,lb,lc,ld,lmin;
double t,ti,tmin;

vector<vector<double> > tab_new;

tab_new.resize(N);
for(int a=0; a<N; a++){
tab_new[a].resize(K);
}

vector<int> pos;
pos.resize(M);

vector<short int> active,lin,lout;
active.resize(M);
lin.resize(M);
lout.resize(M);

///


nd=0;
nar=0;
nac=0;
for(int i=0; i<M; i++){
  Mi[i]=0;
  Li[i]=0;
  active[i]=0;
  pos[i]=Oi[i];
}

R=0;
for(int a=0;a<N;a++){
  for(int la=0;la<ka[a];la++)ro[a][la]=0;
}

    
t=0;
while(nar<M){		/// time evolution of the first-stage movments
     
    nd+=1;
  
    ///
  
    for(int i=0;i<M;i++){
      if(active[i]==0){
        ti=tio[i];
        if((t-Dt<ti)&&(ti<t+Dt)){
          active[i]=1;
	  nac+=1;
	  tid[i]=ti;
	  tiw[i]=+1;
	  tif[i]=-1;
        }
      }
    }

    for(int a=0;a<N;a++){
      for(int la=0;la<ka[a];la++){
	fin[a][la]=0;
	fout[a][la]=0;
      }
    }  
    
    ///
    
    for(int i=0; i<M; i++){
      if(active[i]==1){	
	
      if(tiw[i]>tif[i]){
	
        ai=pos[i];
	di=Di[i];
	if(ai==di){
	  
	  active[i]=2;
	  nar+=1;
	  nac-=1;
	  ti=tid[i];
          if(ti<TD_min)TD_min=ti;
          if(ti>TD_max)TD_max=ti;
          if(ti<td_min[di])td_min[di]=ti;
          if(ti>td_max[di])td_max[di]=ti;

          la=lin[i];
          ci=va[ai][la];
          lc=vl[ai][la];
	  fout[ci][lc]+=1;
	  
	}else{

	  if(tif[i]>0){  
	    la=lin[i];
            ci=va[ai][la];
            lc=vl[ai][la];
	    fout[ci][lc]+=1;
	  }	  
	  tmin=1e+10;
          for(int li=0; li<ka[ai]; li++){
	    bi=va[ai][li];
	    ti=tD[ai][li]+te[bi][di];
	    if(ti<tmin){
	      lmin=li;
	      tmin=ti;
	    }
	  }
	  if(rang(seed)<alpha)lmin=rang(seed)*ka[ai];
	  bi=va[ai][lmin];
	  lb=vl[ai][lmin];
	  lin[i]=lb;
	  lout[i]=lmin;
	  fin[ai][lmin]+=1;
	    
	}
	
	tiw[i]=Dt/2.0;
      }else{
	tiw[i]+=Dt;
      }
      
	tid[i]+=Dt;
      }
    }
    
    ///
        
    for(int a=0; a<N; a++){
      for(int la=0; la<ka[a]; la++){
	ro[a][la]+=fin[a][la]-fout[a][la];
	FD+=fin[a][la];
	fD[a][la]+=fin[a][la];
	if(fin[a][la]>0)ND+=1;
        tab_new[a][la]=to[a][la]*(1+g*pow(fin[a][la]/fo[a][la],mu));
        TD[a][la]+=fin[a][la]*tab_new[a][la];
        TOD+=fin[a][la]*tab_new[a][la];
      }
    }
    
    ///

    for(int i=0; i<M; i++){
      if((active[i]==1)&&(tiw[i]<Dt)){
	
        ai=pos[i];
	la=lout[i];
	bi=va[ai][la];
	tif[i]=tab_new[ai][la];
        Mi[i]+=fin[ai][la];	
        Li[i]+=to[ai][la];	
	pos[i]=bi;

      }
    }
    

    ///
    
	    
    t+=Dt;
    
}    


  
}
/////////////////// make_moves_D




/////////////////// make_path
void make_path(int a)
{

int c0,cl,l0,n;  
double dmin;

vector<double> dist;
vector<int> list,point;
vector<short int> visit;

dist.resize(N);
list.resize(N);
point.resize(N);
visit.resize(N);



for(int b=0; b<N; b++){
  list[b]=b;
  visit[b]=0;
  dist[b]=1e+10;
  point[b]=-1;
}

n=N;
dist[a]=0;
point[a]=a;
while(n>0){  
  
  l0=0;
  c0=list[l0];
  dmin=dist[c0];
  for(int l=1;l<n;l++){
    cl=list[l];
    if(dist[cl]<dmin){
      l0=l;
      c0=cl;
      dmin=dist[cl];
    }
  }
  list[l0]=list[n-1];
  visit[c0]=1;
  n=n-1;
  
  for(int l=0; l<ka[c0]; l++){
    cl=va[c0][l];
    if(visit[cl]==0){
      
       if(dist[c0]+tab[c0][l]<dist[cl]){
        point[cl]=c0;
	dist[cl]=dist[c0]+tab[c0][l];
      }
      
    }
  }
 

}
for(int b=0; b<N; b++)te[a][b]=dist[b];


}
/////////////////// make_path





////////////////////// compute_score

void compute_score()
{

int noo,nxx,ndd,ntt;
int ai,bi,n1,n2,na;
double ti,xi,li,v,m,q;
double sum,sum1,sum2,sum3,sum12;

vector<int> P1,P2,Q1,Q2;
vector<vector<int> > P12,Q12;

noo=0;
nxx=0;
ndd=0;
ntt=0;
for(int i=0;i<M;i++){
  
  ti=tid[i]-tio[i];
  ai=Oi[i];
  bi=Di[i];
  xi=sqrt(pow(X[ai]-X[bi],2.0)+pow(Y[ai]-Y[bi],2.0));

  n1=int(tio[i]/Dt);
  n2=int((tid[i]-TD_min)/Dt);

  noo=max(noo,n1);
  ndd=max(ndd,n2);

  n1=int(xi);
  n2=int(ti/Dt);
  
  nxx=max(nxx,n1);
  ntt=max(ntt,n2);
  
}
noo+=1;
ndd+=1;
nxx+=1;
ntt+=1;
P1.resize(noo);
P2.resize(ndd);
Q1.resize(nxx);
Q2.resize(ntt);
P12.resize(noo);
for(int lo=0;lo<noo;lo++)P12[lo].resize(ndd);
Q12.resize(nxx);
for(int lx=0;lx<nxx;lx++)Q12[lx].resize(ntt);


///

for(int t1=0;t1<noo;t1++){
  for(int t2=0;t2<ndd;t2++){
    P1[t1]=0;
    P2[t2]=0;
    P12[t1][t2]=0;
  }
}

for(int lx=0;lx<nxx;lx++){
  for(int lt=0;lt<ntt;lt++){
    Q1[lx]=0;
    Q2[lt]=0;
    Q12[lx][lt]=0;
  }
}

///

sum1=0;
sum2=0;
sum3=0;
for(int i=0;i<M;i++){
  
  ti=tid[i]-tio[i];
  ai=Oi[i];
  bi=Di[i];
  xi=sqrt(pow(X[ai]-X[bi],2.0)+pow(Y[ai]-Y[bi],2.0));
  v=xi/ti;
  sum1+=v;
  m=Mi[i]/double(M);
  sum2+=m;
  q=xi/Li[i];
  sum3+=q;
    
  ///
  
  n1=int(tio[i]/Dt);
  n2=int((tid[i]-TD_min)/Dt);
  P1[n1]+=1;
  P2[n2]+=1;
  P12[n1][n2]+=1;
    
  n1=int(xi);
  n2=int(ti/Dt);
  Q1[n1]+=1;
  Q2[n2]+=1;
  Q12[n1][n2]+=1;
  
}
Sv+=sum1/M;
Sm+=sum2/M;
Sq+=sum3/M;


sum=0;
for(int a=0;a<N;a++){
  na=1+(td_max[a]-td_min[a])/Dt;
  if(na>0)sum+=log(na);
}
Sd+=sum/N-log(no);
St+=log(TD_max-TD_min)-log(TO);

///

TOD=TOD/M;
FD=FD/M;
etaT+=1.0/TOD;
etaF+=1.0/FD;
eta+=1.0/(TOD*FD);

///

sum1=0;
sum2=0;
sum12=0;
for(int t1=0;t1<noo;t1++){
  for(int t2=0;t2<ndd;t2++){
    if(P1[t1]>0)sum1+=P12[t1][t2]*log(P1[t1]);
    if(P2[t2]>0)sum2+=P12[t1][t2]*log(P2[t2]);
    if(P12[t1][t2]>0)sum12+=P12[t1][t2]*log(P12[t1][t2]);
  }
}
Pt+=(sum12-sum1-sum2)/M+log(M);

///

sum1=0;
sum2=0;
sum12=0;
for(int lx=0;lx<nxx;lx++){
  for(int lt=0;lt<ntt;lt++){
    if(Q1[lx]>0)sum1+=Q12[lx][lt]*log(Q1[lx]);
    if(Q2[lt]>0)sum2+=Q12[lx][lt]*log(Q2[lt]);
    if(Q12[lx][lt]>0)sum12+=Q12[lx][lt]*log(Q12[lx][lt]);
  }
}
Qt+=(sum12-sum1-sum2)/M+log(M);



}
////////////////////// end compute_score




////////////////////// report

void report()
{



}
////////////////////// end report





////////////////////// rang

double rang(int& seed)
{

      int a, m, q, r, l;
      double conv, rand;


      a = 16807;
      m = 2147483647;
      q = 127773;
      r = 2836;

      conv = 1.0 / (m - 1);

      l = seed / q;
      seed = a * (seed - q * l) - r * l;
      if (seed < 0) {
	      seed += m;
      }
      rand = conv * (seed - 1);



      return rand;
} 
////////////////////// end rang


