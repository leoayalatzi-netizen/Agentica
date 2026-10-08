/* Agentes adicionales + tipo "orbit" (mecánica orbital) + presets */
const vdc=(n,b=2)=>{let q=0,f=1/b;while(n>0){q+=(n%b)*f;n=Math.floor(n/b);f/=b}return q};
function psoMin(f,a,b,abs){const N=14,X=[],V=[],B=[],bf=[],F=x=>abs?Math.abs(f(x)):f(x);let G=0,gf=Infinity;for(let i=0;i<N;i++){X[i]=rnd(a,b);V[i]=0;B[i]=X[i];bf[i]=F(X[i]);if(bf[i]<gf){gf=bf[i];G=X[i]}}
 return function*(){for(let k=0;k<1500;k++){for(let i=0;i<N;i++){V[i]=.7*V[i]+1.5*Math.random()*(B[i]-X[i])+1.5*Math.random()*(G-X[i]);X[i]=Math.min(b,Math.max(a,X[i]+V[i]));const v=F(X[i]);if(v<bf[i]){bf[i]=v;B[i]=X[i]}if(v<gf){gf=v;G=X[i]}}yield G}}()}
Object.assign(SOLVE,{
'Híbrido Newton-Bisección':function*(P){const g=P.g,br=yield*scan(g,P.a,P.b);if(!br)return;let[l,r]=br,fl=g(l),x=(l+r)/2;for(let k=0;k<200;k++){const f=g(x),d=dg(g,x);if(f*fl<=0)r=x;else{l=x;fl=f}let nx=x-f/d;if(!isFinite(nx)||nx<=Math.min(l,r)||nx>=Math.max(l,r))nx=(l+r)/2;x=nx;yield x}},
'Multi-inicio Newton':function*(P){const g=P.g;let x=rnd(P.a,P.b);for(let k=0;k<400;k++){if(k%12==11)x=rnd(P.a,P.b);const d=dg(g,x);let nx=x-g(x)/d;if(!isFinite(nx))nx=rnd(P.a,P.b);x=nx;yield x}},
'Enjambre (PSO)':P=>psoMin(P.g,P.a,P.b,true)});
Object.assign(INT,{
'Tanh-Sinh':function*({f,a,b}){const c=(a+b)/2,hw=(b-a)/2;for(let k=0;k<=7;k++){const h=2**-k;let S=0;for(let t=-3.3;t<=3.3;t+=h){const u=Math.PI/2*Math.sinh(t),x=c+hw*Math.tanh(u),w=hw*Math.PI/2*Math.cosh(t)/Math.cosh(u)**2,y=f(x);if(isFinite(y)&&isFinite(w))S+=w*y}yield S*h}},
'Cuasi-Monte Carlo':function*({f,a,b}){let s=0,n=1;for(let k=0;k<400;k++){for(let j=0;j<4000;j++)s+=f(a+(b-a)*vdc(n++));yield(b-a)*s/(n-1)}}});
DER['7 puntos']=function*({f,x0}){for(let k=0;k<30;k++){const h=.2*.7**k;yield(f(x0+3*h)-9*f(x0+2*h)+45*f(x0+h)-45*f(x0-h)+9*f(x0-2*h)-f(x0-3*h))/(60*h)}};
Object.assign(MIN,{
'Recocido simulado':function*({f,a,b}){let x=rnd(a,b),fx=f(x),bx=x,bf=fx,T=Math.abs(fx)+1,r=(b-a)/3;for(let k=0;k<8000;k++){const y=Math.min(b,Math.max(a,x+(Math.random()*2-1)*r)),fy=f(y);if(fy<fx||Math.random()<Math.exp((fx-fy)/T)){x=y;fx=fy}if(fx<bf){bf=fx;bx=x}T*=.998;r=Math.max(r*.9985,1e-13);yield bx}},
'Enjambre (PSO)':({f,a,b})=>psoMin(f,a,b,false)});
const dp=(f,t0,y0,t1,tol)=>{let t=t0,y=y0,h=(t1-t0)/20,it=0;while(t<t1&&it++<2e6){if(t+h>t1)h=t1-t;const k1=f(t,y),k2=f(t+h/5,y+h*k1/5),k3=f(t+3*h/10,y+h*(3*k1/40+9*k2/40)),k4=f(t+4*h/5,y+h*(44*k1/45-56*k2/15+32*k3/9)),k5=f(t+8*h/9,y+h*(19372*k1/6561-25360*k2/2187+64448*k3/6561-212*k4/729)),k6=f(t+h,y+h*(9017*k1/3168-355*k2/33+46732*k3/5247+49*k4/176-5103*k5/18656)),y5=y+h*(35*k1/384+500*k3/1113+125*k4/192-2187*k5/6784+11*k6/84),k7=f(t+h,y5),y4=y+h*(5179*k1/57600+7571*k3/16695+393*k4/640-92097*k5/339200+187*k6/2100+k7/40),e=Math.abs(y5-y4)/(tol*(1+Math.abs(y)))+1e-30;if(e<=1){t+=h;y=y5}h*=Math.min(5,Math.max(.2,.9*e**-.2))}return y};
ODE['Dormand-Prince RK45']=function*({f,t0,y0,t1}){for(let e=4;e<=14;e++)yield dp(f,t0,y0,t1,10**-e)};
ODE['Adams-Bashforth 4']=function*({f,t0,y0,t1}){for(let k=3;k<=16;k++){const n=2**k,h=(t1-t0)/n;let y=y0;const F=[];for(let i=0;i<n;i++){const t=t0+i*h;F.push(f(t,y));if(i<3){const k1=F[i],k2=f(t+h/2,y+h*k1/2),k3=f(t+h/2,y+h*k2/2),k4=f(t+h,y+h*k3);y+=h*(k1+2*k2+2*k3+k4)/6}else{y+=h*(55*F[3]-59*F[2]+37*F[1]-9*F[0])/24;F.shift()}}yield y}};
const heunSys=(f,t0,y0,v0,t1,n)=>{const h=(t1-t0)/n;let y=y0,v=v0;for(let i=0;i<n;i++){const t=t0+i*h,a1=f(t,y,v),v1=v+h*a1,a2=f(t+h,y+h*v,v1);y+=h*(v+v1)/2;v+=h*(a1+a2)/2}return y};
ODE2['Heun']=function*({f,t0,y0,v0,t1}){for(let k=2;k<=17;k++)yield heunSys(f,t0,y0,v0,t1,2**k)};
Object.assign(LIN,{
'Refinamiento iterativo':function*(P){let x=gsolve(P.A,P.b);if(!x)return;yield x.slice();for(let k=0;k<5;k++){const r=mv(P.A,x).map((v,i)=>P.b[i]-v),d=gsolve(P.A,r);if(!d)return;x=x.map((v,i)=>v+d[i]);yield x.slice()}},
'Kaczmarz':function*(P){const n=P.b.length,x=Array(n).fill(0);for(let k=0;k<4000;k++){for(let j=0;j<50;j++){const i=Math.floor(Math.random()*n),A=P.A[i],nn=dot(A,A);if(!nn)continue;const r=(P.b[i]-dot(A,x))/nn;for(let c=0;c<n;c++)x[c]+=r*A[c]}yield x.slice()}}});
Object.assign(LIM,{
'Aitken lateral':function*({f,x0}){const h=[];for(let k=1;k<=40;k++){h.push(f(x0+.1*.5**k));if(h.length>=3){const v=aitken(h);if(isFinite(v))yield v}}},
'Neville derecha':function*({f,x0}){const T=[];for(let i=0;i<10;i++){const h=.1/2**i;T[i]=[f(x0+h)];for(let j=1;j<=i;j++)T[i][j]=(2**j*T[i][j-1]-T[i-1][j-1])/(2**j-1);yield T[i][i]}}});
const wynn=h=>{let prev=h.map(()=>0),cur=h.slice(),best=h[h.length-1];for(let k=0;k<h.length-1;k++){const nx=[];for(let n=0;n<cur.length-1;n++){const d=cur[n+1]-cur[n];nx.push((prev[n+1]??0)+(d==0?1e300:1/d))}prev=cur;cur=nx;if(!cur.length)break;if(k%2==1)best=cur[cur.length-1]}return best};
SER['Wynn ε']=function*(P){const h=[];for(const S of ps(P,'k')){h.push(S);if(h.length<3)continue;const v=wynn(h.slice(-7));if(isFinite(v))yield v}};
/* ---------- Mecánica orbital de dos cuerpos ---------- */
const accO=(mu,x,y)=>{const r2=x*x+y*y,r3=r2*Math.sqrt(r2);return[-mu*x/r3,-mu*y/r3]};
const lf=(mu,s,h)=>{let[x,y,vx,vy]=s,a=accO(mu,x,y);vx+=h/2*a[0];vy+=h/2*a[1];x+=h*vx;y+=h*vy;a=accO(mu,x,y);return[x,y,vx+h/2*a[0],vy+h/2*a[1]]};
const W1=1/(2-Math.cbrt(2)),W0=-Math.cbrt(2)*W1;
function orbStep(mu,s,h,m){
 if(m=='vl')return lf(mu,s,h);if(m=='yo')return lf(mu,lf(mu,lf(mu,s,W1*h),W0*h),W1*h);
 let[x,y,vx,vy]=s,a;
 if(m=='eu'){a=accO(mu,x,y);return[x+h*vx,y+h*vy,vx+h*a[0],vy+h*a[1]]}
 if(m=='sy'){a=accO(mu,x,y);vx+=h*a[0];vy+=h*a[1];return[x+h*vx,y+h*vy,vx,vy]}
 const F=q=>{const a=accO(mu,q[0],q[1]);return[q[2],q[3],a[0],a[1]]},ad=(q,k,c)=>q.map((v,i)=>v+c*k[i]),k1=F(s),k2=F(ad(s,k1,h/2)),k3=F(ad(s,k2,h/2)),k4=F(ad(s,k3,h));
 return s.map((v,i)=>v+h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]))}
const orb=(P,n,m,s0,T)=>{let s=s0||[P.r0,0,0,P.v0];const h=(T??P.t1)/n;for(let i=0;i<n;i++)s=orbStep(P.mu,s,h,m);return s};
const orbE=(mu,s)=>(s[2]**2+s[3]**2)/2-mu/Math.hypot(s[0],s[1]),orbL=s=>s[0]*s[3]-s[1]*s[2];
const mkO=(m,K)=>function*(P){for(let k=5;k<=K;k++){const s=orb(P,2**k,m);yield Math.hypot(s[0],s[1])}};
const ORB={'Euler':mkO('eu',17),'Euler simpléctico':mkO('sy',17),'Leapfrog (Verlet)':mkO('vl',16),'Yoshida 4.º orden':mkO('yo',14),'RK4':mkO('rk',15),
'RK4+Richardson':function*(P){for(let k=5;k<=13;k++){const a=orb(P,2**k,'rk'),b=orb(P,2**(k+1),'rk');yield(16*Math.hypot(b[0],b[1])-Math.hypot(a[0],a[1]))/15}}};
TYPES.orbit={n:'Órbita de dos cuerpos (propagación)',ag:ORB,kind:'cons',tol:1e-9,lab:'r(t₁) ≈ (km)',cal:{P:{mu:1,r0:1,v0:1,t1:1},t:1},
 tests:(P,v,c)=>{const s0=[P.r0,0,0,P.v0],s1=orb(P,20000,'rk'),E0=orbE(P.mu,s0),E1=orbE(P.mu,s1),L0=orbL(s0),L1=orbL(s1),b=orb(P,20000,'rk',s1,-P.t1);
  return[cl2(c),['energía específica conservada (deriva < 1e-8)',Math.abs((E1-E0)/E0)<1e-8],['momento angular conservado (< 1e-8)',Math.abs((L1-L0)/L0)<1e-8],['reversibilidad temporal (vuelve al punto inicial)',Math.hypot(b[0]-P.r0,b[1])<1e-6*P.r0]]},
 extra:(P,v)=>{const s0=[P.r0,0,0,P.v0],E=orbE(P.mu,s0),h=P.r0*P.v0,s1=orb(P,20000,'rk');let t=`Posición en t₁: (${s1[0].toFixed(3)}, ${s1[1].toFixed(3)}) km · |v| = ${Math.hypot(s1[2],s1[3]).toFixed(5)} km/s`;
  if(E<0){const a=-P.mu/(2*E),e=Math.sqrt(Math.max(0,1+2*E*h*h/(P.mu*P.mu)));t+=`<br>Órbita elíptica: a = ${a.toFixed(1)} km · e = ${e.toFixed(5)} · periodo = ${(2*Math.PI*Math.sqrt(a**3/P.mu)).toFixed(1)} s · perigeo = ${(a*(1-e)).toFixed(1)} km · apogeo = ${(a*(1+e)).toFixed(1)} km`}else t+='<br>Órbita abierta (energía ≥ 0): trayectoria de escape.';return t}};
CFG.orbit=['μ (km³/s²) del cuerpo central',['r₀ (km, posición inicial)','v₀ (km/s, tangencial)','t₁ (s)',''],['398600.4418','','']];
const MU='398600.4418';
PRE.push(
['Sistema 2×2 (rectas que se cortan)','lin','1,2;3,-1','5,1'],
['Límite: (1−cos x)/x² en 0','limit','(1-cos(x))/x^2',0],
['Serie alternada Σ(−1)^(n+1)/n = ln 2','series','(-1)^(n+1)/n',1,'inf'],
['Física: trabajo ∫F dx, resorte duro F=50x+200x³ (0→0.5 m)','integral','50*x+200*x^3',0,.5],
['Física: velocidad instantánea, x(t)=5t³−2t en t=2','deriv','5*x^3-2*x',2],
['Física: carga de un condensador RC (τ=0.5 s, 5 V)','ode','(5-y)/0.5',0,0,2],
['Aero: órbita circular baja (400 km, un periodo)','orbit',MU,6778.137,'sqrt('+MU+'/6778.137)','2*pi*sqrt(6778.137^3/'+MU+')'],
['Aero: órbita elíptica (perigeo 6778 km, 1.2·v_circ)','orbit',MU,6778.137,'1.2*sqrt('+MU+'/6778.137)',10000],
['Aero: ecuación de Kepler E−e·sinE=M (e=0.6, M=1.2)','solve','x-0.6*sin(x)=1.2',0,7],
['Aero: Tsiolkovsky, masa final para Δv=9400 m/s (vₑ=3000, m₀=1e5 kg)','solve','3000*ln(100000/x)=9400',1,99999],
['Aero: ascenso vertical de un cohete (60 s de quemado)','ode2','1200000/(50000-400*t)-9.81*(6371000/(6371000+y))^2',0,0,0,60],
['Aero: reentrada, velocidad con arrastre cuadrático','ode','9.81-0.0025*y^2',0,0,20],
['Aero: impulso total ∫2.5e6·e^(−t/40) dt (0→120 s)','integral','2.5e6*exp(-x/40)',0,120]);
