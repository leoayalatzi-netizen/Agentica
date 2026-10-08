/* Arena v3 · motor de agentes autónomos (compilador, agentes base, evolución, run) */

const $=id=>document.getElementById(id);
const FN={sin:Math.sin,cos:Math.cos,tan:Math.tan,exp:Math.exp,ln:Math.log,log:Math.log10,sqrt:Math.sqrt,abs:Math.abs,asin:Math.asin,acos:Math.acos,atan:Math.atan,sinh:Math.sinh,cosh:Math.cosh,tanh:Math.tanh};
function compile(src){
 const s=src.toLowerCase().replace(/\s+/g,'').replace(/\*\*/g,'^').replace(/π/g,'pi');let i=0;
 const ex=()=>{let n=tm();while(s[i]=='+'||s[i]=='-'){const o=s[i++],l=n,r=tm();n=o=='+'?v=>l(v)+r(v):v=>l(v)-r(v)}return n};
 const tm=()=>{let n=un();while(s[i]=='*'||s[i]=='/'){const o=s[i++],l=n,r=un();n=o=='*'?v=>l(v)*r(v):v=>l(v)/r(v)}return n};
 const un=()=>{if(s[i]=='-'){i++;const u=un();return v=>-u(v)}if(s[i]=='+'){i++;return un()}return pw()};
 const pw=()=>{const b=at();if(s[i]=='^'){i++;const e=un();return v=>Math.pow(b(v),e(v))}return b};
 const at=()=>{
  if(s[i]=='('){i++;const n=ex();if(s[i++]!=')')throw Error('falta )');return n}
  let m=/^(\d+\.?\d*|\.\d+)(e[+-]?\d+)?/.exec(s.slice(i));
  if(m){i+=m[0].length;const c=parseFloat(m[0]);return()=>c}
  m=/^[a-z]+/.exec(s.slice(i));if(!m)throw Error('símbolo inesperado: '+s.slice(i,i+5));
  const id=m[0];i+=id.length;
  if(s[i]=='('){if(!FN[id])throw Error('función desconocida: '+id);i++;const a=ex();if(s[i++]!=')')throw Error('falta )');const f=FN[id];return v=>f(a(v))}
  if(id=='pi')return()=>Math.PI;if(id=='e')return()=>Math.E;
  return v=>v[id]};
 const n=ex();if(i<s.length)throw Error('sobra: '+s.slice(i));return n}
const rnd=(a,b)=>a+Math.random()*(b-a),nv=v=>Array.isArray(v)?v.every(isFinite):isFinite(v);
const aitken=h=>{const[x0,x1,x2]=h.slice(-3),d=x2-2*x1+x0;return Math.abs(d)>1e-300?x2-(x2-x1)**2/d:NaN};
/* ===== solve ===== */
function*scan(g,a,b){let px=a,pf=g(a);for(let k=1;k<=2000;k++){const x=a+(b-a)*k/2000,fx=g(x);if(isFinite(pf)&&isFinite(fx)&&pf*fx<=0)return[px,x];px=x;pf=fx;if(k%100==0)yield null}return null}
const dg=(g,x)=>{const h=1e-6*(1+Math.abs(x));return(g(x+h)-g(x-h))/(2*h)};
const SOLVE={
'Bisección':function*(P){const g=P.g,br=yield*scan(g,P.a,P.b);if(!br)return;let[l,r]=br,fl=g(l);for(let k=0;k<200;k++){const m=(l+r)/2,fm=g(m);if(fl*fm<=0)r=m;else{l=m;fl=fm}yield m}},
'Newton':function*(P){const g=P.g;let x=(P.a+P.b)/2;for(let k=0;k<400;k++){const d=dg(g,x);let nx=x-g(x)/d;if(!isFinite(nx)||d==0)nx=rnd(P.a,P.b);x=nx;yield x}},
'Halley':function*(P){const g=P.g;let x=rnd(P.a,P.b);for(let k=0;k<300;k++){const h=1e-4*(1+Math.abs(x)),f=g(x),d=(g(x+h)-g(x-h))/(2*h),dd=(g(x+h)-2*f+g(x-h))/(h*h);let nx=x-2*f*d/(2*d*d-f*dd);if(!isFinite(nx))nx=rnd(P.a,P.b);x=nx;yield x}},
'Secante':function*(P){const g=P.g;let x0=rnd(P.a,P.b),x1=x0+.5;for(let k=0;k<400;k++){const f0=g(x0),f1=g(x1);let x2=x1-f1*(x1-x0)/(f1-f0);if(!isFinite(x2))x2=rnd(P.a,P.b);x0=x1;x1=x2;yield x2}},
'Steffensen':function*(P){const g=P.g;let x=rnd(P.a,P.b);for(let k=0;k<300;k++){const f=g(x),d=(g(x+f)-f)/f;let nx=x-f*f/(g(x+f)-f);if(!isFinite(nx)||f==0)nx=f==0?x:rnd(P.a,P.b);x=nx;yield x}},
'Regula Falsi':function*(P){const g=P.g,br=yield*scan(g,P.a,P.b);if(!br)return;let[a,b]=br,fa=g(a),fb=g(b),sd=0;for(let k=0;k<300;k++){const c=(a*fb-b*fa)/(fb-fa);if(!isFinite(c))return;const fc=g(c);if(fc*fb>0){b=c;fb=fc;if(sd==-1)fa/=2;sd=-1}else{a=c;fa=fc;if(sd==1)fb/=2;sd=1}yield c}},
'Ridders':function*(P){const g=P.g,br=yield*scan(g,P.a,P.b);if(!br)return;let[a,b]=br,fa=g(a),fb=g(b);for(let k=0;k<100;k++){const c=(a+b)/2,fc=g(c),s=Math.sqrt(fc*fc-fa*fb);if(!s||!isFinite(s))return;const x=c+(c-a)*Math.sign(fa-fb)*fc/s,fx=g(x);yield x;if(fc*fx<0){a=c;fa=fc;b=x;fb=fx}else if(fa*fx<0){b=x;fb=fx}else{a=x;fa=fx}}},
'Evolutiva (1+λ)':function*(P){const g=P.g;let x=rnd(P.a,P.b),best=Math.abs(g(x)),r=(P.b-P.a)/4;for(let k=0;k<4000;k++){let ok=0;for(let j=0;j<4;j++){const y=x+(Math.random()*2-1)*r,fy=Math.abs(g(y));if(fy<best){best=fy;x=y;ok=1}}r=ok?r*1.1:Math.max(r*.9,1e-14);yield x}}};
/* ===== lineal ===== */
const mv=(M,v)=>M.map(r=>r.reduce((s,a,j)=>s+a*v[j],0)),dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
function*gj(P,jd){const n=P.b.length,M=P.A.map((r,i)=>[...r,P.b[i]]);for(let c=0;c<n;c++){let p=c;for(let r=c+1;r<n;r++)if(Math.abs(M[r][c])>Math.abs(M[p][c]))p=r;if(Math.abs(M[p][c])<1e-14)return;[M[c],M[p]]=[M[p],M[c]];for(let r=jd?0:c+1;r<n;r++)if(r!=c){const f=M[r][c]/M[c][c];for(let k=c;k<=n;k++)M[r][k]-=f*M[c][k]}yield null}const x=[];for(let i=n-1;i>=0;i--){let s=M[i][n];for(let k=i+1;k<n;k++)s-=M[i][k]*x[k];x[i]=s/M[i][i]}yield x}
const gsolve=(A,b)=>{let x=null;for(const v of gj({A,b}))if(v)x=v;return x};
const itr=(w,jac)=>function*(P){const n=P.b.length;let x=Array(n).fill(0);for(let k=0;k<3000;k++){const y=jac?x.slice():x;for(let i=0;i<n;i++){let s=P.b[i];for(let j=0;j<n;j++)if(j!=i)s-=P.A[i][j]*y[j];y[i]=jac?s/P.A[i][i]:(1-w)*y[i]+w*s/P.A[i][i]}x=y;if(!x.every(isFinite))return;yield x.slice()}};
const LIN={'Gauss':P=>gj(P,0),'Gauss-Jordan':P=>gj(P,1),'Jacobi':itr(1,1),'Gauss-Seidel':itr(1,0),'SOR (ω=1.25)':itr(1.25,0),
'Grad. conjugado':function*(P){const n=P.b.length,A=P.A,AT=A[0].map((_,j)=>A.map(r=>r[j])),x=Array(n).fill(0);let r=mv(AT,P.b),p=r.slice();for(let k=0;k<500;k++){const Ap=mv(AT,mv(A,p)),pa=dot(p,Ap);if(!pa)return;const al=dot(r,r)/pa;for(let i=0;i<n;i++)x[i]+=al*p[i];const r2=r.map((v,i)=>v-al*Ap[i]),be=dot(r2,r2)/dot(r,r);p=r2.map((v,i)=>v+be*p[i]);r=r2;yield x.slice()}}};
/* ===== integral ===== */
const tn=(f,a,b,T,n)=>{const h=(b-a)/n;let m=0;for(let j=0;j<n;j++)m+=f(a+(j+.5)*h);return T/2+h/2*m};
const GX=[0,.5384693101056831,.906179845938664],GW=[.5688888888888889,.4786286704993665,.2369268850561891];
const glq=(f,a,b,m)=>{const h=(b-a)/m;let S=0;for(let j=0;j<m;j++){const c=a+(j+.5)*h;S+=GW[0]*f(c);for(let q=1;q<3;q++)S+=GW[q]*(f(c+GX[q]*h/2)+f(c-GX[q]*h/2))}return S*h/2};
const INT={
'Trapecio':function*({f,a,b}){let T=(b-a)*(f(a)+f(b))/2,n=1;for(let k=0;k<20;k++){T=tn(f,a,b,T,n);n*=2;yield T}},
'Simpson':function*({f,a,b}){let T=(b-a)*(f(a)+f(b))/2,n=1;for(let k=0;k<20;k++){const T2=tn(f,a,b,T,n);n*=2;yield(4*T2-T)/3;T=T2}},
'Romberg':function*({f,a,b}){let R=[(b-a)*(f(a)+f(b))/2],n=1;for(let k=1;k<14;k++){const row=[tn(f,a,b,R[0],n)];n*=2;for(let j=1;j<=k;j++)row[j]=(4**j*row[j-1]-R[j-1])/(4**j-1);R=row;yield R[k]}},
'Gauss-Legendre':function*({f,a,b}){for(let m=1;m<=1<<15;m*=2)yield glq(f,a,b,m)},
'Simpson adaptativo':function*({f,a,b}){for(let e=3;e<=13;e++){const rec=(a,b,fa,fm,fb,w,t,d)=>{const m=(a+b)/2,l=(a+m)/2,r=(m+b)/2,fl=f(l),fr=f(r),sl=(m-a)/6*(fa+4*fl+fm),sr=(b-m)/6*(fm+4*fr+fb);if(d>40||Math.abs(sl+sr-w)<=15*t)return sl+sr+(sl+sr-w)/15;return rec(a,m,fa,fl,fm,sl,t/2,d+1)+rec(m,b,fm,fr,fb,sr,t/2,d+1)};const fa=f(a),fb=f(b),fm=f((a+b)/2);yield rec(a,b,fa,fm,fb,(b-a)/6*(fa+4*fm+fb),10**-e,0)}},
'Monte Carlo':function*({f,a,b}){let s=0,n=0;for(let k=0;k<500;k++){for(let j=0;j<4000;j++){s+=f(a+Math.random()*(b-a));n++}yield(b-a)*s/n}}};
/* ===== derivada ===== */
const DER={
'Adelante':function*({f,x0}){for(let k=1;k<=14;k++){const h=10**-k;yield(f(x0+h)-f(x0))/h}},
'Central':function*({f,x0}){for(let k=1;k<=12;k++){const h=10**-k;yield(f(x0+h)-f(x0-h))/(2*h)}},
'5 puntos':function*({f,x0}){for(let k=0;k<30;k++){const h=.1*.7**k;yield(-f(x0+2*h)+8*f(x0+h)-8*f(x0-h)+f(x0-2*h))/(12*h)}},
'Richardson':function*({f,x0}){const D=[];for(let i=0;i<9;i++){const h=.1/2**i;D[i]=[(f(x0+h)-f(x0-h))/(2*h)];for(let j=1;j<=i;j++)D[i][j]=(4**j*D[i][j-1]-D[i-1][j-1])/(4**j-1);yield D[i][i]}}};
const D5=(f,x,h)=>(-f(x+2*h)+8*f(x+h)-8*f(x-h)+f(x-2*h))/(12*h);
/* ===== EDO 1er orden ===== */
const eu=(f,t0,y0,t1,n)=>{const h=(t1-t0)/n;let y=y0;for(let i=0;i<n;i++)y+=h*f(t0+i*h,y);return y};
const mp=(f,t0,y0,t1,n)=>{const h=(t1-t0)/n;let y=y0;for(let i=0;i<n;i++){const t=t0+i*h;y+=h*f(t+h/2,y+h/2*f(t,y))}return y};
const he=(f,t0,y0,t1,n)=>{const h=(t1-t0)/n;let y=y0;for(let i=0;i<n;i++){const t=t0+i*h,k1=f(t,y),k2=f(t+h,y+h*k1);y+=h*(k1+k2)/2}return y};
const rk=(f,t0,y0,t1,n)=>{const h=(t1-t0)/n;let y=y0;for(let i=0;i<n;i++){const t=t0+i*h,k1=f(t,y),k2=f(t+h/2,y+h*k1/2),k3=f(t+h/2,y+h*k2/2),k4=f(t+h,y+h*k3);y+=h*(k1+2*k2+2*k3+k4)/6}return y};
const mk=(m,K)=>function*({f,t0,y0,t1}){for(let k=2;k<=K;k++)yield m(f,t0,y0,t1,2**k)};
const ODE={'Euler':mk(eu,18),'Punto medio':mk(mp,16),'Heun (RK2)':mk(he,17),'RK4':mk(rk,15),
'Euler+Richardson':function*({f,t0,y0,t1}){for(let k=2;k<=17;k++)yield 2*eu(f,t0,y0,t1,2**(k+1))-eu(f,t0,y0,t1,2**k)},
'RK4+Richardson':function*({f,t0,y0,t1}){for(let k=1;k<=13;k++)yield(16*rk(f,t0,y0,t1,2**(k+1))-rk(f,t0,y0,t1,2**k))/15}};
/* ===== EDO 2º orden (y''=f(t,y,v)) ===== */
const sysE=(f,t0,y0,v0,t1,n,m)=>{const h=(t1-t0)/n;let y=y0,v=v0;for(let i=0;i<n;i++){const t=t0+i*h;
 if(m=='eu'){const a=f(t,y,v);y+=h*v;v+=h*a}
 else if(m=='sy'){v+=h*f(t,y,v);y+=h*v}
 else if(m=='vl'){const a=f(t,y,v);y+=h*v+h*h*a/2;const a2=f(t+h,y,v+h*a);v+=h*(a+a2)/2}
 else{const k1v=f(t,y,v),k1y=v,k2v=f(t+h/2,y+h/2*k1y,v+h/2*k1v),k2y=v+h/2*k1v,k3v=f(t+h/2,y+h/2*k2y,v+h/2*k2v),k3y=v+h/2*k2v,k4v=f(t+h,y+h*k3y,v+h*k3v),k4y=v+h*k3v;y+=h/6*(k1y+2*k2y+2*k3y+k4y);v+=h/6*(k1v+2*k2v+2*k3v+k4v)}}
 return[y,v]};
const mk2=(m,K)=>function*({f,t0,y0,v0,t1}){for(let k=2;k<=K;k++)yield sysE(f,t0,y0,v0,t1,2**k,m)[0]};
const ODE2={'Euler':mk2('eu',18),'Euler simpléctico':mk2('sy',18),'Velocity-Verlet':mk2('vl',16),'RK4':mk2('rk',15),
'RK4+Richardson':function*({f,t0,y0,v0,t1}){for(let k=1;k<=13;k++)yield(16*sysE(f,t0,y0,v0,t1,2**(k+1),'rk')[0]-sysE(f,t0,y0,v0,t1,2**k,'rk')[0])/15}};
/* ===== mínimo ===== */
const MIN={
'Sección áurea':function*({f,a,b}){const g=.6180339887498949;let c=b-g*(b-a),d=a+g*(b-a),fc=f(c),fd=f(d);for(let k=0;k<200;k++){if(fc<fd){b=d;d=c;fd=fc;c=b-g*(b-a);fc=f(c)}else{a=c;c=d;fc=fd;d=a+g*(b-a);fd=f(d)}yield(a+b)/2}},
'Parabólica':function*({f,a,b}){let[x0,x1,x2]=[a,(a+b)/2,b];for(let k=0;k<100;k++){const f0=f(x0),f1=f(x1),f2=f(x2),num=(x1-x0)**2*(f1-f2)-(x1-x2)**2*(f1-f0),den=2*((x1-x0)*(f1-f2)-(x1-x2)*(f1-f0));if(!isFinite(num/den)||Math.abs(den)<1e-300)return;const u=Math.min(b,Math.max(a,x1-num/den));yield u;const arr=[[x0,f0],[x1,f1],[x2,f2],[u,f(u)]].sort((p,q)=>p[1]-q[1]).slice(0,3).sort((p,q)=>p[0]-q[0]);[x0,x1,x2]=arr.map(p=>p[0])}},
'Newton (f′=0)':function*({f,a,b}){let x=(a+b)/2;for(let k=0;k<200;k++){const h=1e-4*(1+Math.abs(x)),d1=(f(x+h)-f(x-h))/(2*h),d2=(f(x+h)-2*f(x)+f(x-h))/(h*h);let s=d2>0?-d1/d2:-Math.sign(d1)*.1*(b-a);x=Math.min(b,Math.max(a,x+s));yield x}},
'Descenso de gradiente':function*({f,a,b}){let x=(a+b)/2,lr=(b-a)/10;for(let k=0;k<3000;k++){const h=1e-5*(1+Math.abs(x)),d=(f(x+h)-f(x-h))/(2*h);let y=Math.min(b,Math.max(a,x-lr*d));while(f(y)>f(x)&&lr>1e-15){lr/=2;y=Math.min(b,Math.max(a,x-lr*d))}x=y;lr*=1.5;yield x}},
'Aleatoria adaptativa':function*({f,a,b}){let x=rnd(a,b),best=f(x),r=(b-a)/2;for(let k=0;k<6000;k++){const y=Math.min(b,Math.max(a,x+(Math.random()*2-1)*r)),fy=f(y);if(fy<best){best=fy;x=y}else r=Math.max(r*.995,1e-14);yield x}}};
/* ===== límite ===== */
const LIM={
'Lateral derecha':function*({f,x0}){for(let k=1;k<=45;k++)yield f(x0+.1*.5**k)},
'Lateral izquierda':function*({f,x0}){for(let k=1;k<=45;k++)yield f(x0-.1*.5**k)},
'Bilateral':function*({f,x0}){for(let k=1;k<=45;k++){const h=.1*.5**k;yield(f(x0+h)+f(x0-h))/2}},
'Richardson bilateral':function*({f,x0}){const T=[];for(let i=0;i<10;i++){const h=.1/2**i;T[i]=[(f(x0+h)+f(x0-h))/2];for(let j=1;j<=i;j++)T[i][j]=(4**j*T[i][j-1]-T[i-1][j-1])/(4**j-1);yield T[i][i]}}};
/* ===== series ===== */
function*ps(P,m){const{f,a,b}=P,inf=!isFinite(b);for(let k=4;k<=(inf?20:6);k++){const N=inf?2**k:b-a+1;let S=0,c=0;if(m=='r'){for(let n=a+N-1;n>=a;n--)S+=f(n)}else for(let n=a;n<a+N;n++){if(m=='k'){const y=f(n)-c,t=S+y;c=(t-S)-y;S=t}else S+=f(n)}yield S}}
const SER={'Suma directa':P=>ps(P,'d'),'Suma inversa':P=>ps(P,'r'),'Kahan':P=>ps(P,'k'),
'Aitken Δ²':function*(P){const h=[];for(const S of ps(P,'k')){h.push(S);if(h.length>=3){const v=aitken(h);if(isFinite(v))yield v}}},
'Richardson 1/N':function*(P){const h=[];for(const S of ps(P,'k')){h.push(S);if(h.length>=2)yield 2*h[h.length-1]-h[h.length-2]}}};
/* ===== definición de tipos ===== */
const vecErr=(x,t)=>Math.max(...t.map((q,i)=>Math.abs(q-x[i])));
const allRoots=(g,a,b)=>{const R=[];let px=a,pf=g(a);for(let k=1;k<=4000;k++){const x=a+(b-a)*k/4000,fx=g(x);if(isFinite(pf)&&isFinite(fx)&&pf*fx<=0){let l=px,r=x,fl=pf;for(let i=0;i<70;i++){const m=(l+r)/2,fm=g(m);if(fl*fm<=0)r=m;else{l=m;fl=fm}}const z=(l+r)/2;if(Math.abs(g(z))<1e-6&&!R.some(q=>Math.abs(q-z)<1e-9))R.push(z)}px=x;pf=fx}return R};
const cl2=c=>['≥2 métodos independientes coinciden',c.cl>=2];
const TYPES={
solve:{n:'Ecuación (todas las raíces)',ag:SOLVE,kind:'res',rt:1e-10,lab:'x ≈',ev:(P,x)=>Math.abs(P.g(x)),pol:(P,x)=>x-P.g(x)/dg(P.g,x),cal:{P:{g:x=>x*x-2,a:0,b:3},t:Math.SQRT2},
 tests:(P,v,c)=>{const g=P.g,d=1e-7*(1+Math.abs(v)),dd=dg(g,v);return[['residual < 1e-9',Math.abs(g(v))<1e-9],['cambio de signo en el entorno',g(v-d)*g(v+d)<=0],['paso de Newton estable',g(v)==0||Math.abs(g(v)/dd)<1e-8*(1+Math.abs(v))],['Newton desde un punto cercano vuelve a la misma raíz',(()=>{let x=v+1e-3*(1+Math.abs(v));for(let i=0;i<30;i++){const d=dg(g,x);if(!d)break;x-=g(x)/d}return Math.abs(x-v)<1e-6*(1+Math.abs(v))})()]]},
 extra:(P,v)=>{const R=allRoots(P.g,P.a,P.b);return`Raíces en [${P.a}, ${P.b}] (barrido de 4000 puntos + refinado): <b>${R.length?R.slice(0,8).map(x=>x.toPrecision(10)).join(' · '):'ninguna con cambio de signo'}</b>${R.length>8?' …':''}`}},
lin:{n:'Sistema lineal Ax=b',ag:LIN,kind:'res',rt:1e-11,lab:'x =',ev:(P,x)=>Math.max(...mv(P.A,x).map((v,i)=>Math.abs(v-P.b[i]))),pol:(P,x)=>{const r=mv(P.A,x).map((v,i)=>P.b[i]-v),d=gsolve(P.A,r);return d?x.map((v,i)=>v+d[i]):null},cal:{P:{A:[[2,1],[1,3]],b:[3,5]},t:[.8,1.4]},
 tests:(P,v,c)=>{const g=gsolve(P.A,P.b),b2=P.b.map(q=>q+1e-8),x2=gsolve(P.A,b2),nx=Math.hypot(...v)||1;return[['residual ∞ < 1e-9',c.r<1e-9],['Gauss independiente coincide',!!g&&vecErr(v,g)<1e-8*(1+nx)],['bien condicionado (δb=1e-8)',!!x2&&Math.hypot(...x2.map((q,i)=>q-v[i]))/nx<1e-4],['Gauss-Jordan independiente coincide',(()=>{let x=null;for(const q of gj(P,1))if(q)x=q;return!!x&&vecErr(v,x)<1e-8*(1+nx)})()]]}},
integral:{n:'Integral definida',ag:INT,kind:'cons',tol:1e-10,lab:'∫ ≈',cal:{P:{f:x=>x*x,a:0,b:3},t:9},
 tests:(P,v,c)=>{const m=(P.a+P.b)/2,e=1e-8*(1+Math.abs(v));return[cl2(c),['aditividad ∫[a,m]+∫[m,b]',Math.abs(glq(P.f,P.a,m,300)+glq(P.f,m,P.b,300)-v)<e],['Gauss-Legendre 4000 paneles coincide',Math.abs(glq(P.f,P.a,P.b,800)-v)<e],['valor finito',isFinite(v)]]}},
deriv:{n:'Derivada en un punto',ag:DER,kind:'cons',tol:1e-8,lab:'f′(x₀) ≈',cal:{P:{f:x=>x**3,x0:2},t:12},
 tests:(P,v,c)=>[cl2(c),['5 puntos con otro h coincide',Math.abs(D5(P.f,P.x0,1e-2)-v)<1e-6*(1+Math.abs(v))],['continuidad de f′ en x₀±δ',Math.abs(D5(P.f,P.x0+1e-4,1e-3)-D5(P.f,P.x0-1e-4,1e-3))<1e-2*(1+Math.abs(v))]]},
ode:{n:'EDO 1.er orden y′=f(t,y)',ag:ODE,kind:'cons',tol:1e-9,lab:'y(t₁) ≈',cal:{P:{f:(t,y)=>y,t0:0,y0:1,t1:1},t:Math.E},
 tests:(P,v,c)=>[cl2(c),['RK4 con n y 2n coincide',Math.abs(rk(P.f,P.t0,P.y0,P.t1,4096)-rk(P.f,P.t0,P.y0,P.t1,8192))<1e-8*(1+Math.abs(v))],['reversibilidad temporal',Math.abs(rk(P.f,P.t1,v,P.t0,8192)-P.y0)<1e-6*(1+Math.abs(P.y0))]]},
ode2:{n:"EDO 2.º orden y″=f(t,y,v)",ag:ODE2,kind:'cons',tol:1e-9,lab:'y(t₁) ≈',cal:{P:{f:(t,y,v)=>-y,t0:0,y0:1,v0:0,t1:1},t:Math.cos(1)},
 tests:(P,v,c)=>{const a=sysE(P.f,P.t0,P.y0,P.v0,P.t1,4096,'rk'),b=sysE(P.f,P.t0,P.y0,P.v0,P.t1,8192,'rk'),r=sysE(P.f,P.t1,b[0],b[1],P.t0,8192,'rk');return[cl2(c),['RK4 con n y 2n coincide',Math.abs(a[0]-b[0])<1e-8*(1+Math.abs(v))],['reversibilidad temporal (y y v)',Math.abs(r[0]-P.y0)<1e-6*(1+Math.abs(P.y0))&&Math.abs(r[1]-P.v0)<1e-5*(1+Math.abs(P.v0))]]},
 extra:(P,v)=>{const b=sysE(P.f,P.t0,P.y0,P.v0,P.t1,8192,'rk');return`v(t₁) ≈ <b>${b[1].toPrecision(10)}</b>`}},
min:{n:'Mínimo de f(x) en [a,b]',ag:MIN,kind:'cons',tol:1e-7,lab:'x* ≈',cal:{P:{f:x=>(x-1.5)**2+1,a:-5,b:5},t:1.5},
 tests:(P,v,c)=>{const f=P.f,d=1e-4*(1+Math.abs(v)),h=1e-4*(1+Math.abs(v)),bd=v<=P.a+1e-9||v>=P.b-1e-9;return[cl2(c),['f(x*) ≤ f(x*±δ)',f(v)<=f(v-d)+1e-12&&f(v)<=f(v+d)+1e-12],['|f′(x*)|≈0 o en el borde',bd||Math.abs((f(v+h)-f(v-h))/(2*h))<1e-4*(1+Math.abs(f(v)))],['f″(x*) ≥ 0',(f(v+h)-2*f(v)+f(v-h))/(h*h)>=-1e-5]]},
 extra:(P,v)=>`f(x*) = <b>${P.f(v).toPrecision(10)}</b> (mínimo local encontrado; los agentes parten de puntos distintos)`},
limit:{n:'Límite x→x₀',ag:LIM,kind:'cons',tol:1e-8,lab:'lím ≈',cal:{P:{f:x=>Math.sin(x)/x,x0:0},t:1},
 tests:(P,v,c)=>[cl2(c),['límites laterales coinciden',Math.abs(P.f(P.x0+1e-4)-P.f(P.x0-1e-4))<1e-3*(1+Math.abs(v))],['f(x₀±h) se aproxima al valor',Math.abs(P.f(P.x0+1e-4)-v)<1e-3*(1+Math.abs(v))]]},
series:{n:'Serie Σ f(n)',ag:SER,kind:'cons',tol:1e-8,lab:'Σ ≈',cal:{P:{f:n=>1/(n*(n+1)),a:1,b:Infinity},t:1},
 tests:(P,v,c)=>{const inf=!isFinite(P.b),h=[];if(inf)for(const S of ps(P,'k'))h.push(S);return[cl2(c),['Aitken independiente coincide',!inf||Math.abs(aitken(h.slice(-3))-v)<1e-6*(1+Math.abs(v))],['término general → 0',!inf||Math.abs(P.f(1e6))<1e-3]]}}};
const CFG={solve:['Ecuación (usa =; ej: cos(x)=x)',['a (inicio)','b (fin)','',''],[-10,10]],lin:['Matriz A (filas con ; y columnas con ,)',['b (vector con comas)','','',''],['']],integral:['f(x)',['a','b','',''],[0,'pi']],deriv:['f(x)',['x₀','','',''],[1]],ode:['f(t,y) en y′=f(t,y)',['t₀','y(t₀)','t₁',''],[0,0,3]],ode2:['f(t,y,v) en y″=f(t,y,v)',['t₀','y(t₀)','v(t₀)','t₁'],[0,1,0,3]],min:['f(x)',['a','b','',''],[-5,5]],limit:['f(x)',['x₀','','',''],[0]],series:['f(n)',['a (n inicial)','b (n final o inf)','',''],[1,'inf']]};
const PRE=[['Raíz: x³−2x−5=0','solve','x^3-2*x-5=0',-10,10],['Raíces múltiples: sin(x)=x/3','solve','sin(x)=x/3',-10,10],['Física: tiempo de vuelo (v₀=20, 45°)','solve','20*sin(pi/4)*x-4.905*x^2=0',.1,10],['Sistema 3×3','lin','2,1,-1;-3,-1,2;-2,1,2','8,-11,-3'],['Integral: sin²(x), 0→π','integral','sin(x)^2',0,'pi'],['Integral: e^(−x²), 0→2','integral','exp(-x^2)',0,2],['Derivada: x²·sin(x) en 1','deriv','x^2*sin(x)',1],['EDO física: caída con rozamiento v′=9.81−0.2v²','ode','9.81-0.2*y^2',0,0,3],['Oscilador: y″=−4y, y(0)=1','ode2','-4*y',0,1,0,3],['Péndulo no lineal (L=1.5 m)','ode2','-9.81/1.5*sin(y)',0,1,0,2],['Oscilador amortiguado','ode2','-4*y-0.5*v',0,1,0,5],['Mínimo: x²+3sin(x)','min','x^2+3*sin(x)',-3,3],['Física: equilibrio resorte (k=50, m=2)','min','0.5*50*x^2-9.81*2*x',-2,2],['Límite: sin(x)/x en 0','limit','sin(x)/x',0],['Límite: (1+x)^(1/x) en 0','limit','(1+x)^(1/x)',0],['Serie: Σ1/n² (n=1→∞)','series','1/n^2',1,'inf'],['Serie: Σ1/(n(n+1))','series','1/(n*(n+1))',1,'inf']];
/* ===== memoria ===== */
let S={};try{S=JSON.parse(localStorage.getItem('arena_v3')||'{}')}catch(e){}
const sv=()=>{try{localStorage.setItem('arena_v3',JSON.stringify(S))}catch(e){}};
const mem=(ty,n)=>S[ty+':'+n]||(S[ty+':'+n]={bank:0,cp:[],cw:0,sw:0,rel:1});
const lvl=m=>1+Math.min(5,Math.floor(m.bank/50));
const fmt=v=>Array.isArray(v)?'['+v.map(q=>q.toPrecision(10)).join(', ')+']':Number(v).toPrecision(12);
function calib(L,G){try{const g=G(L.cal.P),t=L.cal.t;for(let k=0;k<40;k++){const r=g.next();if(r.done)return false;const v=r.value;if(v===null)continue;const d=Array.isArray(t)?vecErr(v,t):Math.abs(t-v);if(d<=5e-3*(1+(Array.isArray(t)?0:Math.abs(t))))return true}}catch(e){}return false}
/* ===== UI ===== */
const ty=$('ty');
function cfg(){const c=CFG[ty.value];$('el').textContent=c[0];[1,2,3,4].forEach(i=>{$('l'+i).textContent=c[1][i-1];$('p'+i).style.visibility=c[1][i-1]?'visible':'hidden';$('p'+i).value=c[2][i-1]??''})}
function pre(i){const p=PRE[i];ty.value=p[1];cfg();$('ex').value=p[2];[1,2,3,4].forEach(j=>$('p'+j).value=p[2+j]??'')}

let token=0;$('st').onclick=()=>{token++};
$('rs').onclick=()=>{S={};sv();$('ags').innerHTML='';$('sol').style.display='none'};
$('go').onclick=run;$('ex').onkeydown=e=>{if(e.key=='Enter')run()};
const err=m=>{$('sol').style.display='block';$('sol').innerHTML='<b>Error:</b> '+m};
const pv=i=>{const s=$('p'+i).value.trim();if(/^[+-]?inf/i.test(s))return s[0]=='-'?-Infinity:Infinity;return compile(s)({})};

async function run(){
 const my=++token,type=ty.value,L=TYPES[type],V={},src=$('ex').value,P={};
 try{
  const c4=()=>{const c=compile(src);return[c,V]};
  if(type=='solve'){const q=src.split('=');const l=compile(q[0]),r=q[1]&&q[1].trim()?compile(q[1]):()=>0;P.g=x=>{V.x=x;return l(V)-r(V)};P.a=pv(1);P.b=pv(2);if(!(P.a<P.b))throw Error('a debe ser menor que b');P.g((P.a+P.b)/2)}
  else if(type=='orbit'){P.mu=compile(src)({});P.r0=pv(1);P.v0=pv(2);P.t1=pv(3);if(!(P.mu>0&&P.r0>0&&P.t1>0))throw Error('mu, r0 y t₁ deben ser positivos')}
  else if(type=='lin'){P.A=src.split(';').map(r=>r.split(',').map(s=>compile(s)({})));P.b=$('p1').value.split(',').map(s=>compile(s)({}));if(P.A.some(r=>r.length!=P.A.length)||P.b.length!=P.A.length)throw Error('A debe ser cuadrada y b del mismo tamaño')}
  else{const c=compile(src);
   if(type=='ode'){P.f=(t,y)=>{V.t=V.x=t;V.y=y;return c(V)};P.t0=pv(1);P.y0=pv(2);P.t1=pv(3);if(!(P.t1>P.t0))throw Error('t₁ debe ser mayor que t₀');P.f(P.t0,P.y0)}
   else if(type=='ode2'){P.f=(t,y,v)=>{V.t=V.x=t;V.y=y;V.v=v;return c(V)};P.t0=pv(1);P.y0=pv(2);P.v0=pv(3);P.t1=pv(4);if(!(P.t1>P.t0))throw Error('t₁ debe ser mayor que t₀');P.f(P.t0,P.y0,P.v0)}
   else if(type=='series'){P.f=n=>{V.n=V.x=n;return c(V)};P.a=pv(1);P.b=pv(2);if(!(P.b>=P.a))throw Error('b debe ser ≥ a');P.f(P.a)}
   else{P.f=x=>{V.x=V.t=x;return c(V)};P.a=pv(1);P.x0=P.a;if(type=='integral'||type=='min'){P.b=pv(2);if(!(P.a<P.b))throw Error('a debe ser menor que b')}if(isNaN(P.a))throw Error('parámetro inválido');if(type=='limit')P.f(P.x0+1e-3);else P.f(P.a)}}
 }catch(e){return err(e.message)}
 window.onStart&&onStart({type,P});
 const md=$('md').value,res_=L.kind=='res';
 const ag=Object.entries(L.ag).map(([name,G])=>{const m=mem(type,name);let coop=md=='coop'?1:md=='solo'?0:(Math.random()<(m.cw+1)/(m.cw+m.sw+2));
  const apt=calib(L,G);m.rel=.7*m.rel+.3*(apt?1:0);
  const e0=(70+Math.min(330,m.bank*.5))*(apt?1:.7);
  return{name,m,gen:G(P),coop,apt,e0,e:e0,raw:[],est:NaN,delta:Infinity,conv:false,alive:true,steps:0,lv:lvl(m)}});
 const t0=performance.now(),tol=L.tol;let tick=0,res=null,BB={r:Infinity,est:null};
 $('sol').style.display='none';
 const upd=c=>{if(res_&&nv(c)){const r=L.ev(P,c);if(r<BB.r)BB={r,est:c};return r}return Infinity};
 const step=a=>{
  if(a.e<=0){a.alive=false;return}
  a.e--;a.steps++;let r;try{r=a.gen.next()}catch(e){a.alive=false;return}
  if(r.done){a.alive=false;return}if(r.value===null)return;
  const est=r.value;if(!nv(est))return;a.raw.push(est);
  let acc=NaN;if((a.coop||a.m.cp.length)&&!Array.isArray(est)&&a.raw.length>=3)acc=aitken(a.raw);
  if(res_){for(const c of[est,acc])if(nv(c)&&upd(c)<=L.rt){a.est=c;return a}a.est=est;return}
  let ch=est;const prev=a.raw[a.raw.length-2];
  if(isFinite(acc)&&prev!==undefined&&Math.abs(acc-est)<10*Math.abs(est-prev)+1e-300)ch=acc;
  a.delta=isNaN(a.est)?Infinity:Math.abs(ch-a.est);a.est=ch;a.conv=a.raw.length>=3&&a.delta<=tol*(1+Math.abs(ch))};
 const cons=()=>{const c=ag.filter(a=>a.conv&&isFinite(a.est));for(const a of c){const cl=c.filter(b=>Math.abs(b.est-a.est)<=tol*10*(1+Math.abs(a.est)));if(cl.length>=2)return cl}return null};
 loop:while(my==token){
  tick++;
  for(const a of ag){for(let k=0;k<a.lv+(a.coop?0:1)&&a.alive;k++){const w=step(a);
   if(w){res={win:[a],val:a.est};break loop}
   if(!res_){const cl=cons();if(cl){res={win:cl};break loop}}}}
  if(res_&&BB.est)for(const a of ag){if(!a.coop||!a.alive||a.e<=0)continue;a.e--;a.steps++;let c=null;try{c=L.pol(P,BB.est)}catch(e){}if(nv(c)&&upd(c)<=L.rt){a.est=c;res={win:[a],val:c,team:1};break loop}}
  render(ag,tick,null);window.onTick&&onTick({type,P,ag,tick});
  if(!ag.some(a=>a.alive))break;
  await new Promise(r=>setTimeout(r,30));
 }
 if(my!=token)return;
 const ms=performance.now()-t0;let val,wins=[],p,tests=[],ex='';
 if(res){
  wins=res.win;
  let ctx={cl:wins.length,agree:0,r:0};
  if(res_){val=res.val;ctx.r=L.ev(P,val);const ok=ag.filter(a=>!wins.includes(a)&&nv(a.est)&&L.ev(P,a.est)<1e-6);ctx.agree=ok.length;ok.forEach(a=>wins.push(a));
   const k0=Math.min(3,-Math.log10(ctx.r+1e-17)/5)+.5*ok.length;ctx.k=k0}
  else{wins=res.win.slice().sort((a,b)=>a.delta-b.delta);val=wins[0].est;const sp=Math.max(...wins.map(a=>Math.abs(a.est-val)));ctx.k=Math.min(3,-Math.log10(sp/(1+Math.abs(val))+1e-17)/4)+(wins.length-1)*.6;ctx.sp=sp}
  try{tests=L.tests(P,val,ctx)}catch(e){tests=[['micro pruebas ejecutables',false]]}
  const f=tests.length?tests.filter(t=>t[1]).length/tests.length:1;
  p=Math.min(.9999,Math.max(.02,(1-10**-(ctx.k+(res_?1.5*f:0)))*(.55+.45*f*f)));
  try{if(L.extra)ex=L.extra(P,val)}catch(e){}
  const w0=wins[0];
  wins.forEach((a,i)=>{a.m.bank+=(i==0?(a.coop?20:35):(a.coop?12:8))+Math.round(10*Math.max(0,a.e)/a.e0)});
  if(w0.coop)w0.m.cw++;else w0.m.sw++;
  const rest=ag.filter(a=>!wins.includes(a)&&nv(a.est)).sort((a,b)=>res_?L.ev(P,a.est)-L.ev(P,b.est):Math.abs(a.est-val)-Math.abs(b.est-val))[0];
  if(rest){if(w0.m.cp.length<3&&!w0.m.cp.includes(rest.name))w0.m.cp.push(rest.name);if(rest.m.cp.length<3&&!rest.m.cp.includes(w0.name))rest.m.cp.push(w0.name)}
  sv();
 }else{const c=ag.filter(a=>nv(a.est)).sort((a,b)=>a.delta-b.delta)[0];val=c?c.est:NaN;p=.05}
 window.afterRun&&afterRun({type,P,L,val,res,ag,wins,p,tests,ms});
 render(ag,tick,wins);
 const s=$('sol');s.style.display='block';
 s.innerHTML=`<b>Solución</b><div class="big">${L.lab} ${nv(val)?fmt(val):'—'}</div>
 <div class="m">Probabilidad de acierto: <b>${(p*100).toFixed(res?2:0)}%</b></div><div class="pb"><i style="width:${p*100}%"></i></div>
 <div class="m" style="margin:8px 0 4px">${res?'Resuelto por: <b>'+wins.map(a=>a.name).join(', ')+'</b>'+(res.team?' (cooperación: pulido del mejor candidato compartido)':'')+' en ':'Sin resolver tras '}${tick} ciclos · ${ms.toFixed(0)} ms · energía gastada ${ag.reduce((x,a)=>x+a.steps,0)}</div>
 ${ex?`<div class="m">${ex}</div>`:''}
 ${tests.length?'<div style="margin-top:8px"><b>Micro pruebas de verificación</b></div>'+tests.map(t=>`<div class="t ${t[1]?'ok':'no'}">${t[1]?'✔':'✘'} ${t[0]}</div>`).join(''):'<div class="m">Los agentes agotaron su energía sin verificación ni consenso. Prueba otro intervalo o expresión.</div>'}
 <div class="m" style="margin-top:6px">p = (1−10^−k)·(0.55+0.45·f²): k según residual/dispersión y acuerdo entre agentes, f = fracción de micro pruebas superadas.</div>`;
}
function render(ag,tick,wins){
 $('tk').textContent='· ciclo '+tick;
 $('ags').innerHTML=ag.map(a=>{const w=wins&&wins.includes(a),pct=Math.max(0,a.e/a.e0*100),e=a.est;
  return`<div class="ag ${w?'win':''}"><span><b>${a.name}</b> ${a.coop?'🤝':'🐺'} <span class="m">nv ${a.lv} · banco ${a.m.bank}⚡ · ${a.apt?'✔ apto':'✘ micro-prueba fallida'}${a.m.cp.length?' · acoplado: '+a.m.cp.join('+'):''}</span></span>
  <span class="m">${wins?(w?'🏆 +'+(wins[0]===a?(a.coop?20:35):(a.coop?12:8))+'⚡':a.alive?'detenido':'sin energía'):a.alive?'activo':'fuera'}</span>
  <div class="bar"><i style="width:${pct}%"></i></div>
  <span class="m" style="grid-column:1/3;font-family:ui-monospace,monospace">${nv(e)&&e!==undefined?'est = '+(Array.isArray(e)?'['+e.slice(0,4).map(q=>q.toPrecision(7)).join(', ')+(e.length>4?'…':'')+']':e.toPrecision(10)):'est = …'} · energía ${Math.max(0,a.e)}/${Math.round(a.e0)}</span></div>`}).join('')}
