/* Visuales: gráfica del problema + iterados en vivo de cada agente + convergencia */
const COL=['#2dd4bf','#f59e0b','#60a5fa','#f472b6','#a3e635','#c084fc','#fb7185','#34d399','#fbbf24','#38bdf8','#e879f9','#f87171','#4ade80','#facc15'];
const PL={k:{}};
const cs=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
function cvs(id){const c=$(id),d=window.devicePixelRatio||1,w=c.clientWidth,h=c.clientHeight;if(c.width!=Math.round(w*d)){c.width=Math.round(w*d);c.height=Math.round(h*d)}const x=c.getContext('2d');x.setTransform(d,0,0,d,0,0);x.clearRect(0,0,w,h);return{x,w,h}}
const nice=r=>{const e=10**Math.floor(Math.log10(r)),f=r/e;return(f<1.5?1:f<3.5?2:f<7.5?5:10)*e};
function mkView(c,x0,x1,y0,y1,eq){const{x,w,h}=c,m={l:50,r:12,t:12,b:26};if(!(x1>x0)){x0-=1;x1+=1}if(!(y1>y0)){y0-=1;y1+=1}
 if(eq){const pw=w-m.l-m.r,ph=h-m.t-m.b,s=Math.min(pw/(x1-x0),ph/(y1-y0)),cx=(x0+x1)/2,cy=(y0+y1)/2;x0=cx-pw/s/2;x1=cx+pw/s/2;y0=cy-ph/s/2;y1=cy+ph/s/2}
 const V={x,w,h,m,x0,x1,y0,y1,X:v=>m.l+(v-x0)/(x1-x0)*(w-m.l-m.r),Y:v=>h-m.b-(v-y0)/(y1-y0)*(h-m.t-m.b)};
 x.save();x.font='11px system-ui';x.fillStyle=cs('--mut');x.strokeStyle=cs('--grid');x.lineWidth=1;
 const sx=nice((x1-x0)/6),sy=nice((y1-y0)/5);
 for(let t=Math.ceil(x0/sx)*sx;t<=x1;t+=sx){x.beginPath();x.moveTo(V.X(t),m.t);x.lineTo(V.X(t),h-m.b);x.stroke();x.textAlign='center';x.fillText(+t.toPrecision(4),V.X(t),h-8)}
 for(let t=Math.ceil(y0/sy)*sy;t<=y1;t+=sy){x.beginPath();x.moveTo(m.l,V.Y(t));x.lineTo(w-m.r,V.Y(t));x.stroke();x.textAlign='right';x.fillText(+t.toPrecision(4),m.l-5,V.Y(t)+4)}
 x.strokeStyle=cs('--mut');if(y0<0&&y1>0){x.beginPath();x.moveTo(m.l,V.Y(0));x.lineTo(w-m.r,V.Y(0));x.stroke()}if(x0<0&&x1>0){x.beginPath();x.moveTo(V.X(0),m.t);x.lineTo(V.X(0),h-m.b);x.stroke()}
 x.restore();x.save();x.beginPath();x.rect(m.l,m.t,w-m.l-m.r,h-m.t-m.b);x.clip();return V}
function path(V,pts,col,wd=2,fill){const x=V.x;x.save();x.strokeStyle=col;x.lineWidth=wd;x.beginPath();let pen=false;for(const[a,b]of pts){if(!isFinite(a)||!isFinite(b)){pen=false;continue}const px=V.X(a),py=V.Y(b);if(!pen){x.moveTo(px,py);pen=true}else x.lineTo(px,py)}x.stroke();if(fill){x.fillStyle=fill;x.fill()}x.restore()}
const dot=(V,a,b,col,r=4,hollow)=>{if(!isFinite(a)||!isFinite(b))return;const x=V.x;x.beginPath();x.arc(V.X(a),V.Y(b),r,0,7);if(hollow){x.strokeStyle=col;x.lineWidth=2;x.stroke()}else{x.fillStyle=col;x.fill()}};
const samp=(f,a,b,n=500)=>{const o=[];for(let i=0;i<=n;i++){const t=a+(b-a)*i/n;let y;try{y=f(t)}catch(e){y=NaN}o.push([t,y])}return o};
function yrange(pts,extra=[]){const ys=pts.map(p=>p[1]).concat(extra).filter(isFinite).sort((a,b)=>a-b);if(!ys.length)return[-1,1];const lo=ys[Math.floor(ys.length*.01)],hi=ys[Math.ceil(ys.length*.99)-1],p=(hi-lo)*.12||1;return[lo-p,hi+p]}
const lab=(V,t,xx,yy,col)=>{V.x.save();V.x.fillStyle=col||cs('--fg');V.x.font='12px system-ui';V.x.fillText(t,xx,yy);V.x.restore()};
PL.title=t=>{$('ptitle').textContent=t};
function rk4Traj(f,t0,y0,t1,n){const o=[[t0,y0]],h=(t1-t0)/n;let y=y0;for(let i=0;i<n;i++){const t=t0+i*h,k1=f(t,y),k2=f(t+h/2,y+h*k1/2),k3=f(t+h/2,y+h*k2/2),k4=f(t+h,y+h*k3);y+=h*(k1+2*k2+2*k3+k4)/6;o.push([t+h,y])}return o}
function drawMain(o){
 const{type,P,ag,val,fin}=o,c=cvs('plot'),ok=fin&&nv(val);let V;
 const iter=(f)=>ag&&ag.forEach((a,i)=>{a.raw.slice(-25).forEach(e=>{if(!Array.isArray(e))dot(V,e,f(e),COL[i%COL.length],3)})});
 try{switch(type){
 case'solve':{const p=samp(P.g,P.a,P.b),[y0,y1]=yrange(p);V=mkView(c,P.a,P.b,y0,y1);path(V,p,cs('--ac'),2.2);iter(P.g);if(fin){try{allRoots(P.g,P.a,P.b).forEach(r=>dot(V,r,0,'#fff',6,1))}catch(e){}if(ok)dot(V,val,0,cs('--ac'),6)}PL.title('Gráfica de f(x)=lado izq − lado der · ○ raíces · puntos = iterados de cada agente');break}
 case'integral':{const w=(P.b-P.a)*.08,p=samp(P.f,P.a-w,P.b+w),[y0,y1]=yrange(p,[0]);V=mkView(c,P.a-w,P.b+w,Math.min(y0,0),y1);const q=samp(P.f,P.a,P.b,300);const poly=[[P.a,0],...q,[P.b,0]];const g=V.x;g.save();g.beginPath();poly.forEach(([a,b],i)=>{if(isFinite(b))i?g.lineTo(V.X(a),V.Y(b)):g.moveTo(V.X(a),V.Y(b))});g.closePath();g.fillStyle=cs('--ac')+'55';g.fill();g.restore();path(V,p,cs('--ac'),2.2);if(ok)lab(V,'∫ = '+val.toPrecision(8),V.m.l+10,V.m.t+18);PL.title('Área bajo f(x) entre a y b');break}
 case'deriv':{const s=Math.max(1,Math.abs(P.x0)),p=samp(P.f,P.x0-s,P.x0+s),[y0,y1]=yrange(p);V=mkView(c,P.x0-s,P.x0+s,y0,y1);path(V,p,cs('--ac'),2.2);const f0=P.f(P.x0);if(ok)path(V,[[P.x0-s,f0-val*s],[P.x0+s,f0+val*s]],'#f59e0b',1.8);dot(V,P.x0,f0,'#fff',5);if(ok)lab(V,"pendiente f′ = "+val.toPrecision(8),V.m.l+10,V.m.t+18);PL.title('Recta tangente en x₀ (su pendiente es la derivada)');break}
 case'min':{const p=samp(P.f,P.a,P.b),[y0,y1]=yrange(p);V=mkView(c,P.a,P.b,y0,y1);path(V,p,cs('--ac'),2.2);iter(P.f);if(ok){dot(V,val,P.f(val),'#fff',7,1);lab(V,'x* = '+val.toPrecision(8),V.m.l+10,V.m.t+18)}PL.title('Función a minimizar · puntos = iterados de cada agente');break}
 case'limit':{const p=samp(P.f,P.x0-1,P.x0+1,600),[y0,y1]=yrange(p,ok?[val]:[]);V=mkView(c,P.x0-1,P.x0+1,y0,y1);path(V,p,cs('--ac'),2.2);if(ok){dot(V,P.x0,val,'#f59e0b',6,1);lab(V,'límite = '+val.toPrecision(8),V.m.l+10,V.m.t+18)}PL.title('f(x) cerca de x₀ · círculo = valor del límite');break}
 case'series':{const N=Math.min(60,isFinite(P.b)?P.b-P.a+1:60),pts=[];let S=0;for(let i=0;i<N;i++){S+=P.f(P.a+i);pts.push([P.a+i,S])}const[y0,y1]=yrange(pts,ok?[val]:[]);V=mkView(c,P.a,P.a+N-1,y0,y1);if(ok)path(V,[[P.a,val],[P.a+N-1,val]],'#f59e0b',1.5);path(V,pts,cs('--ac'),2);pts.forEach(q=>dot(V,q[0],q[1],cs('--ac'),2.5));if(ok)lab(V,'Σ = '+val.toPrecision(10),V.m.l+10,V.m.t+18);PL.title('Sumas parciales Sₙ → valor de la serie (línea ámbar)');break}
 case'ode':{const k=PL.k.ode||(PL.k.ode=rk4Traj(P.f,P.t0,P.y0,P.t1,600)),[y0,y1]=yrange(k);V=mkView(c,P.t0,P.t1,y0,y1);path(V,k,cs('--ac'),2.2);if(ok)dot(V,P.t1,val,'#fff',6,1);PL.title('Solución y(t) de la EDO (referencia RK4)');break}
 case'ode2':{let k=PL.k.o2;if(!k){k=PL.k.o2=[];let y=P.y0,v=P.v0;const n=500,h=(P.t1-P.t0)/n;for(let i=0;i<=n;i++){const t=P.t0+i*h;k.push([t,y,v]);if(i<n)[y,v]=sysE(P.f,t,y,v,t+h,1,'rk')}}const A=k.map(q=>[q[0],q[1]]),B=k.map(q=>[q[0],q[2]]),[y0,y1]=yrange(A.concat(B));V=mkView(c,P.t0,P.t1,y0,y1);path(V,A,cs('--ac'),2.2);path(V,B,'#f59e0b',1.6);if(ok)dot(V,P.t1,val,'#fff',6,1);lab(V,'— y(t)',V.m.l+10,V.m.t+16,cs('--ac'));lab(V,'— v(t)',V.m.l+70,V.m.t+16,'#f59e0b');PL.title('Posición y(t) y velocidad v(t)');break}
 case'lin':{const n=P.b.length;if(n==2){const x=ok?val:[0,0],R=Math.max(2,Math.abs(x[0])*2,Math.abs(x[1])*2);V=mkView(c,x[0]-R,x[0]+R,x[1]-R,x[1]+R,1);P.A.forEach((r,i)=>{const f=u=>r[1]?(P.b[i]-r[0]*u)/r[1]:NaN;const pts=r[1]?samp(f,x[0]-R,x[0]+R,2):[[P.b[i]/r[0],x[1]-R],[P.b[i]/r[0],x[1]+R]];path(V,pts,COL[i],2)});if(ok)dot(V,x[0],x[1],'#fff',6,1);PL.title('Cada ecuación es una recta; la solución es su intersección')}
  else{V=mkView(c,-.5,n-.5,Math.min(0,...(ok?val:[0])),Math.max(1e-9,...(ok?val:[1])));const g=V.x;(ok?val:[]).forEach((v,i)=>{g.fillStyle=COL[i%COL.length];g.fillRect(V.X(i-.35),Math.min(V.Y(v),V.Y(0)),V.X(i+.35)-V.X(i-.35),Math.abs(V.Y(v)-V.Y(0)))});PL.title('Componentes de la solución x₁…xₙ')}break}
 case'orbit':{let k=PL.k.orb;if(!k){k=PL.k.orb=[];let s=[P.r0,0,0,P.v0];const n=1500,h=P.t1/n;for(let i=0;i<=n;i++){k.push(s);s=orbStep(P.mu,s,h,'rk')}}const xs=k.map(q=>q[0]),ys=k.map(q=>q[1]),R=Math.max(...xs.map(Math.abs),...ys.map(Math.abs))*1.15;V=mkView(c,-R,R,-R,R,1);path(V,k.map(q=>[q[0],q[1]]),cs('--ac'),2);dot(V,0,0,'#60a5fa',7);dot(V,P.r0,0,'#f59e0b',5);const e=k[k.length-1];dot(V,e[0],e[1],'#fff',6,1);lab(V,'● cuerpo central  ● inicio  ○ final',V.m.l+10,V.m.t+16);PL.title('Trayectoria orbital (referencia RK4, ejes iguales)');break}
 }}catch(e){}
 if(V)V.x.restore()}
function drawConv(ag,val){const c=cvs('conv');let mx=3;const E=(e)=>Array.isArray(e)?vecErr(e,val)/(1+Math.max(...val.map(Math.abs))):Math.abs(e-val)/(1+Math.abs(val));
 ag.forEach(a=>mx=Math.max(mx,a.raw.length));const V=mkView(c,0,mx,-16,1);
 ag.forEach((a,i)=>{const pts=a.raw.map((e,k)=>[k+1,Math.log10(Math.max(1e-16,E(e)))]);path(V,pts,COL[i%COL.length],1.8)});V.x.restore();
 const g=c.x;g.fillStyle=cs('--mut');g.font='11px system-ui';g.fillText('log₁₀ error relativo',54,24)}
window.onStart=()=>{PL.k={}};
window.onTick=o=>{if(o.tick<4||o.tick%2==0){drawMain({...o,fin:false});$('legend').innerHTML=o.ag.map((a,i)=>`<span><i style="background:${COL[i%COL.length]}"></i>${a.name}</span>`).join('')}};
window.afterRun=o=>{const{type,P,ag,val,res}=o;drawMain({type,P,ag,val,fin:!!res});if(res&&nv(val))drawConv(ag,val);else cvs('conv')};
