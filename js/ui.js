/* Interfaz: secciones Matemáticas / Física / Aeroespacial, conceptos, ranking, calculadora */
const DOM={math:{n:'Matemáticas',t:['solve','lin','integral','deriv','limit','series','min']},phys:{n:'Física',t:['ode','ode2','min','solve','integral','deriv']},aero:{n:'Aeroespacial',t:['orbit','ode2','solve','ode','integral']}};
const pdom=l=>/^Aero/.test(l)?'aero':/^(Física|EDO física|Oscilador|Péndulo)/.test(l)?'phys':'math';
const CON={
solve:{t:'Búsqueda de raíces',f:'f(x) = 0  ·  Newton: xₙ₊₁ = xₙ − f(xₙ)/f′(xₙ)',d:'Un cambio de signo garantiza una raíz (Bolzano): bisección y falsa posición la atrapan sin fallar; Newton, Halley y Steffensen convergen mucho más rápido pero pueden perderse. Los híbridos combinan ambas virtudes.'},
lin:{t:'Sistemas lineales',f:'A·x = b  ·  residual r = b − A·x',d:'Los métodos directos (Gauss, Jordan) eliminan variables; los iterativos (Jacobi, Seidel, SOR, gradiente conjugado, Kaczmarz) refinan una estimación. El refinamiento iterativo corrige el error de redondeo de una solución directa.'},
integral:{t:'Integración numérica',f:'∫ₐᵇ f(x) dx ≈ Σ wᵢ f(xᵢ)',d:'Trapecio y Simpson aproximan por polinomios; Romberg extrapola; Gauss-Legendre elige nodos óptimos; Tanh-Sinh domina con singularidades en los extremos; Monte Carlo y cuasi-Monte Carlo muestrean al azar o con secuencias de baja discrepancia.'},
deriv:{t:'Derivación numérica',f:"f′(x) ≈ [f(x+h) − f(x−h)] / 2h",d:'Reducir h baja el error de truncamiento pero sube el de redondeo: hay un h óptimo. Richardson combina varios h para cancelar términos de error.'},
ode:{t:'Ecuaciones diferenciales (1.er orden)',f:"y′ = f(t, y)  ·  RK4: y += h/6·(k₁+2k₂+2k₃+k₄)",d:'Euler es de orden 1, Heun/punto medio de orden 2, RK4 de orden 4. Dormand-Prince ajusta el paso solo según el error estimado.'},
ode2:{t:'Dinámica de 2.º orden',f:"y″ = f(t, y, y′)  ⇒  sistema (y, v)",d:'Toda ecuación de Newton m·a=F es de este tipo. Los integradores simplécticos conservan mejor la energía en oscilaciones largas; RK4 es más preciso a corto plazo.'},
min:{t:'Optimización',f:"f′(x*) = 0, f″(x*) > 0",d:'Sección áurea y parabólica son seguras; Newton es rápido; recocido simulado y enjambre exploran para evitar mínimos locales.'},
limit:{t:'Límites',f:'lím f(x) cuando x→x₀',d:'Las secuencias h→0 se aceleran con Aitken o Neville/Richardson para extrapolar el valor aunque f no esté definida en x₀.'},
series:{t:'Series infinitas',f:'S = Σ f(n)  ·  Aitken: S ≈ Sₙ − (ΔSₙ)²/Δ²Sₙ',d:'Las sumas parciales convergen lento; Aitken, Richardson y Wynn-ε aceleran. Kahan reduce el error de redondeo al sumar.'},
orbit:{t:'Mecánica orbital (dos cuerpos)',f:'r̈ = −μ·r/|r|³  ·  E = v²/2 − μ/r',d:'La energía y el momento angular se conservan. Los integradores simplécticos (Leapfrog, Yoshida) casi no derivan en energía; Euler hace espiralar la órbita.'}};
const GLO={math:[['Teorema de Bolzano','Si f es continua y cambia de signo en [a,b], hay una raíz dentro.'],['Orden de convergencia','Un método de orden p reduce el error como eₙ₊₁ ≈ C·eₙᵖ.'],['Extrapolación','Combinar resultados con distintos pasos para cancelar términos de error.']],
phys:[['Conservación de energía','En un sistema sin rozamiento E = T + V es constante: sirve de micro prueba.'],['Oscilador armónico','y″ = −ω²y → y = A·cos(ωt); periodo T = 2π/ω.'],['Velocidad terminal','Con arrastre cuadrático v′ = g − kv², v → √(g/k).'],['Reversibilidad','Las leyes de Newton sin disipación son simétricas en el tiempo.']],
aero:[['Ecuación de Kepler','M = E − e·sinE relaciona el tiempo con la posición en una elipse.'],['Vis-viva','v² = μ(2/r − 1/a) da la velocidad en cualquier punto de la órbita.'],['Tsiolkovsky','Δv = vₑ·ln(m₀/m_f), con vₑ = Isp·g₀.'],['Transferencia de Hohmann','Dos impulsos tangenciales: la transferencia de menor Δv entre órbitas circulares coplanares.']]};
let dom='math';
function concept(){const c=CON[ty.value]||{t:'',f:'',d:''};$('concept').innerHTML=`<b>${c.t}</b><div class="fml">${c.f}</div><p>${c.d}</p><details><summary>Glosario · ${DOM[dom].n}</summary>${GLO[dom].map(g=>`<p><b>${g[0]}.</b> ${g[1]}</p>`).join('')}</details>`}
function setDom(d,keep){dom=d;document.body.dataset.d=d;document.querySelectorAll('#tabs button').forEach(b=>b.classList.toggle('on',b.dataset.d==d));
 ty.innerHTML='';DOM[d].t.forEach(k=>ty.add(new Option(TYPES[k].n,k)));const sel=$('pre');sel.innerHTML='';
 PRE.forEach((p,i)=>{if(pdom(p[0])==d&&DOM[d].t.includes(p[1]))sel.add(new Option(p[0],i))});
 $('calc').style.display=d=='aero'?'block':'none';if(!keep&&sel.options.length){sel.value=sel.options[0].value;pre(+sel.value)}concept()}
document.querySelectorAll('#tabs button').forEach(b=>b.onclick=()=>setDom(b.dataset.d));
$('pre').onchange=e=>{pre(+e.target.value);concept()};ty.onchange=()=>{cfg();$('ex').value='';concept()};
$('th').onclick=()=>{const l=document.documentElement.dataset.theme=='light'?'dark':'light';document.documentElement.dataset.theme=l;try{localStorage.setItem('arena_theme',l)}catch(e){}};
try{const t=localStorage.getItem('arena_theme');if(t)document.documentElement.dataset.theme=t}catch(e){}
/* ranking global */
function rank(){const r=Object.entries(S).map(([k,m])=>({k,m})).sort((a,b)=>b.m.bank-a.m.bank).slice(0,8);
 $('rank').innerHTML=r.length?r.map((x,i)=>{const[t,n]=x.k.split(':');return`<div class="rk"><span>${['🥇','🥈','🥉'][i]||i+1+'.'} <b>${n}</b> <span class="m">${TYPES[t]?TYPES[t].n:t}</span></span><span>${x.m.bank}⚡ · nv ${lvl(x.m)} · 🤝${x.m.cw} 🐺${x.m.sw}</span></div>`}).join(''):'<div class="m">Aún no hay victorias. Lanza un problema.</div>'}
const _ar=window.afterRun;window.afterRun=o=>{_ar&&_ar(o);setTimeout(rank,0)};
$('rs').addEventListener('click',()=>setTimeout(rank,0));
/* calculadora orbital */
const BODY={Tierra:[398600.4418,6378.137],Luna:[4902.8,1737.4],Marte:[42828.37,3396.2]};
function calc(){const[mu,R]=BODY[$('cb').value],h1=+$('c1').value,h2=+$('c2').value,r1=R+h1,r2=R+h2,a=(r1+r2)/2;
 const vc=Math.sqrt(mu/r1),T=2*Math.PI*Math.sqrt(r1**3/mu),ve=Math.sqrt(2*mu/r1),d1=Math.sqrt(mu/r1)*(Math.sqrt(2*r2/(r1+r2))-1),d2=Math.sqrt(mu/r2)*(1-Math.sqrt(2*r1/(r1+r2))),tt=Math.PI*Math.sqrt(a**3/mu);
 const isp=+$('ci').value,m0=+$('cm').value,mf=+$('cf').value,dv=isp*9.80665*Math.log(m0/mf);
 $('cout').innerHTML=`<div>v circular: <b>${vc.toFixed(4)} km/s</b> · periodo: <b>${(T/60).toFixed(2)} min</b> · v escape: <b>${ve.toFixed(4)} km/s</b></div>
 <div>Hohmann: Δv₁ <b>${d1.toFixed(4)}</b> + Δv₂ <b>${d2.toFixed(4)}</b> = <b>${Math.abs(d1)+Math.abs(d2)>0?(Math.abs(d1)+Math.abs(d2)).toFixed(4):'—'} km/s</b> · tiempo: <b>${(tt/3600).toFixed(2)} h</b></div>
 <div>Tsiolkovsky: Δv = <b>${isFinite(dv)?(dv/1000).toFixed(4)+' km/s':'—'}</b></div>`}
['cb','c1','c2','ci','cm','cf'].forEach(i=>$(i).oninput=calc);calc();
setDom('math');rank();
window.addEventListener('resize',()=>{if(window.PL&&PL.last)drawMain(PL.last)});
