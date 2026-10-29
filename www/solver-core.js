(function(g){'use strict';const NS=['N','P','K','Ca','Mg','S'],TARGETS={N:150,P:35,K:200,Ca:150,Mg:40,S:50},FERTS=[
{id:'calciumNitrate',name:'نیترات کلسیم',N:15.5,P:0,K:0,Ca:19,Mg:0,S:0,stock:'B',enabled:false},
{id:'potassiumNitrate',name:'نیترات پتاسیم',N:13.7,P:0,K:38.7,Ca:0,Mg:0,S:0,stock:'A',enabled:true},
{id:'mkp',name:'مونوپتاسیم فسفات (MKP)',N:0,P:22.7,K:28.7,Ca:0,Mg:0,S:0,stock:'A',enabled:true},
{id:'potassiumSulfate',name:'سولفات پتاسیم',N:0,P:0,K:44.9,Ca:0,Mg:0,S:18.4,stock:'A',enabled:true},
{id:'magnesiumSulfate',name:'سولفات منیزیم هفت‌آبه',N:0,P:0,K:0,Ca:0,Mg:9.86,S:13,stock:'A',enabled:true},
{id:'map',name:'مونوآمونیوم فسفات (MAP)',N:12.2,P:26.9,K:0,Ca:0,Mg:0,S:0,stock:'A',enabled:false},
{id:'urea',name:'اوره (غیرفعال)',N:46,P:0,K:0,Ca:0,Mg:0,S:0,stock:'A',enabled:false,prohibited:true},
{id:'ammoniumSulfate',name:'سولفات آمونیوم',N:21,P:0,K:0,Ca:0,Mg:0,S:24,stock:'A',enabled:false}];
function coef(p,v,d){p=Number(p);v=Number(v);d=Number(d);return p>=0&&v>0&&d>0?10*p/(v*d):0}
function calcB(t,c){let caC=coef(c.caPct,c.volB,c.denom),nC=coef(c.nPct,c.volB,c.denom);if(caC<=0)return{mass:0,Ca:0,N:0,error:'تنظیمات مخزن B یا درصد Ca معتبر نیست.'};let mass=(+t.Ca||0)/caC;return{mass,Ca:mass*caC,N:mass*nC}}
function solve(t,fs,c){let b=calcB(t,c),act=fs.filter(f=>f.enabled&&f.stock==='A'&&!f.prohibited),r={};NS.forEach(n=>r[n]=Math.max(0,(+t[n]||0)-(n==='Ca'?b.Ca:n==='N'?b.N:0)));let C=act.map(f=>NS.map(n=>coef(f[n]||0,c.volA,c.denom))),x=act.map(()=>0),w=NS.map(n=>1/Math.pow(Math.max(1,+t[n]||0),2));
for(let pass=0;pass<4000;pass++){let delta=0;for(let j=0;j<act.length;j++){let a=0,z=0;for(let i=0;i<NS.length;i++){let pred=0;for(let k=0;k<act.length;k++)pred+=C[k][i]*x[k];let without=pred-C[j][i]*x[j];a+=w[i]*C[j][i]*(r[NS[i]]-without);z+=w[i]*C[j][i]*C[j][i]}let nx=z>0?Math.max(0,a/z):0;delta=Math.max(delta,Math.abs(nx-x[j]));x[j]=nx}if(delta<1e-7)break}
let out={},total={},err={};NS.forEach((n,i)=>{out[n]=C.reduce((s,row,j)=>s+row[i]*x[j],0);total[n]=out[n]+(n==='Ca'?b.Ca:n==='N'?b.N:0);err[n]=total[n]-(+t[n]||0)});let warnings=[];if(b.N>(+t.N||0)+.01)warnings.push('نیتروژن واردشده از نیترات کلسیم ('+b.N.toFixed(2)+' mg/L) از هدف N بیشتر است.');if(act.some(f=>f.id==='magnesiumSulfate')&&(+t.Mg||0)*13/9.86>(+t.S||0)+.01)warnings.push('تأمین Mg از سولفات منیزیم، گوگرد بیشتری از هدف S وارد می‌کند.');warnings.push('حل‌گر تقریبی است؛ کمینه‌سازی جرم مرحلهٔ جداگانه ندارد.');warnings.push('انحلال‌پذیری و سازگاری شیمیایی تأیید نشده‌اند؛ پیش از مصرف کنترل فنی لازم است.');return{b,active:act.map((f,i)=>({name:f.name,mass:x[i]})).filter(q=>q.mass>1e-5),total,err,warnings,aConc:x.reduce((s,v)=>s+v,0)/c.volA,bConc:b.mass/c.volB}}
g.NutritionSolver={NS,TARGETS,FERTS,coef,calcB,solve}})(typeof window!=='undefined'?window:globalThis);
if(typeof module!=='undefined'&&module.exports) module.exports=globalThis.NutritionSolver;
