# -*- coding: utf-8 -*-
"""
methodology_upgrade.py — reproducibility script for the methodology-upgrade analyses
of the HPV post-pandemic time-series study (Wu et al.).

Reproduces, directly from the raw genotyping export (available from the corresponding
author on reasonable request; same input as qc_verification_18185.py):
  * External age-standardisation to the WHO World Standard Population (Table S9 / Fig S6a)
  * Genotype-specific trend by two methods: primary log-linear AAPC (from the canonical
    bootstrap table) vs an independent harmonic seasonal Poisson/negative-binomial GLM
    (Table S10 / Fig S6b)
  * Supplementary Figure S6 (2-panel, house style)

Self-check: reproduces the canonical yearly any-STABLE15 prevalence before running.
Deterministic: numpy default_rng seeded at 20260520.

Usage: set HPV_INPUT (raw export) and HPV_CANON (canonical_aapc_table.json), then run.
"""
import os, numpy as np, pandas as pd, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import statsmodels.api as sm, statsmodels.formula.api as smf

INPUT  = os.environ.get("HPV_INPUT", "2021_2025HPV _12_31.xlsx")
CANON_JSON = os.environ.get("HPV_CANON", "data/canonical_aapc_table.json")
OUTDIR = "outputs_methodology_upgrade"; os.makedirs(OUTDIR, exist_ok=True)
rng = np.random.default_rng(20260520)
CANON_PCT = [9.68, 13.34, 10.87, 15.08, 15.79]
_CANON = json.load(open(CANON_JSON, encoding="utf-8"))["aapc_pr"]
_CKEY = {'HPV16_p':'HPV16','HPV52_p':'HPV52','HPV58_p':'HPV58','HPV51_p':'HPV51'}
def canon_aapc(col):
    d=_CANON[_CKEY[col]]; return d['aapc'], d['aapc_lo'], d['aapc_hi'], d['aapc_p']

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
 'font.size':9,'axes.titlesize':9,'axes.titleweight':'bold','axes.titlelocation':'left','axes.titlepad':4,
 'axes.labelsize':8.5,'axes.linewidth':0.7,'axes.spines.top':False,'axes.spines.right':False,
 'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':7.5,'legend.frameon':False,
 'savefig.dpi':600,'savefig.bbox':'tight','pdf.fonttype':42,'lines.linewidth':1.4,'lines.markersize':4})
NAT={'aggregate':'#1f1f1f','pre':'#4575b4','post':'#d73027','increase':'#d73027','null':'#bdbdbd'}
def panel_label(ax,L,x=-0.075,y=1.06):
    ax.text(x,y,L,transform=ax.transAxes,fontsize=11,fontweight='bold',va='top',ha='left')

STABLE15=['HPV16','HPV18','HPV31','HPV33','HPV35','HPV39','HPV45','HPV51','HPV52','HPV53','HPV56','HPV58','HPV59','HPV66','HPV68']
df=pd.read_excel(INPUT); df.columns=[c.strip() for c in df.columns]
g=df['Gender'].astype(str).str.strip().str.lower()
df=df[g.eq('female')]; df=df[df['Age'].notna()]; df=df[df['Age']>0]; df=df[df['Age']>=21].copy()
ispos=lambda v:'positive' in str(v).strip().lower()
for c in STABLE15: df[c+'_p']=df[c].map(ispos)
df['Year']=pd.to_datetime(df['Detection time']).dt.year
df['ym']=pd.to_datetime(df['Detection time']).dt.to_period('M').astype(str)
df['anyS15']=df[[c+'_p' for c in STABLE15]].any(axis=1)
assert len(df)==18185
for i,y in enumerate(range(2021,2026)):
    assert abs(100*df[df.Year==y]['anyS15'].mean()-CANON_PCT[i])<0.05
print("self-check PASS: N=18185, canonical yearly prevalence reproduced")

# ===== U1: WHO / internal / crude age-standardisation (Table S9 / Fig S6a) =====
WHO={'20-24':8220,'25-29':7930,'30-34':7610,'35-39':7150,'40-44':6590,'45-49':6040,
     '50-54':5370,'55-59':4550,'60-64':3720,'65-69':2960,'70+':5275}
def abin(a):
    if a<25: return '20-24'
    if a>=70: return '70+'
    return "%d-%d"%(int(a//5*5),int(a//5*5)+4)
df['bin']=df['Age'].map(abin); bins=list(WHO)
w2021={b:int(df[(df.Year==2021)&(df.bin==b)].shape[0]) for b in bins}
def wilson(k,n):
    z=1.96; ph=k/n; d=1+z*z/n; c=(ph+z*z/(2*n))/d; h=z*np.sqrt(ph*(1-ph)/n+z*z/(4*n*n))/d
    return 100*(c-h),100*(c+h)
def std(year,wt):
    s=df[df.Year==year]; W=sum(wt[b] for b in bins if len(s[s.bin==b])>0); a=v=0.0
    for b in bins:
        sb=s[s.bin==b]; n=len(sb)
        if n==0: continue
        p=sb['anyS15'].mean(); a+=wt[b]*p; v+=(wt[b]/W)**2*p*(1-p)/n
    a/=W; se=np.sqrt(v); return 100*a,(100*(a-1.96*se),100*(a+1.96*se))
yrs=list(range(2021,2026)); crude=[];cci=[];who=[];wci=[];intn=[];ici=[]
for y in yrs:
    s=df[df.Year==y]; k=int(s['anyS15'].sum()); n=len(s)
    crude.append(100*k/n); cci.append(wilson(k,n))
    r,c=std(y,WHO); who.append(r); wci.append(c)
    r2,c2=std(y,{b:float(w2021[b]) for b in bins}); intn.append(r2); ici.append(c2)
aapc=lambda s:(np.exp(np.polyfit(np.arange(5),np.log(np.array(s)/100),1)[0])-1)*100
strata={y:{b:(int(df[(df.Year==y)&(df.bin==b)]['anyS15'].sum()),len(df[(df.Year==y)&(df.bin==b)])) for b in bins} for y in yrs}
def boot_std(wt):
    out=[]
    for _ in range(2000):
        ser=[]
        for y in yrs:
            st=strata[y]; W=sum(wt[b] for b in bins if st[b][1]>0); a=0.0
            for b in bins:
                if st[b][1]==0: continue
                a+=wt[b]*(rng.binomial(st[b][1],st[b][0]/st[b][1])/st[b][1])
            ser.append(100*a/W)
        out.append(aapc(ser))
    return np.percentile(out,[2.5,97.5])
btc=[]
for _ in range(2000):
    ser=[]
    for y in yrs:
        st=strata[y]; nt=sum(st[b][1] for b in bins); kt=sum(rng.binomial(st[b][1],st[b][0]/st[b][1]) for b in bins)
        ser.append(100*kt/nt)
    btc.append(aapc(ser))
ci_w=boot_std(WHO); ci_i=boot_std({b:float(w2021[b]) for b in bins})
ca=_CANON['Any STABLE15']  # crude AAPC from canonical

# ===== U2: harmonic GLM (4 genotypes; aggregate intentionally excluded) =====
mon=df.groupby('ym').size().rename('N').to_frame()
mc=lambda col: df.groupby('ym')[col].sum().reindex(mon.index).fillna(0).astype(int)
mon['t']=np.arange(len(mon))
for k,mm in [(1,1),(2,2)]: mon['s%d'%k]=np.sin(2*np.pi*mm*mon.t/12); mon['c%d'%k]=np.cos(2*np.pi*mm*mon.t/12)
mon['logN']=np.log(mon.N)
def glm(col):
    d=mon.copy(); d['y']=mc(col).values
    M=smf.glm("y~t+s1+c1+s2+c2",d,family=sm.families.Poisson(),offset=d.logN).fit()
    disp=M.pearson_chi2/M.df_resid; fam='Poisson'
    if disp>1.5:
        M=smf.glm("y~t+s1+c1+s2+c2",d,family=sm.families.NegativeBinomial(alpha=1.0),offset=d.logN).fit(); fam='NB'
    b=M.params['t']; ci=M.conf_int().loc['t']; f=lambda x:(np.exp(x*12)-1)*100
    return f(b),f(ci[0]),f(ci[1]),M.pvalues['t'],fam,disp
gtl=[('HPV16_p','HPV-16'),('HPV52_p','HPV-52'),('HPV58_p','HPV-58'),('HPV51_p','HPV-51')]
u2=[]
for col,lab in gtl:
    bp,bl,bh,bpv=canon_aapc(col); gp,gl,gh,gpv,fam,disp=glm(col)
    u2.append((lab,bp,bl,bh,bpv,gp,gl,gh,gpv,fam,disp))

# ===== FIGURE S6 (2-panel, house style) =====
fig=plt.figure(figsize=(174/25.4,5.6)); gs=GridSpec(2,1,height_ratios=[1.0,0.85],hspace=0.5,figure=fig)
axa=fig.add_subplot(gs[0]); axb=fig.add_subplot(gs[1])
for ser,ci,col,lab,a in [(crude,cci,NAT['aggregate'],'Crude',aapc(crude)),
                         (intn,ici,NAT['pre'],'Internal-2021 standard',aapc(intn)),
                         (who,wci,NAT['post'],'WHO World standard',aapc(who))]:
    axa.fill_between(yrs,[c[0] for c in ci],[c[1] for c in ci],color=col,alpha=0.15,lw=0)
    axa.plot(yrs,ser,'-o',color=col,ms=4,label="%s (AAPC %+.1f%%/yr)"%(lab,a))
axa.set_xticks(yrs); axa.set_xlabel("Year"); axa.set_ylabel("Any HR-HPV prevalence (%)")
axa.set_title("External age-standardisation"); axa.legend(loc='upper left'); panel_label(axa,'a')
rows=sorted(u2,key=lambda r:r[1]); yv=np.arange(len(rows))
for i,(lab,bp,bl,bh,bpv,gp,gl,gh,gpv,fam,disp) in enumerate(rows):
    axb.errorbar(bp,i+0.16,xerr=[[bp-bl],[bh-bp]],fmt='o',color=NAT['increase'] if (bl>0 or bh<0) else NAT['null'],ms=4.5,lw=1.2,capsize=2)
    axb.errorbar(gp,i-0.16,xerr=[[gp-gl],[gh-gp]],fmt='s',color=NAT['increase'] if (gl>0 or gh<0) else NAT['null'],ms=4.5,lw=1.2,capsize=2)
axb.axvline(0,color='#888',lw=0.8,ls='--'); axb.set_yticks(yv); axb.set_yticklabels([r[0] for r in rows])
axb.set_ylim(-0.6,len(rows)-0.4); axb.set_xlabel("Annualised trend (%/yr, 95% CI)"); axb.set_title("Genotype-specific trend by two methods")
axb.legend(handles=[Line2D([0],[0],marker='o',color=NAT['aggregate'],lw=0,ms=6,label='Log-linear AAPC'),
                    Line2D([0],[0],marker='s',color=NAT['aggregate'],lw=0,ms=6,label='Harmonic GLM'),
                    Line2D([0],[0],marker='o',color=NAT['increase'],lw=0,ms=6,label='95% CI excludes 0'),
                    Line2D([0],[0],marker='o',color=NAT['null'],lw=0,ms=6,label='includes 0')],loc='lower right',fontsize=6.6)
panel_label(axb,'b')
for ext in ['pdf','png']:
    fig.savefig(os.path.join(OUTDIR,"FigS6_methodology_upgrade."+ext),dpi=600 if ext=='png' else None,bbox_inches='tight',facecolor='white')
plt.close(fig)

# ===== Tables S9 + S10 to xlsx =====
fmt=lambda p,l,h,d=2:"%.*f (%.*f to %.*f)"%(d,p,d,l,d,h)
s9=pd.DataFrame([["Crude"]+[fmt(crude[i],cci[i][0],cci[i][1]) for i in range(5)]+[fmt(ca['aapc'],ca['aapc_lo'],ca['aapc_hi'],1)],
                 ["Internal 2021 standard"]+[fmt(intn[i],ici[i][0],ici[i][1]) for i in range(5)]+[fmt(aapc(intn),ci_i[0],ci_i[1],1)],
                 ["WHO World standard"]+[fmt(who[i],wci[i][0],wci[i][1]) for i in range(5)]+[fmt(aapc(who),ci_w[0],ci_w[1],1)]],
                columns=["Standardisation"]+[str(y) for y in yrs]+["AAPC %/yr (95% CI)"])
s10=pd.DataFrame([[lab,fmt(bp,bl,bh,1),("<0.001" if bpv<0.001 else "%.3f"%bpv),fmt(gp,gl,gh,1),("%.3f"%gpv),"%s (%.2f)"%(fam,disp)] for lab,bp,bl,bh,bpv,gp,gl,gh,gpv,fam,disp in u2],
                columns=["Genotype","Log-linear AAPC %/yr (95% CI)","p","Harmonic GLM %/yr (95% CI)","p","GLM family (dispersion)"])
with pd.ExcelWriter(os.path.join(OUTDIR,"TablesS9_S10.xlsx")) as xw:
    s9.to_excel(xw,"S9",index=False); s10.to_excel(xw,"S10",index=False)
print("Fig S6 (2-panel) + Tables S9/S10 written to",OUTDIR)
print("S10:",[(r[0],round(r[1],1),round(r[5],1)) for r in u2])
