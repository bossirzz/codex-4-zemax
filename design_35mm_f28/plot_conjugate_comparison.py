"""Publication-quality comparison using native curves and exact native operands."""
import argparse,json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from joint_design_model import R

def run(tag):
 rows=json.loads((R/(tag+'_native_acceptance.json')).read_text(encoding='utf8'))
 font=FontProperties(fname=r'C:\Windows\Fonts\msyh.ttc')
 colors=['#202124','#4285f4','#15966e','#d23d35']
 fig,axes=plt.subplots(2,2,figsize=(14,10),gridspec_kw={'width_ratios':[1.1,1]})
 for j,(wd,suffix,fc,fe) in enumerate([(500,'_review_mtf_1024',225,180),(200,'_near_mtf_1024',150,120)]):
  df=pd.read_csv(R/(tag+suffix+'.csv'),header=[0,1],index_col=0);frequency=np.asarray(df.index,float)
  row=next(x for x in rows if x['wd']==wd);ax,detail=axes[j]
  for field,h in enumerate([0.,3.,6.3,9.08]):
   for k,style in enumerate(['-','--']):
    ax.plot(frequency,df.iloc[:,field*2+k],style,c=colors[field],lw=1.5,label=f'{h:g} mm '+['T','S'][k])
  for col,color,label in [(0,colors[0],'中心'),(6,colors[3],'边缘 T'),(7,'#e7891d','边缘 S')]:
   detail.plot(frequency,df.iloc[:,col],c=color,lw=2,label=label)
  for freq,value,color in [(fc,row['mtf_sampling6'][0],colors[0]),(fe,row['mtf_sampling6'][-2],colors[3]),(fe,row['mtf_sampling6'][-1],'#e7891d')]:
   detail.scatter([freq],[value],s=40,c=color,zorder=5)
  detail.scatter([fc,fe],[.30,.20],s=65,marker='x',c='#777',label='验收门槛')
  detail.text(.04,.96,f'中心 {fc} lp/mm：{row["mtf_sampling6"][0]:.4f}\n边缘 {fe} lp/mm：T {row["mtf_sampling6"][-2]:.4f} / S {row["mtf_sampling6"][-1]:.4f}',
   transform=detail.transAxes,va='top',fontproperties=font,fontsize=11,bbox={'facecolor':'white','edgecolor':'#ddd','alpha':.9})
  ax.set_title(f'WD {wd} mm · 名义 F/{row["nominal_fno"]:.3f} · 全部视场',fontproperties=font,fontsize=13)
  detail.set_title('中心与边缘：精确频率原生结果',fontproperties=font,fontsize=13)
  ax.set_xlim(0,300);ax.set_ylim(0,1);detail.set_xlim(fe-55,fc+35);detail.set_ylim(0,.55)
  for a in [ax,detail]:
   a.grid(alpha=.22);a.set_xlabel('空间频率（lp/mm）',fontproperties=font);a.set_ylabel('光学 MTF',fontproperties=font)
   a.legend(prop=font,loc='upper right',fontsize=9)
 fig.suptitle('Zemax 原生 FFT MTF · 1024×1024 · F/d/C 等权 · 裸像面',fontproperties=font,fontsize=16)
 fig.text(.5,.02,'曲线保留所有视场与 T/S；标注来自精确频率操作数。中间视场的下降必须与中心、边缘一起审查。',ha='center',fontproperties=font,fontsize=10)
 fig.tight_layout(rect=[0,.045,1,.95]);target=R/(tag+'_best_near_comparison.png')
 fig.savefig(target,dpi=230);fig.savefig(R/(tag+'_best_near_comparison.pdf'))
 print(str(target),flush=True)

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('tag');q=a.parse_args();run(q.tag)
