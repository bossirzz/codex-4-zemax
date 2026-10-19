"""Render original Zemax FFT MTF data as a readable, exportable plot."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
ROOT=Path(__file__).resolve().parent
font=FontProperties(fname=r'C:\Windows\Fonts\msyh.ttc')

def run(tag,suffix='_native_mtf'):
 source=ROOT/(tag+suffix+'.csv')
 df=pd.read_csv(source,header=[0,1],index_col=0);freq=np.asarray(df.index,float)
 fig,(ax,detail)=plt.subplots(1,2,figsize=(13.5,5.5),gridspec_kw=dict(width_ratios=[1.15,1]))
 colors=['#202124','#4285f4','#15966e','#d23d35'];heights=[0,3,6.3,9.08]
 for j,h in enumerate(heights):
  for k,style in enumerate(['-','--']):
   ax.plot(freq,df.iloc[:,j*2+k],style,color=colors[j],lw=1.6,label=f'{h:g} mm / '+['T','S'][k])
 for col,color,label in [(0,colors[0],'中心 T/S'),(6,colors[3],'边缘 T'),(7,'#e7891d','边缘 S')]:
  detail.plot(freq,df.iloc[:,col],color=color,lw=2,label=label)
 points=[]
 for col,x,target,color in [(0,225,.30,colors[0]),(6,180,.20,colors[3]),(7,180,.20,'#e7891d')]:
  y=float(np.interp(x,freq,df.iloc[:,col]));detail.scatter([x],[y],color=color,s=35,zorder=4)
  points.append(dict(column=col,frequency=x,mtf_interpolated=y,target=target))
 detail.scatter([225,180],[.30,.20],marker='x',s=60,c='#777777',label='暂定验收门槛')
 detail.annotate(f'225 lp/mm: {points[0]["mtf_interpolated"]:.3f}',(225,points[0]['mtf_interpolated']),xytext=(184,.48),arrowprops=dict(arrowstyle='-',color=colors[0]),fontproperties=font,fontsize=10)
 detail.annotate(f'180 lp/mm: T {points[1]["mtf_interpolated"]:.3f} / S {points[2]["mtf_interpolated"]:.3f}',(180,points[1]['mtf_interpolated']),xytext=(130,.075),arrowprops=dict(arrowstyle='-',color=colors[3]),fontproperties=font,fontsize=9)
 for a in [ax,detail]:
  a.set_xlabel('空间频率（lp/mm）',fontproperties=font);a.set_ylabel('光学 MTF',fontproperties=font)
  a.grid(alpha=.22);a.legend(prop=font,fontsize=9,loc='upper right');a.set_ylim(0,1)
 ax.set_xlim(0,300);detail.set_xlim(125,260);detail.set_ylim(0,.55)
 ax.set_title('完整曲线：所有视场',fontproperties=font,fontsize=12)
 detail.set_title('中心与边缘：验收频率放大',fontproperties=font,fontsize=12)
 fig.suptitle('Zemax 原生 FFT MTF · WD 500 mm · F/2.8 · 真实光阑瞄准',fontproperties=font,fontsize=15)
 fig.text(.5,.015,'F/d/C 等权可见光；固定物理光阑；T=子午、S=弧矢。数字由原始 CSV 插值，最终验收另用精确频率计算。',ha='center',fontproperties=font,fontsize=9,color='#555555')
 fig.tight_layout(rect=[0,.04,1,.92])
 image=ROOT/(tag+suffix+'.png');fig.savefig(image,dpi=230);fig.savefig(ROOT/(tag+suffix+'.pdf'))
 (ROOT/(tag+suffix+'_plot_metadata.json')).write_text(json.dumps(dict(source=str(source),points=points),indent=2),encoding='utf8')
 print(str(image));print(json.dumps(points))

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('tag');a.add_argument('--suffix',default='_native_mtf');q=a.parse_args();run(q.tag,q.suffix)
