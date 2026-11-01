import sys
sys.path.insert(0,'tools')
import compare as C, subprocess
from pathlib import Path
p=sorted(Path('corpus/'+sys.argv[1]).rglob(sys.argv[2]+'.gv'))[0]
b=C.binary_path()
out=subprocess.run([b],input=C.script(C.parse(p)),capture_output=True,text=True).stdout.splitlines()
m=C.read_fdot(out); pa=C.parse(p); t=C.read_dot(p,'.table')
names=list(pa[1]); ids={n:i+1 for i,n in enumerate(names)}
dx=m['nodes'][1][2]-t['nodes'][names[0]][0]; dy=m['nodes'][1][3]-t['nodes'][names[0]][1]
print(m['bb'],t['bb'],dx,dy,len(names),len(pa[2]))
for n in names: print(n,m['nodes'][ids[n]][:2],round(m['nodes'][ids[n]][2]-dx,1),round(m['nodes'][ids[n]][3]-dy,1),t['nodes'][n][:2],pa[1][n].get('label','')[:14])
w={}
for a,b_,q in t['curves']: w[(a,b_)]=q
for k,(a,b_,at) in enumerate(pa[2],1):
    h=m['curves'].get(k,[]);q=w.get((a,b_),[])
    print(k,a,b_,at.get('label',''))
    print('  mine',[[(round(u[0]-dx,1),round(u[1]-dy,1)) for u in x] for x in h])
    print('  dot ',[[(round(u[0],1),round(u[1],1)) for u in x] for x in q])
