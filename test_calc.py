"""Tests for my_calculator_fixed.py  (run:  python test_calc.py)
Keep this file in the same folder as my_calculator_fixed.py. No extra libraries needed."""
import sys, os, types, math, cmath
try:
    import tkinter  # noqa: F401  (real tkinter is fine; we never open a window)
except ImportError:  # machines without Tkinter can still run the maths tests
    tk = types.ModuleType("tkinter"); tk.messagebox = types.ModuleType("m"); tk.simpledialog = types.ModuleType("s")
    sys.modules["tkinter"] = tk; sys.modules["tkinter.messagebox"] = tk.messagebox; sys.modules["tkinter.simpledialog"] = tk.simpledialog
    tk.Tk = object; tk.Frame = object; tk.Label = object
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import my_calculator_fixed as mc
ok=bad=0
def chk(name, got, exp, tol=1e-7):
    global ok,bad
    try:
        if isinstance(exp,str): good = (got==exp)
        else: good = abs(complex(got)-complex(exp))<=tol*max(1,abs(complex(exp)))
    except Exception: good=False
    if good: ok+=1
    else: bad+=1; print("FAIL",name,"got",got,"expected",exp)
def E(mode="COMP",angle="DEG",base="DEC"):
    e=mc.Engine(); e.mode=mode; e.angle=angle; e.base=base; return e
def val(e,t):
    try: return e.evaluate(t)
    except Exception as ex: return "ERR:"+type(ex).__name__+":"+str(ex)

e=E()
cases=[
("2+3×4",14),("(2+3)×4",20),("10÷4",2.5),("2^10",1024),("2^(−1)",0.5),
("√(144)",12),("∛(27)",3),("∛(−8)",-2),("ˣ√(3,27)",3),("58ˣ√(8)",8**(1/58)),("3ˣ√(27)",3),("4ˣ√(81)",3),("5!",120),("0!",1),("10nPr3",720),("10nCr3",120),
("sin(30)",0.5),("cos(60)",0.5),("tan(45)",1),("sin⁻¹(0.5)",30),("cos⁻¹(0.5)",60),("tan⁻¹(1)",45),
("sinh(1)",math.sinh(1)),("cosh(1)",math.cosh(1)),("tanh(1)",math.tanh(1)),
("sinh⁻¹(1)",math.asinh(1)),("cosh⁻¹(2)",math.acosh(2)),("tanh⁻¹(0.5)",math.atanh(.5)),
("log(1000)",3),("ln(ℯ)",1),("logb(2,8)",3),("10^(2)",100),("ℯ^(1)",math.e),
("π",math.pi),("2π",2*math.pi),("π×2²",4*math.pi),("4²",16),("3³",27),("4⁻¹",0.25),
("Abs(−5)",5),("Abs(3−7)",4),("∫(X^2,0,3)",9),("∫(sin(X),0,180)",None),("d/dx(X^3,2)",12),("Σ(X,1,10)",55),("Σ(X^2,1,5)",55),
("−5+2",-3),("2×−3",-6),("5−−3",8),("1÷3×3",1),("0.1+0.2",0.3),("1.5×10^3",1500),("6×10^(−3)",0.006),
("Rnd(1.23456789012)",1.23456789),("2(3+4)",14),("(1+1)(2+2)",8),("3π",3*math.pi),
("30°",None),("sin(30°)",0.5),
]
for t,exp in cases:
    v=val(e,t)
    if exp is None: print("info",t,"->",v); continue
    chk(t,v,exp)
# Pol/Rec
v=val(e,"Pol(3,4)"); chk("Pol r",v,5); chk("Pol theta",e.vars["Y"],math.degrees(math.atan2(4,3)))
val(e,"Rec(2,60)"); chk("Rec x",e.vars["X"],1); chk("Rec y",e.vars["Y"],math.sqrt(3))
# angle modes
e2=E(angle="RAD"); chk("sin rad",val(e2,"sin(π÷2)"),1); chk("asin rad",val(e2,"sin⁻¹(1)"),math.pi/2)
e3=E(angle="GRA"); chk("sin gra",val(e3,"sin(100)"),1); chk("cos gra",val(e3,"cos(200)"),-1)
# errors expected
for t in ["1÷0","√(−1)","ln(0)","(−8)^(0.5)","(","5!!","−5!","70!","2.5!","log(−1)","sin⁻¹(2)","1+","×3"]:
    print("errcase",t,"->",val(e,t))
# format
for v in [1/3,2/3,1e10,1.5e-5,123456789012,0.1+0.2,1e-12,-0.0,100,2**0.5]:
    print("fmt",v,e.fmt(v))
e.disp=("fix",3); print("fix3",e.fmt(2/3)); e.disp=("sci",4); print("sci4",e.fmt(12345.678)); e.disp=("norm",2); print("norm2",e.fmt(0.001))
e.disp=("norm",1)
# exact
for v in [1/3,0.75,math.pi/2,2*math.pi,math.sqrt(2),3*math.sqrt(2),0.5,1.0,7/3]:
    print("exact",v,e.exact(v))
# complex
ec=E("CMPLX")
chk("cmul",val(ec,"(3+4i)×(1−2i)"),11-2j); chk("csqrt",val(ec,"√(−4)"),2j); chk("cabs",val(ec,"Abs(3+4i)"),5)
chk("cdiv",val(ec,"(1+i)÷(1−i)"),1j); chk("i^2",val(ec,"i^2"),-1); chk("conj",val(ec,"Conjg(3+4i)"),3-4j)
chk("arg",val(ec,"Arg(1+i)"),45); chk("Re",val(ec,"Re(3+4i)"),3); chk("Im",val(ec,"Im(3+4i)"),4)
chk("ln(-1)",val(ec,"ln(−1)"),1j*math.pi); chk("sin(i) deg",val(ec,"sin(i)"),cmath.sin(1j*math.pi/180))
print("cfmt",ec.fmt(11-2j),ec.fmt(2j),ec.fmt(3+0j)); ec.polar=True; print("polar",ec.fmt(1+1j))
print("cplx in COMP ->",val(e,"i×i"))
# matrix
A=mc.Matrix([[2,1],[5,3]]); B=mc.Matrix([[1,2],[3,4]])
chk("det",A.det(),1); chk("inv",A.inv().m[0][0],3); chk("inv12",A.inv().m[0][1],-1)
chk("mul",(A*B).m[1][1],5*2+3*4); chk("add",(A+B).m[0][0],3); chk("sub",(A-B).m[0][0],1); chk("scalar",(2*A).m[0][1],2)
chk("pow2",(A**2).m[0][0],2*2+1*5); chk("trn",A.trn().m[0][1],5)
M3=mc.Matrix([[1,2,3],[0,1,4],[5,6,0]]); chk("det3",M3.det(),1)
I=M3*M3.inv(); chk("M*Minv",I.m[0][0],1,1e-9); chk("M*Minv off",I.m[0][1],0,1e-9)
for f,nm in [(lambda:A+mc.Matrix([[1,2,3]]),"dim add"),(lambda:mc.Matrix([[1,2]])*mc.Matrix([[1,2]]),"dim mul"),(lambda:mc.Matrix([[1,2],[2,4]]).inv(),"singular inv")]:
    try: f(); print("FAIL no error",nm); bad+=1
    except Exception as ex: ok+=1
em=E("MATRIX"); em.mats["A"]=A; em.mats["B"]=B
r=val(em,"MatA×MatB"); print("MatA×MatB",r.m if hasattr(r,'m') else r)
r=val(em,"MatA^(−1)"); print("MatA^-1",r.m if hasattr(r,'m') else r)
r=val(em,"MatA⁻¹"); print("MatA⁻¹",r.m if hasattr(r,'m') else r)
r=val(em,"det(MatA)"); print("det(MatA)",r)
r=val(em,"Trn(MatA)"); print("Trn",r.m if hasattr(r,'m') else r)
r=val(em,"Ident(3)"); print("Ident",r.m if hasattr(r,'m') else r)
r=val(em,"2MatA"); print("2MatA",r.m if hasattr(r,'m') else r)
r=val(em,"MatA+MatB"); print("A+B",r.m if hasattr(r,'m') else r)
# vector
ev=E("VECTOR"); ev.vcts["A"]=mc.Vec([1,2,3]); ev.vcts["B"]=mc.Vec([4,5,6])
r=val(ev,"VctA×VctB"); print("cross",r.v if hasattr(r,'v') else r)
print("dot",val(ev,"Dot(VctA,VctB)"),"angle",val(ev,"Angle(VctA,VctB)"),"abs",val(ev,"Abs(VctA)"),"unit",getattr(val(ev,"UnitV(VctA)"),'v',None))
print("vec+",getattr(val(ev,"VctA+VctB"),'v',None),"vec-",getattr(val(ev,"VctA−VctB"),'v',None),"scal",getattr(val(ev,"2VctA"),'v',None))
# eqn
chk("lin2",mc.solve_linear([[2,1],[1,3]],[5,10])[0],1); 
chk("lin3",mc.solve_linear([[1,1,1],[0,2,5],[2,5,-1]],[6,-4,27])[0],5)
print("quad",mc.quad_roots(1,-5,6),mc.quad_roots(1,2,5),mc.quad_roots(1,2,1),mc.quad_roots(1,0,-2))
print("cubic",mc.poly_roots([1,-6,11,-6]),mc.poly_roots([1,0,0,-8]),mc.poly_roots([2,-4,-22,24]))
try: mc.solve_linear([[1,2],[2,4]],[1,2]); print("FAIL singular linear"); bad+=1
except ValueError: ok+=1
# solve
e=E()
print("solve X^2=2:",e.solve("X^2=2",1),"solve sin X=0.5:",e.solve("sin(X)=0.5",10),"solve X^2+1=0:",end=" ")
try: print(e.solve("X^2+1=0",1))
except Exception as ex: print("ERR (expected)",type(ex).__name__)
# stat
items,pred=mc.stat_calc(0,[4,8,6,5,3,7],[]); d=dict(items)
chk("mean",d["x̄"],5.5); chk("pop sd",d["σx"],math.sqrt(sum((x-5.5)**2 for x in [4,8,6,5,3,7])/6)); chk("sample sd",d["sx"],math.sqrt(17.5/5))
print("1var",items)
for xs in ([1,2,3,4,5,6,7],[1,2,3,4,5,6,7,8],[3,1,4,1,5,9,2,6]):
    it,_=mc.stat_calc(0,xs,[]); dd=dict(it); import statistics as st
    print("quartile",xs,"Q1",dd["Q1"],"Med",dd["Med"],"Q3",dd["Q3"],"| statistics median",st.median(xs))
xs=[1,2,3,4,5]; ys=[2.1,3.9,6.2,7.8,10.1]
it,pr=mc.stat_calc(1,xs,ys); d=dict(it)
def linfit(xs, ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    sxx=sum((x-mx)**2 for x in xs); sxy=sum((x-mx)*(y-my) for x,y in zip(xs,ys)); syy=sum((y-my)**2 for y in ys)
    b=sxy/sxx; return my-b*mx, b, sxy/math.sqrt(sxx*syy)
def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
def quadfit(xs, ys):  # Cramer's rule, independent of the code under test
    S=lambda p,q=0: sum(x**p*y**q for x,y in zip(xs,ys)); n=len(xs)
    M=[[n,S(1),S(2)],[S(1),S(2),S(3)],[S(2),S(3),S(4)]]; R=[S(0,1),S(1,1),S(2,1)]; D=det3(M)
    out=[]
    for k in range(3):
        Mk=[row[:] for row in M]
        for i in range(3): Mk[i][k]=R[i]
        out.append(det3(Mk)/D)
    return out
it,pr=mc.stat_calc(1,xs,ys); d=dict(it)
a,b,r=linfit(xs,ys); chk("A",d["A"],a); chk("B",d["B"],b); chk("r",d["r"],r)
chk("yhat",pr[0](6)[0],a+b*6); chk("xhat",pr[1](10)[0],(10-a)/b)
it,pr=mc.stat_calc(2,xs,ys); d=dict(it); a2,b2,c2=quadfit(xs,ys); chk("quadA",d["A"],a2); chk("quadB",d["B"],b2); chk("quadC",d["C"],c2)
yx=[2*x**2 for x in xs]
it,pr=mc.stat_calc(6,xs,yx); d=dict(it); chk("pow A",d["A"],2); chk("pow B",d["B"],2)
ye=[3*1.5**x for x in xs]
it,pr=mc.stat_calc(5,xs,ye); d=dict(it); chk("A·B^X A",d["A"],3); chk("A·B^X B",d["B"],1.5)
it,pr=mc.stat_calc(4,xs,[2*math.exp(0.3*x) for x in xs]); d=dict(it); chk("e^X A",d["A"],2); chk("e^X B",d["B"],0.3)
it,pr=mc.stat_calc(3,xs,[1+2*math.log(x) for x in xs]); d=dict(it); chk("ln A",d["A"],1); chk("ln B",d["B"],2)
it,pr=mc.stat_calc(7,xs,[1+2/x for x in xs]); d=dict(it); chk("1/x A",d["A"],1); chk("1/x B",d["B"],2)
# constants & conversion
cv=dict(mc.CONV)
chk("in>cm",cv["in ▸ cm"](1),2.54); chk("F>C",cv["°F ▸ °C"](212),100); chk("C>F",cv["°C ▸ °F"](100),212); chk("mile>km",cv["mile ▸ km"](1),1.609344)
chk("km/h>m/s",cv["km/h ▸ m/s"](36),10); chk("lb>kg",cv["lb ▸ kg"](1),0.45359237); chk("J>cal",cv["J ▸ cal"](4.18605),1); chk("kg>lb",cv["kg ▸ lb"](0.45359237),1)
print("const names unique:",len({c[0] for c in mc.CONSTS})==len(mc.CONSTS), [c[0] for c in mc.CONSTS if [x[0] for x in mc.CONSTS].count(c[0])>1])
print("RESULT ok",ok,"bad",bad)
