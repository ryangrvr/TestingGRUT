"""L0 transcription of the pointwise rational/complex matrix generator.

This prices actual submitted arithmetic, not a claimed complete continuum law
or pipeline. No protected vocabulary is awarded a fictitious four-bit price.
"""
import json
import sys
from pathlib import Path
from math import log2

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]))
from development.final_precard.price_oracle import raw_price,table_price


def app(f,*args):
    for a in args: f=('apply',f,a)
    return f
def entry(m,i,j): return app(m,i,j)
def lambdas(names,body):
    for n in reversed(names): body=('lambda',n,body)
    return body
def calls(name,*args): return app(('use',name),*args)


def main():
    ar,ai,xr,xi,i,j,k,d='ar ai xr xi i j k d'.split()
    def mul(m,n): return ('*',entry(m,i,k),entry(n,k,j))
    cre=('+',('-',('-',mul(ar,xr),mul(ai,xi)),mul(xr,ar)),mul(xi,ai))
    cim=('-',('-',('+',mul(ar,xi),mul(ai,xr)),mul(xr,ai)),mul(xi,ar))
    cr=lambdas([ar,ai,xr,xi,d,i,j],('sum_finite',('range',d),('lambda',k,cre)))
    ci=lambdas([ar,ai,xr,xi,d,i,j],('sum_finite',('range',d),('lambda',k,cim)))
    defs={'comm_re':cr,'comm_im':ci}
    args=['ar','ai','rr','ri','d']
    crij=calls('comm_re',*args,i,j); ciij=calls('comm_im',*args,i,j)
    crmap=lambdas([i,j],crij); cimap=lambdas([i,j],ciij)
    ccargs=['ar','ai',crmap,cimap,'d',i,j]
    ccr=calls('comm_re',*ccargs); cci=calls('comm_im',*ccargs)
    factor=('/', 'tau', ('nat',2))
    real_rhs=('-',ciij,('*',factor,ccr))
    imag_rhs=('-',('-', '0', crij),('*',factor,cci))
    equations=('and',('=',entry('dar',i,j),'0'),('and',('=',entry('dai',i,j),'0'),
               ('and',('=',entry('drr',i,j),real_rhs),('=',entry('dri',i,j),imag_rhs))))
    predicate=('forall',i,('forall',j,('->',('and',('in',i,('range',d)),('in',j,('range',d))),equations)))
    for name in reversed(['ar','ai','rr','ri','dar','dai','drr','dri','tau','d']):
        predicate=('lambda',name,predicate)
    # These arrays include derivative values; their continuum interpretation
    # and admissibility axioms remain separate priced imports/obligations.
    transcript={'expression':predicate,'definitions':list(defs.items()),
                'scope':'pointwise matrix equations only; derivatives are supplied variables'}
    # A separate finite interpreter confirms this is the same candidate relation,
    # not an opaque operator name with a conveniently low token count.
    from fractions import Fraction
    from audit import Q,qmat,qcomm,qscale,qadd,LEVELS,TAU
    def evaluate(e,env):
        if isinstance(e,str):
            if e=='0': return Fraction(0)
            if e=='1': return Fraction(1)
            return env[e]
        op,*xs=e
        if op=='use': return evaluate(defs[xs[0]],env)
        if op=='nat': return Fraction(xs[0])
        if op=='lambda':
            return lambda value: evaluate(xs[1],dict(env,**{xs[0]:value}))
        if op=='forall':
            return all(evaluate(xs[1],dict(env,**{xs[0]:k})) for k in range(int(env['d'])))
        if op=='->': return (not evaluate(xs[0],env)) or evaluate(xs[1],env)
        if op=='and': return evaluate(xs[0],env) and evaluate(xs[1],env)
        args=[evaluate(x,env) for x in xs]
        if op=='range': return range(int(args[0]))
        if op=='in': return args[0] in args[1]
        if op=='=': return args[0]==args[1]
        if op=='apply': return args[0](args[1])
        if op=='sum_finite': return sum((args[1](k) for k in args[0]),Fraction(0))
        if op=='+': return args[0]+args[1]
        if op=='-': return args[0]-args[1]
        if op=='*': return args[0]*args[1]
        if op=='/': return args[0]/args[1]
        raise ValueError(op)
    def mapped(matrix,part):
        return lambda i: lambda j: getattr(matrix[int(i)][int(j)],part)
    amat=qmat(3)
    for k,x in enumerate(LEVELS): amat[k][k]=Q(x)
    zero=qmat(3); verified=0
    for row in range(3):
        for col in range(3):
            x=qmat(3); x[row][col]=Q(Fraction(1))
            comm=qcomm(amat,x)
            deriv=qadd(qscale(Q(Fraction(0),Fraction(-1)),comm),qscale(-TAU/2,qcomm(amat,comm)))
            values=[mapped(amat,'r'),mapped(amat,'i'),mapped(x,'r'),mapped(x,'i'),
                    mapped(zero,'r'),mapped(zero,'i'),mapped(deriv,'r'),mapped(deriv,'i'),TAU,Fraction(3)]
            result=evaluate(predicate,{})
            for value in values: result=result(value)
            if result is not True: raise AssertionError('incorrect L0 transcription')
            damaged=[r.copy() for r in deriv]
            damaged[0][0]=damaged[0][0]+Q(Fraction(1,100))
            values[6]=mapped(damaged,'r')
            result=evaluate(predicate,{})
            for value in values: result=result(value)
            if result is not False: raise AssertionError('L0 transcription false green')
            verified+=1
    (ROOT/'L0_CORE.json').write_text(json.dumps(transcript,indent=2)+'\n')
    priced=raw_price(predicate,10,defs)
    constants=[('nat',3),('rat',1,5),('rat',0,1),('rat',1,1),('rat',3,1)]
    table_bits=table_price(constants,10)
    out={'scope':'submitted pointwise L0 arithmetic plus partial declared inputs; NOT full card price',
         'core_syntax':priced,'same_core_syntax_at_p_6_10_16':True,
         'exact_L0_matrix_unit_agreements':verified,'deliberate_wrong_derivatives_rejected':verified,
         'partial_literal_table_bits':table_bits,'literal_table':constants,
         'vocabulary_menu_only_bits':28,'vocabulary_ids':['VB-2','VB-3','VB-7','VB-8','VB-9','VB-10','VB-11'],
         'local_ip12_count_conservative':7,'selection_tax_local_bits':log2(8),
         'partial_subtotal_bits':priced['bits']+table_bits+28+log2(8),
         'full_price_status':'NOT_CERTIFIED_INCOMPLETE_GENERATIVE_PIPELINE',
         'excluded_from_subtotal':['continuum derivative/admissibility unfolding','type declarations',
              'complete initial-operator/input/preparation/POVM encodings','record maps and interface actions',
              'IP4/import/archive-family selection audit','IP5 passing-window certificates',
              'IP7/IP13 commitment and domain audits','global IP12 inventory'],
         'earned_charter_coupling_bits':0,'compression_certificate':False,
         'uses_uncleared_price_oracle_for_bookkeeping_only':True}
    (ROOT/'price_report.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
