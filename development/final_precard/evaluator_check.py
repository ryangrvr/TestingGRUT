"""Compare hash-frozen independent arithmetic with supplied/proposed modules."""
from fractions import Fraction as F
import itertools
from . import price_oracle as P, hb_oracle as H


def attempt(f):
    try: return {'value':f()}
    except Exception as exc: return {'error':type(exc).__name__}


def prices(old,new):
    rows=[]
    def row(name,expected,before,after,kind='EXACT_VALUE'):
        before,after=attempt(before),attempt(after)
        def agrees(v):
            return ('error' in v) if kind=='REJECT' else v.get('value')==expected
        rows.append({'fixture':name,'kind':kind,'expected':expected,'before':before,
            'after':after,'before_agrees':agrees(before),'agrees':agrees(after)})
    expressions=[]
    for op,n in dict(P.ARITY,**P.EXPLICIT_SIGNATURES).items():
        args=['x','y','z'][:n]
        if op in P.BINDERS: args=['x',('=','x','y')]
        expressions.append((op,*args))
    expressions.extend(('=', 'x',('nat',n)) for n in (0,1,2,7,15,127,1023,65535))
    expressions.extend(('=', 'x',('rat',n,d)) for n,d in ((1,2),(-1,2),(7,3),(-127,64)))
    expressions.extend(('=', 'x',('real',e)) for e in (-64,-3,0,1,7))
    for i,e in enumerate(expressions):
        for p in (6,10,16):
            row(f'raw-{i}-p{p}',P.raw_price(e,p,signatures=P.EXPLICIT_SIGNATURES)['bits'],
                lambda e=e,p=p:old.L_stmt(e,p),
                lambda e=e,p=p:new.L_stmt(e,p,arities=P.EXPLICIT_SIGNATURES))
    row('alphabet-64',list(P.ALPHABET),lambda:old.ALPHABET,lambda:new.ALPHABET)
    for k in range(1,257):
        for offset in (-1,0,1):
            n=2**k+offset
            row(f'Elias-2^{k}{offset:+}',P.delta(n),lambda n=n:old.elias_delta_len(n),lambda n=n:new.elias_delta_len(n))
    invalid=[(),[],None,3,False,('unpriced',),('nat',0,'extra'),('def',),('end',),
        ('var','x'),('=',),('forall',('nat',2),'top'),('nat',-1),('nat',1.5),
        ('nat',True),('rat',1,0),('rat',1,-1),('rat',0.5,2),('real',1.5),
        ('real',1,'ignored'),('and','top','top','bot')]
    invalid.extend(P.RESERVED)
    for i,e in enumerate(invalid):
        row('invalid-AST-'+str(i),'ERROR',lambda e=e:old.L_stmt(e),lambda e=e:new.L_stmt(e), 'REJECT')
    for op,n in P.ARITY.items():
        for count in (n+1,n-1):
            if count<0: continue
            e=(op,*(['x']*count))
            row(f'arity-{op}-{count}','ERROR',lambda e=e:old.L_stmt(e),lambda e=e:new.L_stmt(e),'REJECT')
    for name,defs in (
        ('one',{'d':('=', 'x',('nat',3))}),
        ('chain',{'d':('use','e'),'e':('=', 'x',('nat',3))}),
    ):
        expr=('and',('use','d'),('use','d'))
        row('definitions-'+name,P.raw_price(expr,definitions=defs)['bits'],
            lambda:old.L_stmt(expr),lambda expr=expr,defs=defs:new.L_stmt(expr,definitions=defs))
    for defs in ({'d':('use','d')},{'d':('use','e'),'e':('use','d')},[('d','top'),('d','bot')]):
        row('bad-definitions-'+repr(defs),'ERROR',lambda:old.L_stmt('top'),
            lambda defs=defs:new.L_stmt('top',definitions=defs),'REJECT')
    atoms=[('=','x','q0'),('!=','x','q1'),('=','q0','q1'),'A','B']
    forms=[(op,('forall','x',a),('exists','x',b))
           for op in ('and','or','->','<->') for a,b in itertools.product(atoms,repeat=2)]
    forms += [('not',f) for f in forms[:10]]
    forms += [('ite','A',('forall','x',atoms[0]),('exists','x',atoms[1])),
              ('forall','x',('=','q0',('apply',('lambda','x','x'),'x'))),
              ('and',('=',('apply',('lambda','x','x'),'y'),'y'),
                     ('=',('apply',('lambda','x','x'),'z'),'z')),
              ('and',('=',('apply',('lambda','a','a'),'y'),'y'),
                     ('=',('apply',('lambda','b','b'),'z'),'z'))]
    semantics=0; failures=[]
    def semantic(module,expr):
        nonlocal semantics
        norm=module.normalize(expr)
        if P.free_variables(norm)!=P.free_variables(expr): return False
        for values in itertools.product((0,1),repeat=7):
            env=dict(zip(('x','q0','q1','A','B','y','z'),values)); semantics+=1
            if bool(P.evaluate(expr,env))!=bool(P.evaluate(norm,env)): return False
        return True
    for i,e in enumerate(forms):
        row('normal-semantics-'+str(i),True,lambda e=e:semantic(old,e),lambda e=e:semantic(new,e))
        for p in (6,10,16):
            row(f'normal-price-{i}-p{p}',P.normal_price(e,p)['bits'],
                lambda e=e,p=p:old.L_stmt(old.normalize(e),p),lambda e=e,p=p:new.statement_price(e,p))
    row('quantified-term-rejection','ERROR',lambda:old.normalize(('=',('forall','x','A'),'B')),
        lambda:new.normalize(('=',('forall','x','A'),'B')),'REJECT')
    for p in (-1,0,1.5,True):
        row('invalid-precision-'+str(p),'ERROR',lambda p=p:old.L_stmt(('real',0),p),lambda p=p:new.L_stmt(('real',0),p),'REJECT')
    for L,m in ((6,1024),(0,1),(30,3)):
        row('selection-'+str((L,m)),P.selection_price(L,m),lambda L=L,m=m:old.ip4_selection(L,m),lambda L=L,m=m:new.ip4_selection(L,m))
    for name,before,after in (
        ('negative-ledger',lambda:old.price([{'item':'x','class':'input','bits':-1}]),lambda:new.price([{'item':'x','class':'input','bits':-1}])),
        ('fractional-family',lambda:old.ip4_selection(10,1.5),lambda:new.ip4_selection(10,1.5)),
        ('fractional-draft-count',lambda:old.ip12_tax(0.5),lambda:new.ip12_tax(0.5)),
        ('negative-draft-count',lambda:old.ip12_tax(-1),lambda:new.ip12_tax(-1)),
        ('short-kappa-list',lambda:old.coupling_credit([1,0],1,1,1,kappa=([0],0,0,0)),lambda:new.coupling_credit([1,0],1,1,1,kappa=([0],0,0,0))),
        ('too-many-fibers',lambda:old.coupling_credit([1],1,1,2),lambda:new.coupling_credit([1],1,1,2)),
        ('empty-credit-list',lambda:old.coupling_credit([],1,1,0),lambda:new.coupling_credit([],1,1,0)),
    ): row(name,'ERROR',before,after,'REJECT')
    # No candidate is scored: abstract ledger rows exercise supplied-bit hooks.
    def synthetic_ledger(rule,encoding,bits,**extra):
        ledger=[{'item':f'empty-{i}','class':'input','bits':0,'encoding':None,
            'hostile_family_size':1,'rule_applied':f'IP-{i}',
            'not_applicable_reason':'ABSTRACT KIT FIXTURE','not_applicable_certificate_id':'fixture'} for i in range(1,15)]
        ledger[int(rule[3:])-1]=dict(item='tested',**{'class':'input'},bits=bits,encoding=encoding,
            hostile_family_size=extra.pop('hostile_family_size',1),rule_applied=rule,**extra)
        return ledger
    for rule,encoding,bits,extra in (
        ('IP-2',('nat',3),P.raw_price(('nat',3))['bits'],{}),
        ('IP-3',('pair',('nat',0),('nat',1)),P.raw_price(('pair',('nat',0),('nat',1)))['bits'],{}),
        ('IP-4','top',10,{'hostile_family_size':1024,'family_certificate_id':'fixture'}),
        ('IP-4','top',2,{'hostile_family_size':4,'relation_choice_only':True,'family_certificate_id':'fixture','forced_relation_certificate_id':'fixture'}),
        ('IP-2','top',P.input_price(6,24),{'is_supplied_prior':True,'prior_kl_bits':24,'prior_kl_certificate_id':'fixture'}),
        ('IP-3','emptyset',6+P.table_price([('nat',0),('nat',1)]),{'table_entries':[('nat',0),('nat',1)]}),
    ):
        ledger=synthetic_ledger(rule,encoding,bits,**extra)
        tag=rule+'-'+str(bits)+('-explicit-table' if 'table_entries' in extra else '')
        row('ledger-hook-'+tag,bits,lambda ledger=ledger:old.price(ledger),
            lambda ledger=ledger:new.price(ledger,checked_certificate_ids={'fixture'}))
        underpriced=synthetic_ledger(rule,encoding,bits-1,**extra)
        row('ledger-underprice-'+tag,'ERROR',lambda ledger=underpriced:old.price(ledger),
            lambda ledger=underpriced:new.price(ledger,checked_certificate_ids={'fixture'}),'REJECT')
    mismatches=[x for x in rows if not x['agrees']]
    return {'rows':rows,'comparison_count':len(rows),'semantic_evaluations':semantics,
            'mismatches':mismatches,'original_defects':sum(not x['before_agrees'] for x in rows)}


def baths(old,new):
    rows=[]
    def row(name,expected,before,after,reject=False):
        b,a=attempt(before),attempt(after)
        rows.append({'fixture':name,'expected':expected,'before':b,'after':a,
            'agrees':('error' in a) if reject else a.get('value')==expected,
            'before_agrees':('error' in b) if reject else b.get('value')==expected})
    identity=[[1,0],[0,1]]
    for G1,G2,M1,M2 in ((identity,[[2,1],[-1,3]],[0,0],[2,-3]),
                        ([[2,0],[0,3]],[[1,1],[0,2]],[1,-2],[0,3])):
        p=H.law([(0,0),(1,0),(0,1)],[1,2,3])
        def verified(subject):
            A,b,exact=subject.tlin_witness_map(M1,G1,M2,G2)
            source=H.push(p,G1,M1); target=H.push(p,G2,M2)
            return bool(exact and H.affine_certificate(source,target,A,b))
        row('affine-proof-'+str(len(rows)),True,lambda:verified(old),lambda:verified(new))
    for singular in ([[1,0],[0,0]],[[0,0],[0,0]],[[1,2],[2,4]]):
        row('singular-target-'+str(singular),'ERROR',
            lambda singular=singular:old.eps_tlin_affine_entry([0,0],identity,[0,0],singular),
            lambda singular=singular:new.eps_tlin_affine_entry([0,0],identity,[0,0],singular),True)
    for A in ([],[[1,0],[0]],[[1,2,3],[0,1,2]]):
        row('degenerate-input-'+str(A),'ERROR',lambda A=A:old.inv(A),lambda A=A:new.inv(A),True)
    for invalid in ([],[(0,F(1,2))],[(0,-1),(1,2)],[(1,1)]):
        row('invalid-law-'+str(invalid),'ERROR',lambda invalid=invalid:old.gamma1_sq(invalid),lambda invalid=invalid:new.gamma1_sq(invalid),True)
    xi=[(F(-1),F(1,3)),(F(0),F(1,3)),(F(1),F(1,3))]
    p=H.law([(x,) for x,w in xi],[w for x,w in xi])
    first=H.push(p,[[1]],[0]); second=H.push(p,[[2]],[F(3,10)])
    tail=lambda x:x+F(1,2)*x**3
    l1=H.law([(tail(x[0]),) for x in first],list(first.values()))
    l2=H.law([(tail(x[0]),) for x in second],list(second.values()))
    expected={'with_static_tail':H.skew_square(l2)-H.skew_square(l1),'tail_removed':F(0)}
    row('common-globally-invertible-static-tail',expected,lambda:old.static_map_control(xi),lambda:new.static_map_control(xi))
    tiny=F(1,10**40)
    near=H.law([(0,),(1,),(2+tiny,)],[1,1,1])
    row('near-zero-skew-is-not-zero',H.skew_square(near),
        lambda:old.gamma1_sq([(x[0],w) for x,w in near.items()]),
        lambda:new.gamma1_sq([(x[0],w) for x,w in near.items()]))
    return {'rows':rows,'comparison_count':len(rows),'mismatches':[x for x in rows if not x['agrees']],
            'original_defects':sum(not x['before_agrees'] for x in rows)}
