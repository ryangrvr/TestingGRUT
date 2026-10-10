"""Independent Appendix-B arithmetic and scoped logic reference.

Written from frozen charter §4.1–4.2 and Appendix B, before subject comparison.
This is a review DSL: ASCII spellings stand for the printed tokens. The two ×
tokens occupy distinct set/arithmetic slots. It is not a type checker or a proof
of every IP gate. Conditional/process arities are explicit fixture contracts;
Appendix B does not freeze their AST signatures. Unsupported contracts reject.
No kit implementation, price function, or candidate is imported.
"""
from fractions import Fraction
from math import log2, isfinite

GROUPS = (
    'not and or -> <-> forall exists = != top bot ite',
    'in subseteq emptyset singleton pair times powerset card union inter minus range',
    'lambda apply compose id image preimage inverse restrict',
    '0 1 succ + - * / <= < sum_finite',
    'kernel cond_prob E tensor marginal do law support',
    'sequential parallel contract wire',
    'def var nat rat real end',
    'reserved1 reserved2 reserved3 reserved4',
)
ALPHABET = tuple(t for group in GROUPS for t in group.split())
RESERVED = frozenset(ALPHABET[-4:])
ARITY = dict(zip(('not','and','or','->','<->','forall','exists','=','!=',
    'top','bot','ite','in','subseteq','emptyset','singleton','pair','times',
    'powerset','card','union','inter','minus','range','lambda','apply','compose',
    'id','image','preimage','inverse','restrict','0','1','succ','+','-','*','/',
    '<=','<','sum_finite'), (1,2,2,2,2,2,2,2,2,0,0,3,2,2,0,1,2,2,1,1,2,2,2,1,
    2,2,2,0,2,2,1,2,0,0,1,2,2,2,2,2,2,2)))
EXPLICIT_SIGNATURES = {'kernel':3, 'cond_prob':3, 'E':2, 'tensor':2,
    'marginal':2, 'do':2, 'law':1, 'support':1, 'sequential':2,
    'parallel':2, 'contract':2, 'wire':2}
BINDERS = frozenset(('forall','exists','lambda'))


def integer(n, minimum=0):
    if type(n) is not int or n < minimum:
        raise ValueError('integer outside required domain')
    return n


def delta(n):
    """Elias δ: gamma-coded bit length followed by the integer's binary tail."""
    integer(n, 1)
    bits = bin(n)[2:]
    length_bits = bin(len(bits))[2:]
    return (len(length_bits)-1) + len(length_bits) + (len(bits)-1)


def width(n):
    integer(n)
    return n.bit_length()  # ceil(log2(n+1)), including width(0)=0


def definition_map(entries=()):
    if isinstance(entries, dict):
        entries = list(entries.items())
    out = {}
    for row in entries:
        if not isinstance(row, (list,tuple)) or len(row) != 2:
            raise ValueError('definition must be a name/body pair')
        name, body = row
        if not isinstance(name,str) or not name or name in ALPHABET or name in out:
            raise ValueError('duplicate/invalid definition name')
        out[name] = body
    return out


def contracts(signatures):
    sigs = {} if signatures is None else dict(signatures)
    for name, arity in sigs.items():
        if name not in EXPLICIT_SIGNATURES or integer(arity) != EXPLICIT_SIGNATURES[name]:
            raise ValueError('unsupported/unfrozen operator signature')
    return dict(ARITY, **sigs)


def validate(expr, definitions=(), signatures=None):
    defs, arities = definition_map(definitions), contracts(signatures)
    def visit(e):
        if isinstance(e,str):
            if not e or e in RESERVED or e in ('def','var','nat','rat','real','end','use'):
                raise ValueError('invalid bare token')
            if e in ALPHABET and arities.get(e) != 0:
                raise ValueError('operator requires operands/signature')
            return set()
        if not isinstance(e,(tuple,list)) or not e or not isinstance(e[0],str):
            raise ValueError('malformed AST')
        op, *args = e
        if op == 'use':
            if len(args)!=1 or not isinstance(args[0],str) or args[0] not in defs:
                raise ValueError('undefined/malformed symbol use')
            return {args[0]}
        if op == 'nat':
            if len(args)!=1: raise ValueError('natural payload arity')
            integer(args[0]); return set()
        if op == 'rat':
            if len(args)!=2 or type(args[0]) is not int: raise ValueError('rational payload')
            integer(args[1],1); return set()
        if op == 'real':
            if len(args)!=1 or type(args[0]) is not int: raise ValueError('real exponent payload')
            return set()
        if op not in arities or len(args)!=arities[op]:
            raise ValueError('unsupported symbol or wrong arity')
        if op in BINDERS and (not isinstance(args[0],str) or not args[0] or args[0] in ALPHABET):
            raise ValueError('binder must be a variable name')
        refs=set()
        for arg in args: refs.update(visit(arg))
        return refs
    graph={name:visit(body) for name,body in defs.items()}
    visit(expr)
    done, active = set(), set()
    def acyclic(name):
        if name in active: raise ValueError('definition cycle')
        if name in done: return
        active.add(name)
        for dependency in graph[name]: acyclic(dependency)
        active.remove(name); done.add(name)
    for name in defs: acyclic(name)
    return defs


def free_variables(e, bound=frozenset()):
    if isinstance(e,str):
        return set() if e in ALPHABET or e in bound else {e}
    op,*args=e
    if op in ('nat','rat','real','use'): return set()
    if op in BINDERS: return free_variables(args[1], bound | {args[0]})
    return set().union(*(free_variables(a,bound) for a in args))


def normalize(expr, signatures=None):
    """Capture-free prenex NNF, on nonempty common first-order domains.

    Standardization apart happens twice: first before NNF duplication, then
    after implication/equivalence/ite elimination duplicates subformulas.
    Lambda binders stay scoped in terms and never enter the quantifier prefix.
    Quantifiers buried under a term operator are rejected, not hoisted through
    an arbitrary kernel/map. Such a hoist is not a first-order equivalence.
    """
    validate(expr,signatures=signatures)
    used=set(free_variables(expr)); counter=0
    def fresh():
        nonlocal counter
        while True:
            name='q'+str(counter); counter+=1
            if name not in used: used.add(name); return name
    def apart(e, env):
        if isinstance(e,str): return env.get(e,e)
        op,*args=e
        if op in ('nat','rat','real','use'): return tuple(e)
        if op in BINDERS:
            name=fresh(); extended=dict(env); extended[args[0]]=name
            return (op,name,apart(args[1],extended))
        return (op,*(apart(a,env) for a in args))
    def nnf(e, negative=False):
        if isinstance(e,str):
            if e=='top': return 'bot' if negative else 'top'
            if e=='bot': return 'top' if negative else 'bot'
            return ('not',e) if negative else e
        op,*args=e
        if op=='not': return nnf(args[0],not negative)
        if op=='->': return nnf(('or',('not',args[0]),args[1]),negative)
        if op=='<->':
            return nnf(('and',('->',args[0],args[1]),('->',args[1],args[0])),negative)
        if op=='ite':
            return nnf(('or',('and',args[0],args[1]),('and',('not',args[0]),args[2])),negative)
        if op in ('and','or'):
            head=({'and':'or','or':'and'}[op] if negative else op)
            return (head,nnf(args[0],negative),nnf(args[1],negative))
        if op in ('forall','exists'):
            head=({'forall':'exists','exists':'forall'}[op] if negative else op)
            return (head,args[0],nnf(args[1],negative))
        return ('not',e) if negative else e
    def has_quantifier(e):
        if isinstance(e,str): return False
        if e[0] in ('nat','rat','real','use'): return False
        return e[0] in ('forall','exists') or any(has_quantifier(a) for a in e[1:])
    def extract(e):
        if isinstance(e,str): return [],e
        op,*args=e
        if op in ('forall','exists'):
            prefix,body=extract(args[1]); return [(op,args[0])]+prefix,body
        if op in ('and','or'):
            lp,l=extract(args[0]); rp,r=extract(args[1])
            return lp+rp,(op,l,r)
        if has_quantifier(e): raise ValueError('unsupported quantified term/formula placement')
        return [],e
    normal=apart(nnf(apart(expr,{})),{})
    prefix,body=extract(normal)
    for op,var in reversed(prefix): body=(op,var,body)
    assert free_variables(expr)==free_variables(body)
    return body


def raw_price(expr, p=10, definitions=(), signatures=None):
    """Price a supplied syntax form; normal-form/equivalence proof is separate.

    Signed rational numerator uses its magnitude plus the fixed rational sign
    bit. Definitions' bodies are charged once each; uses cost the menu index
    width, with no extra base-token charge for the review-DSL 'use' marker.
    """
    integer(p,1); defs=validate(expr,definitions,signatures)
    variables=[]; tokens=0; literal=0; uses=0
    def walk(e):
        nonlocal tokens,literal,uses
        if isinstance(e,str):
            tokens+=1
            if e not in ALPHABET: variables.append(e)
            return
        op,*args=e
        if op=='use': uses+=1; return
        tokens+=1
        if op=='nat': literal+=delta(args[0]+1)
        elif op=='rat': literal+=delta(abs(args[0])+1)+delta(args[1]+1)+1
        elif op=='real': literal+=p+delta(abs(args[0])+1)+1
        else:
            for a in args: walk(a)
    walk(expr)
    for body in defs.values(): walk(body)
    return {'bits':6*tokens+len(variables)*width(len(set(variables)))+literal+uses*width(len(defs)),
            'tokens':tokens,'variable_occurrences':len(variables),
            'distinct_variables':len(set(variables)),'literal_bits':literal,
            'definition_uses':uses,'definitions':len(defs)}


def normal_price(expr,p=10,signatures=None):
    return raw_price(normalize(expr,signatures),p,signatures=signatures)


def nonnegative(value):
    try: q=Fraction(value)
    except (ValueError,TypeError,OverflowError,ZeroDivisionError) as exc:
        raise ValueError('finite numeric amount required') from exc
    if q<0: raise ValueError('negative amount')
    return q


def selection_price(statement_bits, family_size, relation_only=False):
    integer(statement_bits); integer(family_size,1)
    choice=log2(family_size)
    return choice if relation_only else max(statement_bits,choice)


def input_price(statement_bits, kl_bits=0):
    integer(statement_bits)
    return max(statement_bits,nonnegative(kl_bits))


def table_price(entries,p=10):
    # Every hand-listed entry supplied through this hook is charged. The caller
    # must additionally price the table's structural declaration, never inferred.
    return sum(raw_price(e,p)['bits'] for e in entries)


def literal_price(kind,value,p=10):
    e=(kind,*value) if isinstance(value,tuple) else (kind,value)
    return raw_price(e,p)['literal_bits']


def ledger(fields):
    if not isinstance(fields,dict) or not fields: raise ValueError('nonempty ledger required')
    return sum(integer(v) for v in fields.values())


def equal_length(*arrays):
    if not arrays or any(not isinstance(a,(list,tuple)) for a in arrays):
        raise ValueError('arrays required')
    if len({len(a) for a in arrays})!=1: raise ValueError('array length mismatch')
    return len(arrays[0])


def evaluate(e,env,domain=(0,1)):
    """Finite semantic control only; no physical state or law under evaluation."""
    if isinstance(e,str):
        if e in ('top','bot','0','1'): return {'top':True,'bot':False,'0':0,'1':1}[e]
        return env[e]
    op,*a=e
    if op in ('forall','exists'):
        values=[evaluate(a[1],dict(env,**{a[0]:x}),domain) for x in domain]
        return all(values) if op=='forall' else any(values)
    if op=='lambda': return lambda x:evaluate(a[1],dict(env,**{a[0]:x}),domain)
    if op=='apply': return evaluate(a[0],env,domain)(evaluate(a[1],env,domain))
    if op=='nat': return a[0]
    if op=='rat': return Fraction(*a)
    if op=='not': return not evaluate(a[0],env,domain)
    if op=='ite': return evaluate(a[1] if evaluate(a[0],env,domain) else a[2],env,domain)
    if not a: return evaluate(op,env,domain)
    x,y=(evaluate(v,env,domain) for v in a)
    return {'and':lambda:bool(x) and bool(y),'or':lambda:bool(x) or bool(y),
        '->':lambda:not bool(x) or bool(y),'<->':lambda:bool(x)==bool(y),
        '=':lambda:x==y,'!=':lambda:x!=y,'<':lambda:x<y,'<=':lambda:x<=y,
        '+':lambda:x+y,'-':lambda:x-y,'*':lambda:x*y,'/':lambda:x/y}[op]()
