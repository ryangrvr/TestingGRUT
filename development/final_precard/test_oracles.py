"""Subject-free controls, frozen before importing the evaluated implementation."""
from fractions import Fraction as F
import itertools
import unittest
from . import price_oracle as P, hb_oracle as H, nr4_oracle as N
from .weight_boundary import PrevalidatedKit


class PriceReference(unittest.TestCase):
    def test_alphabet_and_reserved_slots(self):
        self.assertEqual(len(P.ALPHABET),64)
        self.assertEqual(len(set(P.ALPHABET)),64)
        self.assertEqual(len(P.RESERVED),4)
        for r in P.RESERVED:
            with self.subTest(r=r),self.assertRaises(ValueError): P.raw_price(r)

    def test_elias_exact_power_boundaries(self):
        for k in range(1,257):
            for n in (2**k-1,2**k,2**k+1):
                bits=len(format(n,'b'))
                # Direct length-of-length via base-two string, independent of
                # oracle's bit-length formula and any floating-point logarithm.
                expected=bits+2*len(format(bits,'b'))-2
                self.assertEqual(P.delta(n),expected)

    def test_all_operators_have_explicit_valid_payload(self):
        for op,arity in dict(P.ARITY,**P.EXPLICIT_SIGNATURES).items():
            args=['x','y','z'][:arity]
            if op in P.BINDERS: args=['x',('=','x','y')]
            expr=(op,*args)
            for precision in (6,10,16):
                value=P.raw_price(expr,precision,signatures=P.EXPLICIT_SIGNATURES)
                self.assertGreater(value['bits'],0)

    def test_wrong_argument_counts(self):
        for op,n in dict(P.ARITY,**P.EXPLICIT_SIGNATURES).items():
            for arity in (n+1, n-1):
                if arity<0: continue
                with self.subTest(op=op,arity=arity),self.assertRaises(ValueError):
                    P.raw_price((op,*(['x']*arity)),signatures=P.EXPLICIT_SIGNATURES)

    def test_malformed_unsupported_and_empty_ast(self):
        for expr in ((),[],None,3,False,('unpriced',),('nat',0,'extra'),
                     ('def',),('end',),('var','x'),('kernel','x','y'),
                     ('=',),('forall',('nat',2),'top'),('',)):
            with self.subTest(expr=expr),self.assertRaises(ValueError): P.raw_price(expr)
        self.assertEqual(P.raw_price('emptyset')['bits'],6)
        self.assertEqual(P.raw_price('top')['bits'],6)

    def test_integer_and_literal_domains(self):
        for x in (-1,1.5,True,float('nan'),float('inf')):
            with self.subTest(x=str(x)),self.assertRaises(ValueError): P.raw_price(('nat',x))
        for expr in (('rat',1,0),('rat',1,-2),('rat',0.5,2),('real',1.5)):
            with self.assertRaises(ValueError): P.raw_price(expr)
        for p in (-1,0,1.5,True):
            with self.assertRaises(ValueError): P.raw_price(('real',0),p)

    def test_repricing_exact_literal_controls(self):
        for p in (6,10,16):
            self.assertEqual(P.raw_price(('real',0),p)['bits'],6+p+2)
            self.assertEqual(P.raw_price(('=', 'x',('nat',0)),p)['bits'],20)
            self.assertEqual(P.raw_price(('rat',-1,2),p)['literal_bits'],9)
        self.assertEqual(P.literal_price('real',-7,16)-P.literal_price('real',-7,6),10)

    def test_definition_charge_and_cycles(self):
        defs=[('d',('=', 'x',('nat',3)))]
        row=P.raw_price(('and',('use','d'),('use','d')),definitions=defs)
        self.assertEqual(row['definition_uses'],2)
        self.assertEqual(row['definitions'],1)
        self.assertEqual(row['bits'],P.raw_price(defs[0][1])['bits']+6+2)
        for definitions in ([('d','top'),('d','bot')], [('d',('use','d'))],
                            [('d',('use','e')),('e',('use','d'))]):
            with self.assertRaises(ValueError): P.raw_price('top',definitions=definitions)
        with self.assertRaises(ValueError): P.raw_price(('use','absent'))

    def test_capture_and_free_variable_preservation(self):
        for expression in (
            ('forall','x',('=','x','q0')),
            ('and',('forall','x',('=','x','y')),('=','x','q0')),
            ('forall','x',('=','q0',('apply',('lambda','x','x'),'x'))),
        ):
            result=P.normalize(expression)
            self.assertEqual(P.free_variables(expression),P.free_variables(result))
            for values in itertools.product((0,1),repeat=4):
                env=dict(zip(('x','y','q0','q1'),values))
                self.assertEqual(bool(P.evaluate(expression,env)),bool(P.evaluate(result,env)))

    def test_prenex_nnf_semantics(self):
        atoms=[('=','x','q0'),('!=','x','q1'),('=','q0','q1'),'A','B']
        formulas=[]
        for op in ('and','or','->','<->'):
            formulas.extend((op,('forall','x',a),('exists','x',b))
                            for a,b in itertools.product(atoms,repeat=2))
        formulas.extend(('not',f) for f in formulas[:10])
        formulas.append(('ite','A',('forall','x',atoms[0]),('exists','x',atoms[1])))
        for expression in formulas:
            normal=P.normalize(expression)
            for values in itertools.product((0,1),repeat=5):
                env=dict(zip(('x','q0','q1','A','B'),values))
                self.assertEqual(bool(P.evaluate(expression,env)),bool(P.evaluate(normal,env)))
            self.assertEqual(P.free_variables(expression),P.free_variables(normal))

    def test_alpha_equivalent_normal_prices(self):
        a=('and',('forall','x',('=','x','y')),('exists','x',('=','x','q0')))
        b=('and',('forall','a',('=','a','y')),('exists','b',('=','b','q0')))
        self.assertEqual(P.normal_price(a),P.normal_price(b))

    def test_quantified_term_is_not_silently_hoisted(self):
        with self.assertRaises(ValueError): P.normalize(('=',('forall','x','A'),'B'))

    def test_ledger_and_hooks(self):
        self.assertEqual(P.ledger({'table':3,'inputs':4}),7)
        for v in (-1,0.5,True,float('nan')):
            with self.assertRaises(ValueError): P.ledger({'x':v})
        self.assertEqual(P.selection_price(6,1024),10)
        self.assertEqual(P.selection_price(100,4,relation_only=True),2)
        self.assertEqual(P.input_price(8,F(25,2)),F(25,2))
        self.assertEqual(P.table_price([('nat',0),('nat',1)]),17)
        with self.assertRaises(ValueError): P.equal_length([1,2],[1])
        self.assertEqual(P.equal_length([],[]),0)


class HostileBathReference(unittest.TestCase):
    def setUp(self):
        self.p=H.law([(0,0),(1,0),(0,1)], [1,2,3])

    def test_invertible_affine_exact_zero_certificate(self):
        for A,b in (([[1,0],[0,1]],[0,0]),([[2,1],[-1,3]],[5,-3]),
                    ([[-1,0],[0,F(1,2)]],[0,1])):
            q=H.push(self.p,A,b)
            self.assertTrue(H.affine_certificate(self.p,q,A,b))

    def test_representation_covariance(self):
        A=[[2,1],[-1,3]]; b=[5,-3]
        q=H.push(self.p,A,b)
        # Applying a common invertible coordinate change R to both sides.
        R=[[1,1],[0,1]]; Rp=H.push(self.p,R,[0,0]); Rq=H.push(q,R,[0,0])
        conjugate=[[1,3],[-1,4]]; translated=[2,-3]
        self.assertTrue(H.affine_certificate(Rp,Rq,conjugate,translated))

    def test_singular_and_degenerate_covariance_reject(self):
        for A in ([[1,0],[0,0]],[[0,0],[0,0]],[[1,2],[2,4]]):
            with self.assertRaises(ValueError): H.affine_certificate(self.p,H.push(self.p,A,[0,0]),A,[0,0])
        degenerate=H.law([(0,0),(1,0)],[1,1])
        with self.assertRaises(ValueError): H.affine_certificate(degenerate,degenerate,[[1,0],[0,1]],[0,0])

    def test_domain_and_finite_chart_validation(self):
        with self.assertRaises(ValueError): H.affine_certificate(self.p,self.p,[[1,0],[0,1]],[0,0],source_domain={(0,0)})
        for points,weights in (([],[]),([(0,),(1,)],[1]),([(0,),(1,2)],[1,1]),([(0,)],[0])):
            with self.assertRaises(ValueError): H.law(points,weights)

    def test_static_tail_exact_witness(self):
        p=H.law([(0,),(1,),(2,)],[1,1,1]); shift=H.push(p,[[1]],[3])
        self.assertTrue(H.affine_certificate(p,shift,[[1]],[3]))
        p_squared=H.law([(x[0]**2,) for x in p],list(p.values()))
        q_squared=H.law([(x[0]**2,) for x in shift],list(shift.values()))
        self.assertTrue(H.scalar_control([p_squared,q_squared])['nonzero_witness'])

    def test_near_zero_is_not_exact_zero(self):
        tiny=F(1,10**40)
        p=H.law([(0,),(1,),(2,)],[1,1,1])
        q=H.law([(0,),(1,),(2+tiny,)],[1,1,1])
        self.assertFalse(H.affine_certificate(p,q,[[1]],[0]))
        self.assertTrue(H.scalar_control([p,q])['nonzero_witness'])


class AblationReference(unittest.TestCase):
    def test_full_K_mask_and_missing_image(self):
        mask=N.full_mask([{'x':0},{'x':2}],[{'x':0},{'x':0}],{'x':1})
        self.assertEqual(mask,[True,False])
        result=N.runs([1,1],mask,{'K':[False,True],'empty':[False,False]},F(1,2),['a','a'])
        self.assertEqual(result['K']['fraction'],1)
        self.assertEqual(result['empty']['fraction'],0)

    def test_zero_weight_fiber_and_all_excluded_VOID(self):
        result=N.runs([1,0],[True,False],{'K':[True,True]},F(1,2),['a','b'])['K']
        self.assertIsNone(result['fraction'])
        self.assertEqual(result['per_fiber'],{})
        self.assertIsNone(result['appears'])

    def test_reference_measure_contract(self):
        for weights in ([-1,2],[2,-1],[0,0],[],[1,float('nan')]):
            with self.assertRaises((ValueError,OverflowError)):
                N.measure(weights,[False]*len(weights),[True]*len(weights))
        self.assertEqual(N.measure([0,2],[False,False],[True,False])['fraction'],0)

    def test_unresolved_and_mismatched_inputs(self):
        for images,refs,tol in (([None],[{'x':0}],{'x':1}),([{'x':1}],[],{'x':1}),
                                ([{'y':0}],[{'x':0}],{'x':1})):
            with self.assertRaises(ValueError): N.full_mask(images,refs,tol)

    def test_boundary_order_before_all_callbacks(self):
        calls=[]
        class Subject:
            def appearance(self,*a,**kw): calls.append('appearance')
            def nr4(self,*a,**kw): calls.append('nr4')
        wrapper=PrevalidatedKit(Subject())
        for weights in ([-1,2],[2,-1],[0,0],[],[1,float('nan')],[float('inf'),1]):
            points=[{'xi':i,'w':w} for i,w in enumerate(weights)]
            with self.assertRaises(ValueError): wrapper.appearance('K',points)
            with self.assertRaises(ValueError): wrapper.nr4('K','E','S',{},points)
        self.assertEqual(calls,[])


if __name__=='__main__': unittest.main()
