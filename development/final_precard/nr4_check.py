"""Historical regression, certified-mask audit and prevalidation hostiles."""
from fractions import Fraction as F
import itertools
from . import nr4_oracle as N
from .weight_boundary import PrevalidatedKit


def regression(pre,post,D,cr5):
    points=[{'xi':i,'w':1,'fiber':'a'} for i in range(3)]
    def pipeline(law,i):
        if law=='K': return {'x':i}
        if law in ('empty','deleted'): return {'x':1}
        return None
    holds=lambda image,tol:image['x']>0
    before=pre.nr4('K','empty','standard',{'c':'deleted'},points,pipeline,lambda _:0,holds,{'x':1},2)
    after=post.nr4('K','empty','standard',{'c':'deleted'},points,pipeline,lambda _:0,holds,{'x':1},2)
    binding=D.Binding('fixture-full-K','context','zero','domain','metric','chart','resolution','fiber',('x',))
    registry={'reference':D.ExactZeroProof(binding,'synthetic-reference',{'x':D.Interval(0)})}
    for i in range(3): registry['image-'+str(i)]=D.ImageProof(binding,'synthetic-image',str(i),{'x':D.Interval(i)})
    calls=[]
    def certificate(pt):
        calls.append(pt['xi'])
        return D.CornerEvidence('image-'+str(pt['xi']),'reference',binding,binding,str(pt['xi']))
    masks=[]; original=cr5.appearance
    def observed(*args,**kwargs):
        masks.append(id(args[5]))
        return original(*args,**kwargs)
    cr5.appearance=observed
    try:
        certified=PrevalidatedKit(cr5).nr4('K','empty','standard',{'c':'deleted'},points,pipeline,
            lambda _:(_ for _ in ()).throw(AssertionError('zero-sampling forbidden')),
            holds,{'x':1},2,certificates=certificate,approved_registry=registry,full_pipeline_digest='fixture-full-K')
    finally: cr5.appearance=original
    mask=N.full_mask([pipeline('K',i) for i in range(3)],[{'x':0}]*3,{'x':F(1,1024)})
    reference=N.runs([1]*3,mask,{'K':[False,True,True],'K_empty':[True]*3,
        'K_S':[False]*3,'K-c':[True]*3},F(1,2),['a']*3)
    rows=[{'fixture':'historical-pre-fix-false-negative','expected':'NOT-RELOCATED-BY-NR-4','actual':before['verdict']},
          {'fixture':'CR4-post-fix-relocation','expected':'RELOCATED','actual':after['verdict']},
          {'fixture':'CR5-certified-relocation','expected':'RELOCATED','actual':certified['verdict']},
          {'fixture':'full-K-shared-mask','expected':mask,'actual':certified['shared_mask']},
          {'fixture':'certificate-call-count','expected':list(range(3)),'actual':calls},
          {'fixture':'one-identical-mask-object-for-all-four-runs','expected':True,'actual':len(masks)==4 and len(set(masks))==1},
          {'fixture':'responsibility-not-repaired-by-variant-exclusion','expected':False,'actual':certified['responsibility']['c']}]
    for name,value in reference.items():
        rows.append({'fixture':'independent-fraction-'+name,'expected':value['fraction'],'actual':certified['runs'][name]['fraction']})
    for r in rows: r['agrees']=r['expected']==r['actual']
    return {'rows':rows,'mismatches':[r for r in rows if not r['agrees']],
            'pre_fix_regression_exhibits_defect':before['verdict']!='RELOCATED',
            'synthetic_receipts_only':True,'physical_certificate':False}


def weights(subject):
    safe=PrevalidatedKit(subject)
    invalid=[('signed-positive-total',[-1,2]),('negative-second',[2,-1]),
        ('signed-zero-total',[-1,1]),('zero-reference-measure',[0,0]),
        ('empty-reference-measure',[]),('nan-second',[1,float('nan')]),
        ('positive-infinity-second',[1,float('inf')]),('negative-infinity-second',[1,float('-inf')]),
        ('nonfinite-string-second',[1,'NaN']),('nan-first',[float('nan'),1]),
        ('positive-infinity-first',[float('inf'),1]),('negative-infinity-first',[float('-inf'),1]),
        ('invalid-late-after-visible-callback',[1,1,'invalid']),
        ('negative-last-of-three',[1,2,-1]),('single-zero',[0]),('negative-single',[-1])]
    rows=[]
    for name,ws in invalid:
        pt=[{'xi':i,'w':w} for i,w in enumerate(ws)]
        for entry in ('appearance','nr4'):
            observations=[]
            for implementation in (subject,safe):
                calls=[]; error=None; result=None
                def pipeline(law,i): calls.append(('pipeline',i)); return {'x':i}
                def holds(*_): calls.append(('holds',None)); return True
                def certificate(*_): calls.append(('certificate',None)); return None
                try:
                    if entry=='appearance':
                        result=implementation.appearance('K',pt,pipeline,holds,{'x':1},[False]*len(pt))
                    else:
                        result=implementation.nr4('K','empty','standard',{},pt,pipeline,None,holds,{'x':1},2,
                            certificates=certificate,approved_registry={},full_pipeline_digest='fixture')
                except Exception as exc: error=type(exc).__name__
                observations.append({'error':error,'calls':calls,
                    'fraction':None if result is None else result.get('fraction'),
                    'verdict':None if result is None else result.get('verdict')})
            rows.append({'fixture':name+'-'+entry,'weights':[str(w) for w in ws],
                'before':observations[0],'after':observations[1],
                'agrees':observations[1]['error']=='ValueError' and not observations[1]['calls']})
    valid=[]
    for ws in itertools.product((F(0),F(1),F(2),F(1,3)),repeat=2):
        if sum(ws)==0: continue
        for mask,hits in itertools.product(itertools.product((False,True),repeat=2),repeat=2):
            pt=[{'xi':i,'w':w,'fiber':str(i)} for i,w in enumerate(ws)]
            pipeline=lambda law,i:{'hit':hits[i]}
            actual=safe.appearance('K',pt,pipeline,lambda image,_:image['hit'],{'x':1},list(mask))
            expected=N.measure(ws,list(mask),list(hits))
            if actual['fraction']!=expected['fraction'] or actual['excluded_fraction']!=expected['excluded_fraction']:
                raise AssertionError('validated boundary changed valid weighted arithmetic')
            if actual['fraction'] is not None and not 0<=actual['fraction']<=1:
                raise AssertionError('fraction escaped [0,1]')
            valid.append({'weights':list(ws),'mask':list(mask),'hits':list(hits),'fraction':actual['fraction'],
                          'excluded_fraction':actual['excluded_fraction'],'agrees':True})
    # Fraction zero denominators remain VOID; no unmasking of zero-mass fibers.
    return {'hostile_rows':rows,'invalid_weight_lists':len(invalid),'entrypoint_checks':len(rows),
            'mismatches':[r for r in rows if not r['agrees']],
            'additional_valid_measure_comparisons':len(valid),'valid_rows':valid,
            'private_source_modified':False,'isolation_boundary':'PrevalidatedKit ONLY',
            'direct_original_module_cleared':False,'external_approval':False}
