#!/usr/bin/env python3
"""Verify bibliographic identifiers used in 8986d34's LITERATURE_AUDIT.md / CSTAR_HYPOTHESIS.md
and in this review's NOTES.md.

Independent code (track cstar_lit). Queries the public arXiv export API and the Crossref API,
records title / authors / date / journal reference / DOI exactly as returned, and writes
results/refs.json. Network access is required; nothing is inferred offline.

Usage:  python3 -I code/verify_refs.py
"""
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NS = {'a': 'http://www.w3.org/2005/Atom', 'x': 'http://arxiv.org/schemas/atom'}

# (id, cited_by): 'codex' = cited in 8986d34 LITERATURE_AUDIT.md; 'review' = added by this review
ARXIV = [
    ('2002.08950', 'codex'), ('2402.00098', 'codex'), ('2403.13047', 'codex'),
    ('2501.02359', 'codex'), ('2602.03835', 'codex'), ('hep-th/0010252', 'codex'),
    ('0711.3170', 'codex'), ('1810.05337', 'codex'), ('1912.06993', 'codex'),
    # closed-universe one-state line, observer rules, dissent
    ('2004.06738', 'review'), ('2501.02632', 'review'), ('2503.09681', 'review'),
    ('2602.05939', 'review'), ('2505.20390', 'review'), ('2602.13387', 'review'),
    ('2507.10649', 'review'), ('2509.14327', 'review'), ('2504.14586', 'review'),
    ('2602.02682', 'review'), ('2505.14771', 'review'), ('2610.10667', 'review'),
    ('2207.06536', 'review'),
    # algebras / observers / clocks in closed or cosmological spacetimes
    ('2206.10780', 'review'), ('2308.03663', 'review'), ('2406.02116', 'review'),
    ('2410.00111', 'review'), ('2301.07257', 'review'),
    # finite-dimensional de Sitter Hilbert space / N-bound / covariant bound
    ('hep-th/0007146', 'review'), ('hep-th/0106109', 'review'), ('hep-th/0203101', 'review'),
    ('hep-th/9908070', 'review'), ('1404.5635', 'review'), ('hep-th/0102077', 'review'),
    # canonical GR, q-theory dynamics, symmetry protection, Gauss law in closed space
    ('gr-qc/0405109', 'review'), ('0806.2805', 'review'), ('1011.5120', 'review'),
    ('1810.05338', 'review'), ('hep-th/9502069', 'review'), ('0912.0225', 'review'),
    ('1508.06910', 'review'),
]

CROSSREF = [
    ('10.1103/PhysRevD.27.2885', 'codex', 'Page-Wootters 1983'),
    ('10.1103/PhysRev.160.1113', 'review', 'DeWitt 1967 canonical theory'),
    ('10.1016/0003-4916(74)90404-7', 'review', 'Regge-Teitelboim 1974 boundary terms'),
    ('10.1038/246396a0', 'review', 'Tryon 1973 zero-energy universe'),
    ('10.1103/PhysRevLett.62.1221', 'review', 'Krauss-Wilczek 1989 discrete gauge symmetry'),
]
CROSSREF_QUERIES = [
    ('Higuchi Quantum linearization instabilities of de Sitter spacetime 1991', 'review'),
    ('Moncrief Spacetime symmetries and linearization stability of the Einstein equations', 'review'),
    ('Arms Marsden Moncrief The structure of the space of solutions of Einstein equations II several Killing fields', 'review'),
]


def get(url, timeout=90):
    req = urllib.request.Request(url, headers={'User-Agent': 'grut-review-cstar-lit/1.0 (mailto:none)'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def arxiv_meta(ids):
    out = {}
    for i in range(0, len(ids), 20):
        chunk = ids[i:i + 20]
        url = 'https://export.arxiv.org/api/query?' + urllib.parse.urlencode(
            {'id_list': ','.join(chunk), 'max_results': len(chunk)})
        root = ET.fromstring(get(url))
        for e in root.findall('a:entry', NS):
            aid = e.find('a:id', NS).text.split('/abs/')[1]
            base = aid.rsplit('v', 1)[0]
            j = e.find('x:journal_ref', NS)
            d = e.find('x:doi', NS)
            c = e.find('x:comment', NS)
            out[base] = dict(
                latest_version=aid,
                title=' '.join(e.find('a:title', NS).text.split()),
                authors=[a.find('a:name', NS).text for a in e.findall('a:author', NS)],
                first_submitted=e.find('a:published', NS).text[:10],
                last_updated=e.find('a:updated', NS).text[:10],
                journal_ref=j.text if j is not None else None,
                doi=d.text if d is not None else None,
                comment=' '.join(c.text.split()) if c is not None else None,
            )
        time.sleep(3.5)
    return out


def crossref(doi):
    m = json.loads(get('https://api.crossref.org/works/' + urllib.parse.quote(doi, safe='/()')))['message']
    return dict(doi=doi, title=(m.get('title') or [''])[0],
                authors=[a.get('family', '') for a in m.get('author', [])],
                container=(m.get('container-title') or [''])[0], volume=m.get('volume'),
                page=m.get('page'), year=m.get('issued', {}).get('date-parts', [[None]])[0][0])


def crossref_query(q):
    url = 'https://api.crossref.org/works?' + urllib.parse.urlencode({'query.bibliographic': q, 'rows': 2})
    items = json.loads(get(url))['message']['items']
    return [dict(doi=m.get('DOI'), title=(m.get('title') or [''])[0],
                 authors=[a.get('family', '') for a in m.get('author', [])],
                 container=(m.get('container-title') or [''])[0], volume=m.get('volume'),
                 page=m.get('page'), year=m.get('issued', {}).get('date-parts', [[None]])[0][0]) for m in items]


def main():
    ids = [a for a, _ in ARXIV]
    meta = arxiv_meta(ids)
    rows = []
    for aid, who in ARXIV:
        r = dict(id=aid, cited_by=who, resolved=aid in meta)
        r.update(meta.get(aid, {}))
        rows.append(r)
    cr = []
    for doi, who, label in CROSSREF:
        try:
            x = crossref(doi); x.update(cited_by=who, label=label, resolved=True)
        except Exception as e:  # record failure rather than guess
            x = dict(doi=doi, cited_by=who, label=label, resolved=False, error=str(e))
        cr.append(x)
        time.sleep(1)
    cq = []
    for q, who in CROSSREF_QUERIES:
        try:
            cq.append(dict(query=q, cited_by=who, hits=crossref_query(q)))
        except Exception as e:
            cq.append(dict(query=q, cited_by=who, error=str(e)))
        time.sleep(1)
    out = dict(scope='Bibliographic identity only (title/authors/date/journal/DOI as returned by arXiv and Crossref). '
                     'Content claims are in NOTES.md with read levels.',
               arxiv=rows, crossref=cr, crossref_queries=cq,
               unresolved=[r['id'] for r in rows if not r['resolved']])
    (ROOT / 'results' / 'refs.json').write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n')
    for r in rows:
        print(('OK ' if r['resolved'] else 'MISSING ') + r['id'], '|', r.get('title', '')[:90], '|',
              r.get('journal_ref') or '-', '|', r.get('doi') or '-')
    for x in cr:
        print('CR', x['doi'], '|', x.get('title', x.get('error', ''))[:90], '|', x.get('container', ''), x.get('volume', ''), x.get('page', ''), x.get('year', ''))
    for q in cq:
        for h in q.get('hits', []):
            print('CQ', h['doi'], '|', h['title'][:90], '|', h['container'], h['volume'], h['page'], h['year'])


if __name__ == '__main__':
    main()
