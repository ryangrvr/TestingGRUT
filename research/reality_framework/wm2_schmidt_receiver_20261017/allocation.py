"""A fixed total Schmidt budget, shared across exact charge sectors."""


def global_keep(spectra,total):
    result={q:0 for q in spectra}
    entries=sorted(((float(value),q,j) for q,s in spectra.items() for j,value in enumerate(s)),reverse=True)
    for value,q,j in entries[:total]:result[q]+=1
    return result
