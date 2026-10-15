"""Binary outcome entropy, including the unsuccessful preparation outcome."""
from math import log2


def binary_entropy(p):
    if p<=0 or p>=1:return 0.
    return -p*log2(p)-(1-p)*log2(1-p)


def include_binary_record_entropy(data):
    for row in data['rows']:
        lo,hi=row['formation_interval']
        endpoints=[binary_entropy(lo),binary_entropy(hi)]
        row['record_entropy_bits']=binary_entropy(row['formation_probability'])
        row['record_entropy_interval_bits']=[min(endpoints),
            1. if lo<=.5<=hi else max(endpoints)]
    return data
