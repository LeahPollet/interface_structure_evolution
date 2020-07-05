"""
Given query sequence, subject sequence, and subject annotations,
align query to subject, then transfer annotations to query;
gapped positions in the query will be ignored.

Example:
q = AAAAGTTTT
s = AAAATTTTC
p = 012345679

Alignment:
q = AAAAGTTTT-
s = AAAA-TTTTC
p = 0123-45678

Output:
q = AAAAGTTTT
s = AAAA-TTTT
p = 0123-4567 ***

"""

import os
import sys
import re
import subprocess
import tempfile


def align(query, subject, pattern):
    """Given query sequence, subject sequence, and subject annotations,
    align query to subject, then transfer annotations to query;
    gapped positions in the query will be ignored.

    Args:
        query (): aa sequence to align to.
        subject (): aa sequence of interest.
        pattern (): residue-level annotations for protein of interest.

    Returns:
        : aligned annotations

    """
    
    # Process pattern
    if type(pattern) is not list:
        pattern = re.split('\s+', pattern) if ' ' in pattern else [k for k in pattern]
    pgap = "-----"
    if len(subject) != len(pattern):
        raise ValueError('subject and pattern are not the same length')
    # Write the query and subject sequences to memory
    tmp = tempfile.NamedTemporaryFile(delete=False, prefix='query')
    with open (tmp.name,'w') as f:
        queryPath = f.name
        f.write('>query\n')
        f.write(query + '\n')
    tmp = tempfile.NamedTemporaryFile(delete=False, prefix='subject')    
    with open (tmp.name,'w') as g:
        subjectPath = g.name
        g.write('>subject\n')
        g.write(subject + '\n')

    # Execute blast alignment (bl2seq equivalent)
    blastCmd = ('blastp' +
                ' -query ' + queryPath +
                ' -subject ' + subjectPath +
                ' -outfmt \"6 qlen slen qstart qend sstart send qseq sseq\"')
    cmd = subprocess.run(blastCmd, shell=True, stdout=subprocess.PIPE)
    results = cmd.stdout.decode('utf-8')
    results = results.split('\n')[0].split('\t')
    os.remove(queryPath)
    os.remove(subjectPath)
    # Format blast result
    if results == ['']:
        raise UserWarning('BLAST alignment failed')
    qlen, slen, qstart, qend, sstart, send = [int(k) for k in results[0:-2]]
    qseq, sseq = results[-2:]
    alen = len(qseq)
    # Use alignment result to align pattern
    pseq = []
    qindex = qstart - 1
    sindex = sstart - 1
    for aindex in range(alen):
        # No query, skip subject/pattern info
        if qseq[aindex] == "-":
            qindex += 0
            sindex += 1
        # No subject, gap
        elif sseq[aindex] == "-":
            pseq.append(pgap)
            qindex += 1
            sindex += 0
        # Transfer pattern
        else:
            pseq.append(pattern[sindex])
            qindex += 1
            sindex += 1
    # Add opening / closing gaps
    opening_gaps = [pgap] * (qstart - 1)
    closing_gaps = [pgap] * (qlen - qend)
    pseq = opening_gaps + pseq + closing_gaps
    if len(pseq) != len(query):
        raise UserWarning('query and aligned pattern not of the same length')
    
    return pseq


def main():
    q = 'AAAAGTTTT'
    s = 'AAAATTTTC'
    p = '012345679'
    print (align(q, s, p))


if __name__ == "__main__":
    main()
