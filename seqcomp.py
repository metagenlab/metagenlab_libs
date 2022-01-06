

def hamming(seq1, seq2):
    import Levenshtein as lv
    
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    
    closest = 1000000000
    
    for n in range(0, len_seq2):
        stop = n+len_seq1
        if stop > len_seq2:
            return closest
        subseq = seq2[n:stop]
        dist = lv.hamming(subseq, seq1)
        if dist < closest:
            closest = dist 
