

'''
Primer name pair1
Amplimer 1
	Sequence: CP072761.1  
	Mycobacterium tuberculosis strain CG24 chromosome, complete genome
	CGAACGGAAAGGCCCCTTCG hits forward strand at 1471917 with 2 mismatches
	CCGTCGTCGCCTTGGTAG hits reverse strand at [2987315] with 0 mismatches
	Amplimer length: 219 bp

Primer name pair2
Amplimer 1
	Sequence: CP072761.1  
	Mycobacterium tuberculosis strain CG24 chromosome, complete genome
	CGAACGGAAAGGTCTCTTCG hits forward strand at 1471917 with 0 mismatches
	CCGTCGTCGCCTTGGTAG hits reverse strand at [2987315] with 0 mismatches
	Amplimer length: 219 bp
'''


class Amplifier:
    """Represent a single amplification from a primer."""

    def __init__(self):
        """Initialize the class."""
        self.hit_info = ""
        self.length = 0
        self.forward_primer = None 
        self.reverse_primer = None
        self.forward_n_missmatch = None
        self.reverse_n_missmatch = None
        self.start_position = None
        
    def __str__(self,):
        pstr = f'''{self.hit_info}:
Match start position:\t{self.start_position}
Forward-primer:\t{self.forward_primer}\t{self.forward_n_missmatch}
Reverse-primer:\t{self.reverse_primer}\t{self.reverse_n_missmatch}
        '''
        return pstr
        
        

class PrimersearchOutputRecord:
    """Represent the information from a primersearch.
    """

    def __init__(self):
        """Initialize the class."""
        self.amplifiers = {}
    def __str__(self,):
        return f'Number of matches: {len(self.amplifiers)}: {",".join(self.amplifiers)}'


def read_primersearch(handle):
    """Get output from primersearch into a PrimerSearchOutputRecord."""
    import re
    record = PrimersearchOutputRecord()
    rf = '([A-Z]+) hits forward strand at ([0-9]+) with ([0-9]+) mismatches'
    rr = '([A-Z]+) hits reverse strand at \[[0-9]+\] with ([0-9]+) mismatches'

    for line in handle:
        if not line.strip():
            continue
        elif line.startswith("Primer name"):
            name = line.split()[-1]
            record.amplifiers[name] = []
        elif line.startswith("Amplimer"):
            amplifier = Amplifier()
            record.amplifiers[name].append(amplifier)
        elif line.startswith("\tSequence: "):
            amplifier.hit_info = line.replace("\tSequence: ", "")
        elif line.startswith("\tAmplimer length: "):
            length = line.split()[-2]
            amplifier.length = int(length)
        elif 'forward strand' in line:
            # CGAACGGAAAGGTCTCTTCG hits forward strand at 1471917 with 0 mismatches
            s = re.search(rf, line)
            amplifier.forward_primer = s.group(1)
            amplifier.start_position = s.group(2)
            amplifier.forward_n_missmatch = s.group(3)
        elif 'reverse strand' in line:
            s = re.search(rr, line)
            amplifier.reverse_primer = s.group(1)
            amplifier.reverse_n_missmatch = s.group(2)
        else:
            amplifier.hit_info += line

    for name in record.amplifiers:
        for amplifier in record.amplifiers[name]:
            amplifier.hit_info = amplifier.hit_info.rstrip()
    
    return record