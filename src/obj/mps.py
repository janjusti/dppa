class MyProtSeq:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        if 'consensus sequence' in name:
            self.type = 'consensus'
        else:
            self.type = 'pure'