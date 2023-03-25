import random
import myerror
from collections import defaultdict
from pprint import pformat
import copy


def pairing(strand1):
    '''配对'''
    pw = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    strand2 = ''
    for i in strand1:
        strand2 += pw[i]
    return strand2  # 从3'到5'


def dict_append(item, dict_):
    '''字典添加,已弃用'''
    try:
        dict_[item] += 1
    except:
        dict_[item] = 1


def random_DNA(n):
    '''生成随机DNA'''
    baselist = ''
    m = 'ATCG'
    for i in range(n):
        baselist += m[random.randint(0, 3)]
    return DNA.init_with_single(baselist)


class DNA:
    '''DNA类'''
    type = 'DNA'

    def __init__(self, template_strand, sense_strand='', start_sen=0, been_anneal=False):
        self.template_strand = template_strand  # 模板链
        self.sense_strand = sense_strand  # 模板链
        self.start_sen = start_sen  # 模板链起始位置
        if self.sense_strand == '':
            self.is_single = True
        else:
            self.is_single = False
        self.been_anneal = been_anneal

    @classmethod
    def init_with_single(self, template_strand):
        '使用类方法自动创建双链DNA'
        sense_strand = pairing(template_strand)
        return DNA(template_strand, sense_strand)

    def __str__(self):
        if self.sense_strand != '':
            a = "5'-%s-3'\n3'-%s-5'" % (self.template_strand,
                                        ' '*self.start_sen+self.sense_strand)
        else:
            a = "5'-%s-3'" % (self.template_strand)
        return a

    def __repr__(self) -> str:
        if self.sense_strand != '':
            a = "5'-%s-3' 3'-%s-5'" % (self.template_strand, self.sense_strand)
        else:
            a = "5'-%s-3'" % (self.template_strand)
        return a

    def __len__(self):
        return max(len(self.sense_strand)+self.start_sen, len(self.template_strand))

    def __eq__(self, __o: object) -> bool:
        if type(__o) != type(self):
            return False
        else:
            if (self.template_strand, self.sense_strand, self.start_sen) == (__o.template_strand, __o.sense_strand, __o.start_sen):
                return True

            else:
                return False

    def __hash__(self) -> int:
        c = self.template_strand
        d = self.sense_strand
        if c > d:
            c, d = d, c
        a = hash(c+d)
        return a

    def self_checking(self):
        '''调转DNA链'''
        if self.template_strand > self.sense_strand:
            self.template_strand, self.sense_strand = self.sense_strand, self.template_strand

    def count_hydrogen_bonds(self):
        '''氢键计数'''
        n = 0
        for i in self.template_strand[self.start_sen:min(len(self.template_strand), self.start_sen+len(self.sense_strand))]:
            if i in 'AT':
                n += 2
            else:
                n += 3
        return n

    def generate_RNA(self):
        '''转录形成RNA'''
        pw = {'A': 'T', 'C': 'U', 'T': 'A', 'G': 'C'}
        strand2 = ''
        for i in self.template_strand:
            strand2 += pw[i]
        return RNA(strand2)

    def dna_denature(self):
        '''变性'''
        if self.is_single is not True:
            return DNA(self.template_strand), DNA(self.sense_strand[::-1])
        else:
            return (self,)

    def dna_pairing(self, another_dna):
        '''主动与其他DNA配对'''
        if self.is_single and another_dna.is_single is True:
            position = another_dna.template_strand.find(
                pairing(self.template_strand)[::-1])
        else:
            raise myerror.PairingError('只能暂时配对单链')
        if position == -1:
            raise myerror.PairingError('这两段不能配对')
        return DNA(another_dna.template_strand, self.template_strand[::-1], start_sen=position)


class Primer(DNA):
    def __init__(self, template_strand, sense_strand='', start_sen=0):
        super().__init__(template_strand, sense_strand, start_sen)
        self.type = 'primer'

    def tm_caculater(self):
        # Tm=4℃ (G + C)+ 2℃ (A + T)
        _ = defaultdict(int)
        for i in self.template_strand:
            _[i] += 1
        return 4*(_['G']+_['C'])+2*(_['A']+_['T'])


class RNA:
    '''RNA类'''

    def __init__(self, strand):
        self.strand = strand[::-1]

    def __str__(self):
        return "5'-%s-3'" % self.strand

    def __repr__(self) -> str:
        return "5'-%s-3'" % self.strand

    def count_hydrogen_bonds(self):
        '''氢键计数'''
        n = 0
        for i in self.strand:
            if i in 'AT':
                n += 2
            else:
                n += 3
        return n


class RestrictionEnzyme():
    '''限制酶'''

    def __init__(self, rank, cut):
        self.type = 'enzyme'
        self.rank = rank
        self.cut = cut


class PCRMachine():
    def __init__(self):
        self.enzymelist = []
        self.dnadict = defaultdict(int)
        self.primerlist = []

    def __str__(self):
        return f'PCR仪\nDNA\n{pformat(dict(self.dnadict))}\n共有{len(self.dnadict)}种DNA\n酶{self.enzymelist}\n引物{self.primerlist}'

    def restart(self):
        self.enzymelist.clear()
        self.dnadict.clear()
        self.primerlist.clear()

    def add(self, *args):
        for i in args:
            if i.type == 'enzyme':
                self.enzymelist.append(i)
            elif i.type == 'DNA':
                self.dnadict[i] += 1
            elif i.type == 'primer':
                self.primerlist.append(i)

    def denature(self):
        '''变性'''
        origin_dnadict = self.dnadict.copy()
        self.dnadict.clear()
        for i in origin_dnadict:
            for j in i.dna_denature():
                self.dnadict[j] += origin_dnadict[i]

    def anneal(self):
        iterdict = self.dnadict.copy()
        self.dnadict.clear()
        for i in iterdict:
            for j in self.primerlist:
                try:
                    new = j.dna_pairing(i)
                    new.been_anneal=True
                    if new.template_strand == new.sense_strand[::-1]:
                        new.been_anneal=False
                    self.dnadict[new] += iterdict[i]
                except myerror.PairingError:
                    pass

    def extension(self):
        '''延伸'''
        iterdict = self.dnadict.copy()
        self.dnadict.clear()
        for i in iterdict:
            if i.been_anneal:
                a = copy.copy(i)
                a.sense_strand = pairing(i.template_strand)[
                    :i.start_sen+len(i.sense_strand)]
                a.start_sen = 0
                a.been_anneal = False
            self.dnadict[a] += iterdict[i]

    def check(self):
        for i in self.dnadict:
            print(i,
                  self.dnadict[i],
                  i.template_strand,
                  i.sense_strand,
                  i.start_sen,
                  i.is_single)
