import random

def pairing(strand1):
    '''配对'''
    pw={'A':'T','C':'G','T':'A','G':'C'}
    strand2=''
    for i in strand1:
        strand2+=pw[i]
    return strand2
def dict_append(item,dic):
    '''字典添加'''
    try:
        dic[item]+=1
    except:
        dic[item]=1
def random_DNA(n):
    '''生成随机DNA'''
    baselist=''
    m='ATCG'
    for i in range(n):
        baselist+=m[random.randint(0,3)]
    return DNA(baselist)

def cut(DNA,enzeme):    
    pass


class DNA:
    '''DNA类'''
    def __init__(self,template_strand,sense_strand='',start_sen=0,auto=True):
        self.type='DNA'
        self.template_strand=template_strand    #模板链
        self.sense_strand=sense_strand
        if auto:
            self.sense_strand=pairing(self.template_strand)     #编码链
        self.start_sen=start_sen
        if self.sense_strand=='':
            self.is_single=True
        else:
            self.is_single=False       
    def __str__(self):
        if self.sense_strand!='':
            a="5'-%s-3'\n3'-%s-5'"%(self.template_strand,' '*self.start_sen+self.sense_strand)
        else:
            a="5'-%s-3'"%(self.template_strand)
        return a
    def __repr__(self) -> str:
        if self.sense_strand!='':
            a="5'-%s-3'\n3'-%s-5'"%(self.template_strand,self.sense_strand)
        else:
            a="5'-%s-3'"%(self.template_strand)
        return a
    def __len__(self):
        return max(len(self.sense_strand)+self.start_sen,len(self.template_strand))        
    def count_hydrogen_bonds(self):
        '''氢键计数'''
        n=0
        for i in self.template_strand[self.start_sen:]:
            if i in 'AT' :
                n+=2
            else:
                n+=3
        return n
    
    def generate_RNA(self):
        '''转录形成RNA'''
        pw={'A':'T','C':'U','T':'A','G':'C'}
        strand2=''
        for i in self.template_strand:
            strand2+=pw[i]
        return RNA(strand2)
    
    def denature(self):
        '''变性'''
        return DNA(self.template_strand,auto=False),DNA(self.sense_strand,auto=False)
    def dna_pairing(self,another_dna):
        '''主动与其他DNA配对'''
        if self.is_single and another_dna.is_single is True:
            position=another_dna.template_strand.find(pairing(self.template_strand)[::-1])
        else:
            raise'只能暂时配对单链'
        if position==-1:
            raise'这两段不能配对'
        return DNA(another_dna.template_strand,self.template_strand[::-1],start_sen=position,auto=False)
    
class RNA:
    '''RNA类'''
    def __init__(self,strand):
        self.strand=strand[::-1]
    def __str__(self):
        return "5'-%s-3'"%self.strand
    def __repr__(self) -> str:
        return "5'-%s-3'"%self.strand
    def count_hydrogen_bonds(self): 
        '''氢键计数'''
        n=0
        for i in self.strand:
            if i in 'AT':
                n+=2
            else:
                n+=3
        return n
 
class RestrictionEnzyme():
    '''限制酶'''
    def __init__(self,rank,cut):
        self.type='enzyme'
        self.rank=rank
        self.cut=cut
        
class PCRMachine():
    def __init__(self):
        self.enzymelist=[]
        self.dnadict={}
        self.primerlist=[]
    def __str__(self):
        return 'PCR仪\nDNA%s\n共有%s种DNA\n酶%s\n引物%s'%(self.dnadict,len(self.dnadict),self.enzymelist,self.primerlist)
    def restart(self):
        self.enzymelist=[]
        self.dnalist=[]
    def add(self,*args):
        for i in args:
            if i.type=='enzyme':
                self.enzymelist.append(i)
            elif i.type=='DNA':
                dict_append(i,self.dnadict)
    def denature(self):
        '''变性'''
        origin_dnadict=self.dnadict
        self.dnadict={}
        for i in origin_dnadict:
            for j in i.denature():
                dict_append(j,self.dnadict)
    def anneal(self):
        '''退火'''
        pass


