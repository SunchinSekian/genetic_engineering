from main import *

a=DNA.init_with_single('ATCG')

a=random_DNA(100)

b=PCRMachine()
#b.add(a)


c=DNA.init_with_single('ATCGCTA')
d=Primer('TAGC')
e=Primer('ATCG')
#print(d.dna_pairing(c))
b.add(c)
b.add(e)
b.add(d)
for i in range(10):
    b.denature()
    b.anneal()
    b.extension()
print(b)