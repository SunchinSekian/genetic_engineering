from main import *


def check():
    for i in b.dnadict:
        print(i)
        print()


b = PCRMachine()
# b.add(a)


c = DNA.init_with_single('AGCTAA')

d = Primer('TAGC')
e = Primer('GCTA')
# print(d.dna_pairing(c))
b.add(c)
b.add(e)
b.add(d)
for i in range(30):
    b.denature()
    b.anneal()
    b.extension()
print(b)
