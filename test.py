from main import *
def check():
    for i in b.dnadict:
        print(i)
        print()


b=PCRMachine()
#b.add(a)


c=DNA.init_with_single('AGCTAA')

d=Primer('TAGC')
e=Primer('GCTA')
#print(d.dna_pairing(c))
b.add(c)
b.add(e)
b.add(d)
for i in range(2):
    b.denature()
    print(1)
    check()
    b.anneal()
    print(2)
    check()
    b.extension()
    print(3)
    check()
print(b)