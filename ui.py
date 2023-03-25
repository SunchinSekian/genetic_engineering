import streamlit as st
from main import *
import pandas as pd
st.set_page_config(page_title='你的基因工程实验室', page_icon=None,
                   layout="centered", initial_sidebar_state="auto", menu_items={'About': 'https://github.com/SunchinSekian/genetic_engineering'})

st.title('你的第一台虚拟PCR机')

pcrmachine = PCRMachine()
pcrmachine.name = st.text_input('给你的PCR仪一个好听的名字吧')
st.subheader(f'{pcrmachine.name}')
raw_dna = st.text_input('输入原始DNA单链,将自动生成双链',value='AGCTAA')
dna = DNA.init_with_single(raw_dna)
st.write(str(dna))
pcrmachine.add(dna)
primer1 = Primer(st.text_input('输入第一个引物',value='TAGC'))
primer2 = Primer(st.text_input('输入第二个引物',value='GCTA'))
pcrmachine.add(primer1, primer2)
times = st.slider('选择进行PCR循环轮数', min_value=1, max_value=45)

for i in range(times):
    pcrmachine.denature()
    pcrmachine.anneal()
    pcrmachine.extension()
showdict={}
for i in pcrmachine.dnadict:
    showdict[str(i)]=pcrmachine.dnadict[i]
st.json(showdict)