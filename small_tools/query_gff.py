#!/usr/bin/env python
# coding: utf-8

# In[67]:


#!/usr/bin/env python3
# coding: utf-8

"""
2020/08/26
author:guisen chen
email:thecgs001@foxmail.com
usage:./query_gff.py id.txt file.gff
"""

import sys

def query_gff(id_file, gff_file):
    id_file = open(id_file,'r')
    gff_file = open(gff_file,'r')
    gff = gff_file.read()
    query_file = open('query.gff','a')
    
    list = gff.split('\n\n') #以空行为分割
    query = []
    for id in id_file:
        for line in list:
            if str(line).find(id.strip()) > 0:
                query.append(line)
    query.sort() #按照染色体顺序排序
    for line in query:
        query_file.write(f'{str(line)}\n\n')
    query_file.close()

#执行
query_gff(sys.argv[1],sys.argv[2])
print('--------completed--------')


# In[ ]:




