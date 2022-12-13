import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import itertools
from math import ceil

print("script has begun")
print("reading data")



files = glob.glob("./snek_apply/output/learn/Seq-Annotation-Scores-*")

def chunk_into_n(lst, n):
  size = ceil(len(lst) / n)
  return list(
    map(lambda x: lst[x * size:x * size + size],
    list(range(n)))
  )

# end_ranks = list()
first_time = True
for i,file in enumerate(files):
    data =pd.read_csv(str(file), index_col="__index_level_0__", header=0)

    list_1 = data.index
    list_2 = data.columns

    comparisons = [int(a in b) for (a, b) in itertools.product(list_2, list_1)]

    qq = chunk_into_n(comparisons,len(list_2))

    new = pd.DataFrame(qq).transpose()
    new.columns =list_2
    new.index=list_1
    print(new)



    TF_out = ("./ruo_output/RUO_AUC_TF_data_" + str(file[48:-4]) + ".csv")
    print('Saving results to CSV.')
    new.to_csv(TF_out, sep=',',index=True)

    print(TF_out, ' Saved.')
    print("Files complete: ", i)


print('Now run pROC in Rstudio.')



