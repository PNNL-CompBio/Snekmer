import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt



print("script has begun")
print("reading data")

correct_prediction = 0
false_prediction = 0

#files = glob.glob("./snek_apply/output/learn/kmer-summary-UP000325827_1631477.csv")
files = glob.glob("./snek_apply/output/learn/kmer-summ*")

TF_1 = []
scores = []
TF_10 = []


end_ranks = list()
for file in files:
    print(file)
    data =pd.read_csv(str(file), index_col="__index_level_0__", header=0)
    #print(data)
    for row in data.iterrows():
        index = row[0]
        value = row[1]
        # print("index is", index)
        # print("value is", value)
        # print("value split is", str(value).split(",")[0])
        # print(str(index))
        if "known" in str(index):
            for i,val in enumerate(value):
                # print(str(item).split(",")[0])
                if i == 1000:
                    end_ranks.append(1000)
                    TF_1.append(0)
                    TF_10.append(0)
                    if (float((str(val).split(",")[1]))) > 0:
                        scores.append(0)
                    else:
                        scores.append(float((str(val).split(",")[1])))
                    break
                elif str(val).split(",")[0] in str(index):
                    # print("here1 :", (str(val).split(",")[1]))
                    # print("here2 :", (str(val).split(",")[1])[1:-1])
                    scores.append(float((str(val).split(",")[1])))
                    end_ranks.append(i)
                    if i == 1:
                        TF_1.append(1)
                    elif i != 1:
                        TF_1.append(0)
                    if i <= 10:
                        TF_10.append(1)
                    elif i > 10:
                        TF_10.append(0)
                    break
                    # print("truth: ", row)
    

print('Saving results to CSV.')
df = pd.DataFrame(data={"TF_1": TF_1, "TF_10": TF_10, "Scores":scores})
df.to_csv("./RUO_data.csv", sep=',',index=False)
# print('CSV Saved.')
# print('Now run pROC in Rstudio.')




print(end_ranks)
count_1 = end_ranks.count(0)
count_2 = end_ranks.count(1)
count_3 = end_ranks.count(2)
count_4 = end_ranks.count(3)
count_5 = end_ranks.count(4)
count_6 = end_ranks.count(5)
count_7 = end_ranks.count(6)
count_8 = end_ranks.count(7)
count_9 = end_ranks.count(8)
count_10 = end_ranks.count(9)

total = len(end_ranks)

print("total ", total)
print("Number of correct predictions in top1 ", count_1)
print("Percent of correct predictions in top1 ", count_1/total)

print("Number of correct predictions in top5 ", (count_1+count_2+count_3+count_4+count_5))
print("Percent of correct predictions in top5 ", (count_1+count_2+count_3+count_4+count_5)/total)


print("Number of correct predictions in top10 ", (count_1+count_2+count_3+count_4+count_5+count_6+count_7+count_8+count_9+count_10))
print("Percent of correct predictions in top10 ", (count_1+count_2+count_3+count_4+count_5+count_6+count_7+count_8+count_9+count_10)/total)


plt.hist(end_ranks, bins=7740)
plt.show()

plt.hist(end_ranks, bins=3450)
plt.show() 

plt.hist(end_ranks, bins=1725)
plt.show() 

plt.hist(end_ranks, bins=450)
plt.show() 

plt.hist(end_ranks, bins=12)
plt.show() 





######

# print("Starting ROC PART")






# print("\n\n\n more data for me.. ")
# correct_prediction = 0
# false_prediction = 0

# for row in data.iterrows():
#     index = row[0]
#     value = str(row[1][0]) +str(row[1][1])+str(row[1][2])+str(row[1][3])+str(row[1][4])
#     if "known" in str(value):
#         if str(value) in str(index):
#             correct_prediction += 1
#         else:
#             false_prediction += 1


# total = correct_prediction + false_prediction


# print("Accuracy for Top1-5 Prediction")
# print("Total is: ", total)
# print("Correct is: ", correct_prediction)
# print("Incorrect is: ", false_prediction)
# print("Percentage correct is: ", correct_prediction/total)
# print("Percentage incorrect is: ", false_prediction/total)
