import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt
import requests
import random
f = open("./output/drug_result_v8_final.txt", mode='r', encoding='utf-8')
drug_lines = f.readlines()
pred_list = list()
true_list = list()
true_list2 = list()
all_list2 = list()
f.close()

drug_sets = set()
all_list = list()
url = 'https://api.drugbankplus.com/v1/indications/drugs'
synonyms = ['Malignant Neoplasm of Stomach', 'Adenocarcinoma of the Stomach','Metastatic Adenocarcinoma of the Stomach','Stage 4 gastrointestinal adenocarcinoma', 'Stomach cancer', 'Stomach neoplasm', 'Stomach carcinoma', 'Gastric cancer', 'Gastric neoplasm', 'Gastric carcinoma', 'Cancer of stomach', 'Cancer of gastric']
paramdict = dict()
for disease in synonyms:
    paramdict['q'] = disease
    response = requests.get(url, params=paramdict, headers={'Authorization': '45eae790a64e50c7905d7bc878eb90ed'})
    print(response.status_code)
    res_list =response.json()
    for drug in res_list:
        drug_sets.add(drug["drugbank_id"])
for each_drug in drug_lines:
    items = each_drug.split('\t')
    drug_id = items[0]
    drug_score = items[1].split('\n')[0]
    all_list2.append(float(drug_score)*(10**3))

    if drug_id in drug_sets:
        print(drug_id)
        pred_list.append(float(drug_score)*(10**3))
        true_list.append(1)
        true_list2.append(1)
    else:
        true_list2.append(0)
        all_list.append(float(drug_score)*(10**3))


true_array = np.array(true_list2, dtype='f')
pred_array = np.array(all_list2, dtype='f')
print(roc_auc_score(true_array, pred_array))
#
# sum = 0
# for j in range(10000):
#     temp_true = true_list[:]
#     temp_pred = pred_list[:]
#     fake_list = random.sample(all_list, len(temp_pred))
#     for i in range(len(fake_list)):
#         temp_true.append(0)
#         temp_pred.append(fake_list[i])
# # for i in range(10):
# #     rand = random.randrange(0, len(true_list)-1)
# #     rand_list[rand] = 1
# # y_random = np.array(rand_list, dtype='f')
#
#     y_true = np.array(temp_true, dtype='f')
#     y_scores = np.array(temp_pred, dtype='f')
# # print(roc_auc_score(y_random, y_scores))
#     sum += roc_auc_score(y_true, y_scores)
#     # print(roc_auc_score(y_true , y_scores))
fpr1, tpr1, thresholds1 = roc_curve(true_array, pred_array)
# # fpr2, tpr2, thresholds1 = roc_curve(y_random, y_scores)
# print("average is ", sum/10000)

plt.plot(fpr1, tpr1, 'o-', ms=2, label="my model")
# # plt.plot(fpr2, tpr2, 'o-', ms=2, label="random model")
plt.legend()
plt.plot([0, 1], [0, 1], 'k--', label="random guess")
plt.xlabel('FPR(False Positive Rate)')
plt.ylabel('TPR(True Positive Rate)')
plt.title('ROC curve')
plt.savefig("./output/ROC_curve_final.png")
