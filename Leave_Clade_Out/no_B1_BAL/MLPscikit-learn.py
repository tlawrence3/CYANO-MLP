import numpy as np
import sys
import matplotlib.pyplot as plt
import scikitplot as skplt
from sklearn import preprocessing
from sklearn.pipeline import make_pipeline
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn import metrics
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import permutation_test_score
from itertools import product

y = []
X = []
y_names = []
with open("taxalist_resample_noB1.txt", "r") as taxon_class:
    for line in taxon_class:
        class_name = line.split()[0]
        y_names.append(line.split()[1])
        y.append(class_name)

with open("loocv_resample_noB1.txt", "r") as lcv_scores:
    for line in lcv_scores:
        scores = line.split()[1:-1]
        scores = [float(x) for x in scores]
        X.append(scores)

X = np.array(X)
X_mean = np.mean(X, axis = 0, dtype=np.float64)
X_std = np.std(X, axis = 0, dtype=np.float64)
#X = (X - np.mean(X, axis = 0, dtype=np.float64)) / np.std(X, axis = 0, dtype=np.float64)
y = np.array(y)

#Art testing
#test_layer_list = []
#test_layer_list.extend(list(product([8,9,10,11,12,13,14,15,16], repeat=1)))
#test_layer_list.extend(list(product([8,9,10,11,12,13,14,15,16], repeat=2)))
#test_layer_list.extend(list(product([8,9,10,11,12,13,14,15,16], repeat=3)))
#test_layer_list.extend(list(product([8,9,10,11,12,13,14,15,16], repeat=4)))
#[coef.shape for coef in clf.coefs_]
#for layer in test_layer_list:
#    clf = MLPClassifier(solver="lbfgs", alpha=.01, max_iter=2000,
#                        hidden_layer_sizes = layer,random_state=1)
    #X = np.array(X)
    #X_mean = np.mean(X, axis = 0, dtype=np.float64)
    #X_std = np.std(X, axis = 0, dtype=np.float64)
    #X = (X - np.mean(X, axis = 0, dtype=np.float64)) / np.std(X, axis = 0, dtype=np.float64)
    #y = np.array(y)
#    loo = LeaveOneOut()
#    scores = cross_val_score(clf, X, y, cv = loo, n_jobs = -1)
#    predicted = cross_val_predict(clf, X, y, cv = loo, n_jobs = -1)
#    print("layout:{}\tScore:{}".format(layer, np.average(scores)))
###Comment out for pipeline for CV
#clf = MLPClassifier(solver="lbfgs", alpha=.01, max_iter=4000,
#                    hidden_layer_sizes = (3, 10, 11), random_state=1, tol = 1e-4)


###CV setup
#(11, 14)
clf = make_pipeline(preprocessing.StandardScaler(),MLPClassifier(solver="lbfgs", alpha=.01, max_iter=4000,
                                                                 hidden_layer_sizes =(16, 12, 16, 15),random_state=1, tol = 1e-4))
## 10-fold CV
#skf = StratifiedKFold(n_splits = 10)
#scores = cross_val_score(clf, X, y, cv = skf, n_jobs = -1)
#print("Score:{}".format(np.average(scores)))
#print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std()))

#LOOCV code. Don't need to run everytime
loo = LeaveOneOut()
#scores = cross_val_score(clf, X, y, cv = loo, n_jobs = -1)
#print("Score:{}".format(np.average(scores)))
#predicted = cross_val_predict(clf, X, y, cv = loo, n_jobs = -1, method = "predict_proba")
#Make confusion matrix
predicted = cross_val_predict(clf, X, y, cv = loo, n_jobs = -1)
skplt.metrics.plot_confusion_matrix(y, predicted, normalize=True)
plt.savefig("no_B1_upsample_confusionMatrix.pdf")


#print(scores.mean())
#print("Clade\tPredicted\tSpecies")
#for i, pred in enumerate(predicted):
#    pred_list = pred.tolist()
#    print("{}\t{}\t{}".format(y[i],"\t".join([str(y) for y in pred_list]), y_names[i]))
#    print("{}\t{}\t{}".format(y[i],pred, y_names[i]))

###permutation label swaping accuracy
#n_classes = np.unique(y).size
#score, permutation_scores, pvalue = permutation_test_score(clf, X, y,
#                                                           scoring="accuracy",
#                                                           cv=loo,
#                                                           n_permutations=300,
#                                                           n_jobs=-1,
#                                                           random_state = 1)
#print("Classification score %s (pvalue : %s)" % (score, pvalue))
#sns.distplot(permutation_scores, bins=20, kde=False, rug=True)
#plt.hist(permutation_scores, 20, label='Permutation scores',
#         edgecolor='black')
#ylim = plt.ylim()
#plt.plot(2 * [score], ylim, '--g', linewidth=3,
#         label='Classification Score'
#         ' (pvalue %s)' % pvalue)
#plt.plot(2 * [1. / n_classes], ylim, '--k', linewidth=3)

#plt.ylim(ylim)
#plt.legend()
#plt.xlabel('Score')
#plt.savefig("loocv.acc.perm1000.pdf")
#plt.close()
######################

clf.fit(X, y)

##Code to test acc on full dataset
#predicted = clf.predict(X)
#print(metrics.accuracy_score(y, predicted))
##
#preds = clf.predict_proba(X)
#print("A\tB1\tB23\tC1\tC3\tE\tF\tG\tpredicted\tClade\tID")
#for i, pred in enumerate(preds):
#    pred_list = pred.tolist()
#    print("{}\t{}\t{}\t{}".format("\t".join([str(y) for y in pred_list]),clf.predict(X[i].reshape(1, -1)), y[i], y_names[i]))


####Code for classifying dataset of interest
id_info = []
with open(sys.argv[1], "r") as eudicots:
    eudi_classify = []
    for line in eudicots:
        line_split = line.split()
        scores = line_split[1:-1]
        scores = [float(x) for x in scores]
        id_info.append([line_split[0], line_split[-1]])
        eudi_classify.append(scores)

print("A\tB23\tC1\tC3\tE\tF\tG\tpredicted\tClade\tID")
eudi_classify = np.array(eudi_classify)
#eudi_classify = (eudi_classify - X_mean) / X_std
preds = clf.predict_proba(eudi_classify)
for i, pred in enumerate(preds):
    pred_list = pred.tolist()
    print("{}\t{}\t{}".format("\t".join([str(y) for y in pred_list]),clf.predict(eudi_classify[i].reshape(1, -1)), "\t".join(id_info[i])))
