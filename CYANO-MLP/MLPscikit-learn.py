import numpy as np
import sys
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn import metrics
from sklearn.pipeline import make_pipeline
from sklearn import preprocessing
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import permutation_test_score
from itertools import product

y = []
X = []
y_names = []
with open("taxon_eight_clades_new_alpha.txt", "r") as taxon_class:
    for line in taxon_class:
        class_name = line.split()[0]
        y_names.append(line.split()[1])
        y.append(class_name)

with open("eight_model_LOOCV_scores_alpha.txt", "r") as lcv_scores:
    for line in lcv_scores:
        scores = line.split()[1:-1]
        scores = [float(x) for x in scores]
        X.append(scores)

X = np.array(X)
y = np.array(y)

###############################################################################
#Testing for different network architecture.
#Uncomment lines 33-50 and comment lines 51-end to run tests.
#test_layer_list = []
#test_layer_list.extend(list(product([8,9,10,11,12,13,14,15,16], repeat=1)))
#test_layer_list.extend(list(product([8,9,10,11,12,13,14,15,16], repeat=2)))
#test_layer_list.extend(list(product([8,9,10,11,12,13,14,15,16], repeat=3)))
#test_layer_list.extend(list(product([8,9,10,11,12,13,14,15,16], repeat=4)))
#[coef.shape for coef in clf.coefs_]
#for layer in test_layer_list:
#    clf = make_pipline(preprocessing.StandardScaler(),MLPClassifier(solver="lbfgs", alpha=.01, max_iter=4000,
#                        hidden_layer_sizes = layer,random_state=1))
    #X = np.array(X)
    #y = np.array(y)
#    loo = LeaveOneOut()
#    scores = cross_val_score(clf, X, y, cv = loo, n_jobs = -1)
#    predicted = cross_val_predict(clf, X, y, cv = loo, n_jobs = -1)
#    print("layout:{}\tScore:{}".format(layer, np.average(scores)))
###############################################################################

clf = make_pipeline(preprocessing.StandardScaler(),MLPClassifier(solver="lbfgs", alpha=.01, max_iter=4000,
                    hidden_layer_sizes = (13,),random_state=1, tol = 1e-4))

###############################################################################
#LOOCV code. To run LOOCV code uncomment lines 55-63 
#loo = LeaveOneOut()
#scores = cross_val_score(clf, X, y, cv = loo, n_jobs = -1)
#print("Score:{}".format(np.average(scores)))
#predicted = cross_val_predict(clf, X, y, cv = loo, n_jobs = -1, method = "predict_proba")
#
#print("Clade\tPredicted\tSpecies")
#for i, pred in enumerate(predicted):
#    pred_list = pred.tolist()
#    print("{}\t{}\t{}".format(y[i],"\t".join([str(y) for y in pred_list]), y_names[i]))
###############################################################################

###############################################################################
#permutation label swaping accuracy code. To run uncomment lines 66-87
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
###############################################################################

clf.fit(X, y)

#Code to test acc on full dataset
#predicted = clf.predict(X)
#print(metrics.accuracy_score(y, predicted))
##
#preds = clf.predict_proba(X)
#print("A\tB1\tB23\tC1\tC3\tE\tF\tG\tpredicted\tClade\tID")
#for i, pred in enumerate(preds):
#    pred_list = pred.tolist()
#    print("{}\t{}\t{}\t{}".format("\t".join([str(y) for y in pred_list]),clf.predict(X[i].reshape(1, -1)), y[i], y_names[i]))


################################################################################
#Code to classify list of score vectors contained in the first command line argument
id_info = []
with open(sys.argv[1], "r") as eudicots:
    eudi_classify = []
    for line in eudicots:
        line_split = line.split()
        scores = line_split[1:-1]
        scores = [float(x) for x in scores]
        id_info.append([line_split[0], line_split[-1]])
        eudi_classify.append(scores)

print("A\tB1\tB23\tC1\tC3\tE\tF\tG\tpredicted\tClade\tID")
eudi_classify = np.array(eudi_classify)
preds = clf.predict_proba(eudi_classify)
for i, pred in enumerate(preds):
    pred_list = pred.tolist()
    print("{}\t{}\t{}".format("\t".join([str(y) for y in pred_list]),clf.predict(eudi_classify[i].reshape(1, -1)), "\t".join(id_info[i])))
################################################################################
