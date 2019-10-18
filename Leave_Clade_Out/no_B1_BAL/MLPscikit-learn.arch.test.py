import numpy as np
from sklearn import preprocessing
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from itertools import product
from sklearn.pipeline import make_pipeline
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
#X_mean = np.mean(X, axis = 0, dtype=np.float64)
#X_std = np.std(X, axis = 0, dtype=np.float64)
#X = (X - np.mean(X, axis = 0, dtype=np.float64)) / np.std(X, axis = 0, dtype=np.float64)
y = np.array(y)

test_layer_list = []
test_layer_list.extend(list(product([8,9,10,11,12,13,14,15,16], repeat=1)))
test_layer_list.extend(list(product([8,9,10,11,12,13,14,15,16], repeat=2)))
test_layer_list.extend(list(product([8,9,10,11,12,13,14,15,16], repeat=3)))
test_layer_list.extend(list(product([8,9,10,11,12,13,14,15,16], repeat=4)))
#[coef.shape for coef in clf.coefs_]
for layer in test_layer_list:
    clf = make_pipeline(preprocessing.StandardScaler(),MLPClassifier(solver="lbfgs", alpha=.01, max_iter=4000,
                                                                     hidden_layer_sizes = layer,random_state=1,
                                                                     tol = 1e-4))
    #clf = MLPClassifier(solver="lbfgs", alpha=.01, max_iter=2000,
     #                   hidden_layer_sizes = layer,random_state=1)
    loo = LeaveOneOut()
    scores = cross_val_score(clf, X, y, cv = loo, n_jobs = -1)
    predicted = cross_val_predict(clf, X, y, cv = loo, n_jobs = -1)
    print("layout:{}\tScore:{}".format(layer, scores.mean()))
    print("predicted:{}".format(predicted))
#clf = MLPClassifier(solver="lbfgs", alpha=.01, max_iter=2000,
#                    hidden_layer_sizes = (11,15,15,11),random_state=1)
#loo = LeaveOneOut()
#scores = cross_val_score(clf, X, y, cv = loo, n_jobs = -1)
#predicted = cross_val_predict(clf, X, y, cv = loo, n_jobs = -1)
##print(scores)
#print("Score:{}".format(np.average(scores)))
#print("Clade\tPredicted")
#for i, pred in enumerate(predicted):
#    print("{}\t{}".format(y[i],pred))
#print(y)
#clf.fit(X, y)
##print(scores)
##for i, score in enumerate(scores):
##    if (score == 0):
##        print("{} {}".format(y[i], y_names[i]))
#predicted = clf.predict(X)
##print(metrics.accuracy_score(y, predicted))
##
##preds = clf.predict_proba(X)
##print("A\tB1\tB23\tC1\tC3\tE\tF\tG\tpredicted\tClade\tID")
##for i, pred in enumerate(preds):
##    pred_list = pred.tolist()
##    print("{}\t{}\t{}\t{}".format("\t".join([str(y) for y in pred_list]),clf.predict(X[i].reshape(1, -1)), y[i], y_names[i]))
##
#id_info = []
#with open("CPH.test.txt", "r") as eudicots:
#    eudi_classify = []
#    for line in eudicots:
#        line_split = line.split()
#        scores = line_split[1:-1]
#        scores = [float(x) for x in scores]
#        id_info.append([line_split[0], line_split[-1]])
#        eudi_classify.append(scores)
#
#print("A\tB1\tB23\tC1\tC3\tE\tF\tG\tpredicted\tClade\tID")
#eudi_classify = np.array(eudi_classify)
#eudi_classify = (eudi_classify - X_mean) / X_std
#preds = clf.predict_proba(eudi_classify)
#for i, pred in enumerate(preds):
#    pred_list = pred.tolist()
#    print("{}\t{}\t{}".format("\t".join([str(y) for y in pred_list]),clf.predict(eudi_classify[i].reshape(1, -1)), "\t".join(id_info[i])))
#
##for i, x in enumerate(eudi_cla:ssify):
##    pred = clf.predict_proba(x)
##    pred_list = pred.tolist()
##    print("{}\t{}\t{}".format("\t".join([str(y) for y in pred_list]),clf.predict(x), "\t".join(id_info[i])))
