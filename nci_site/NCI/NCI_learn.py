import os
import pdb
import pickle
import random

import sklearn
from sklearn import svm
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.datasets import load_svmlight_file
import pdb

base_p = "/Users/george/git/nci/nci_site/NCI/915/"
abstracts_path = os.path.join(base_p, "abstract/Cleaned/")
keywords_path = os.path.join(base_p, "keywords/")
titles_path = os.path.join(base_p, "title/Cleaned/")
outpath = "/Users/george/git/nci/nci_site/NCI/output/915"

class Ensemble:

    def __init__(self, abstract_m, titles_m, keywords_m):
        self.abstract_m = abstract_m
        self.title_m = titles_m
        self.keyword_m = keywords_m


    def predict(self, X_abstract, X_title, X_keywords):
        #return max([m.predict(X)])
        #pdb.set_trace()
        ab_pred = self.abstract_m.predict(X_abstract)[0]
        ti_pred = self.title_m.predict(X_title)[0]
        kw_pred = self.keyword_m.predict(X_keywords)[0]
        print "\n--ab: {0}; ti: {1}; keywords: {2}\n".format(ab_pred, ti_pred, kw_pred)
        return sum([ab_pred, ti_pred, kw_pred]) > 0


def train_NCI_models():
    committee_size = 11

    lbl_d = pickle.load(open(os.path.join(base_p, "labels.pickle")))

    #abstracts, titles, keywords, lbls = _load_data(lbl_d)
    ids_to_abstracts, ids_to_titles, ids_to_keywords, abstracts_vectorizer, titles_vectorizer, kw_vectorizer =\
                        _load_data(lbl_d)

    # now fit models
    committee = []
    pos_ids, neg_ids = [], []
    for x_id, lbl in lbl_d.items():
        if lbl == 1:
            pos_ids.append(x_id)
        else:
            neg_ids.append(x_id)

    n_pos = len(pos_ids)
    pos_ids = [str(id_) for id_ in pos_ids]
    for i in xrange(committee_size):
        # sample negative id's for this classifier.
        # note that this is not technically bootstrapping,
        # since we sample without replacement, but eh.
        neg_ids_i = [str(id_) for id_ in random.sample(neg_ids, n_pos)]
        lbls_i = [-1]*n_pos + [1]*n_pos

        abstracts_i = [ids_to_abstracts[id_] for id_ in neg_ids_i] + [ids_to_abstracts[id_] for id_ in pos_ids]
        abstracts_X = abstracts_vectorizer.transform(abstracts_i)
        titles_i = [ids_to_titles[id_] for id_ in neg_ids_i] + [ids_to_titles[id_] for id_ in pos_ids]
        titles_X = titles_vectorizer.transform(titles_i)
        keywords_i = [ids_to_keywords[id_] for id_ in neg_ids_i] + [ids_to_keywords[id_] for id_ in pos_ids]
        keywords_X = kw_vectorizer.transform(keywords_i)
    
        abstracts_svm = _get_svm()
        print "training abstracts SVM..."
        abstracts_svm.fit(abstracts_X, lbls_i)

        titles_svm = _get_svm()
        print "training titles SVM..."
        titles_svm.fit(titles_X, lbls_i)

        keywords_svm = _get_svm()
        print "training keyword SVM..."
        keywords_svm.fit(keywords_X, lbls_i)

        ensemble_i = Ensemble(abstracts_svm, titles_svm, keywords_svm)
        committee.append(ensemble_i)

        #pdb.set_trace()
        #a_vec = abstracts_vectorizer.transform([abstracts_i[0]])
    print "dumping..."
    with open(os.path.join(outpath, "committee.pickle"), 'w') as f_out:
        pickle.dump(committee, f_out)


def _get_svm():
    # in case we want to make settings
    return svm.SVC(kernel="linear")

def _load_data(ids_to_labels):
    ids_to_abstracts, ids_to_titles, ids_to_keywords = {}, {}, {}
    for id_ in ids_to_labels.keys():
        # note that we pickle the ids as numbers, but need them
        # as strings here.
        str_id = str(id_)
        ids_to_abstracts[str_id] = open(os.path.join(abstracts_path, str_id)).readline()
        ids_to_titles[str_id] = open(os.path.join(titles_path, str_id)).readline()
        ids_to_keywords[str_id] = " ".join(open(os.path.join(keywords_path, str_id)).readline().split(","))

    # get feature space mapping functions
    print "encoding abstracts..."
    abstracts_vectorizer = _vectorize(ids_to_abstracts.values())
    print "ok. encoding titles..."
    titles_vectorizer = _vectorize(ids_to_titles.values())
    print "ok. encoding keywords..."
    #pdb.set_trace()
    kw_vectorizer = _vectorize(ids_to_keywords.values())
    print "all done encoding."

    # and pickle/write them out
    vectorizers = [abstracts_vectorizer, titles_vectorizer, kw_vectorizer]
    v_strs = ["abstracts", "titles", "keywords"]
    for field, v in zip(v_strs, vectorizers):
        with open(os.path.join(outpath, "{0}_vectorizer.pickle".format(field)), 'w') as v_out:
            pickle.dump(v, v_out)

    #return abstracts, titles, keywords, labels
    return ids_to_abstracts, ids_to_titles, ids_to_keywords, abstracts_vectorizer, titles_vectorizer, kw_vectorizer


def _vectorize(texts):
    ''' texts is e.g., a list of abstracts '''
    vectorizer = CountVectorizer(binary=True, max_features=50000, min_df=3, stop_words="english")
    vectorizer.fit(texts)
    return vectorizer

'''
.predict(a_vec.transform([an_abstract]), title_vec.transform([a_title]), title_vec.transform([""]))

'''



