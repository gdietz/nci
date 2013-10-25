import pickle
import os
import pdb
import sklearn
import NCI_learn

''' change me '''
path_to_abstract_vectorizer = "output/915/abstracts_vectorizer.pickle"
path_to_keyword_vectorizer = "output/915/keywords_vectorizer.pickle"
path_to_title_vectorizer = "output/915/titles_vectorizer.pickle"

path_to_model = "output/915/committee.pickle"

class NCI_predictor:
    def __init__(self, nci_dir_path="/Users/george/git/nci/nci_site/NCI"):
        print "loading vectorizers..."
        #import pdb; pdb.set_trace()
        self.abstract_vectorizer = pickle.load(open(os.path.join(nci_dir_path, path_to_abstract_vectorizer)))
        self.title_vectorizer    = pickle.load(open(os.path.join(nci_dir_path, path_to_title_vectorizer)))
        self.keyword_vectorizer  = pickle.load(open(os.path.join(nci_dir_path, path_to_keyword_vectorizer)))

        print "ok. loading model.."
        self.committee = pickle.load(open(os.path.join(nci_dir_path, path_to_model)))
        print "ok. all set."

    def predict(self, title, abstract, keywords):
        ti_vec = self.title_vectorizer.transform([title])
        ab_vec = self.abstract_vectorizer.transform([abstract])
        kw_vec = self.keyword_vectorizer.transform([keywords])

        
        votes = [ensemble.predict(ab_vec, ti_vec, kw_vec) for ensemble in self.committee]
        print "ensemble votes: {0}\n".format(votes)
        return votes.count(True) >= float(len(votes))/2.0