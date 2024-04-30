from concurrent.futures import ThreadPoolExecutor
import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(
        self, n_estimators=10, max_depth=None, max_features=None, random_state=None
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def _fit_tree(self, i, random_state, max_features, max_depth, X, y):
        np.random.seed(random_state + i)
        feat_idx = np.random.choice(X.shape[1], max_features, replace=False)

        bootstrap_idx = np.random.choice(X.shape[0], X.shape[0], replace=True)
        bootstrap_X = X[bootstrap_idx]
        bootstrap_y = y[bootstrap_idx]

        tree = DecisionTreeClassifier(
            max_depth=max_depth, random_state=(random_state + i)
        )
        tree.fit(bootstrap_X[:, feat_idx], bootstrap_y)

        return tree, feat_idx

    def fit(self, X, y, n_jobs=1):
        self.classes_ = sorted(np.unique(y))

        args = [
            (i, self.random_state, self.max_features, self.max_depth, X, y)
            for i in range(self.n_estimators)
        ]

        with ThreadPoolExecutor(n_jobs) as pool:
            results = pool.map(self._fit_tree, *zip(*args))

        self.trees, self.feat_ids_by_tree = zip(*results)

        return self

    def _predict_proba(self, args):
        tree, feat_idx, X = args
        return tree.predict_proba(X[:, feat_idx])

    def predict_proba(self, X, n_jobs=1):
        args = []
        for tree, feat_idx in zip(self.trees, self.feat_ids_by_tree):
            args.append((tree, feat_idx, X))

        with ThreadPoolExecutor(n_jobs) as pool:
            probas = list(pool.map(self._predict_proba, args))

        return np.mean(probas, axis=0)

    def predict(self, X, n_jobs=1):
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)

        return predictions
