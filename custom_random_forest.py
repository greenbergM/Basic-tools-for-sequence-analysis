from concurrent.futures import ThreadPoolExecutor
import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier


class RandomForestClassifierCustom(BaseEstimator):
    """Random Forest classifier with custom implementation."""

    def __init__(
        self, n_estimators: int = 10, max_depth: int = None, max_features: int = None, random_state: int = None
    ) -> None:
        """
        Initialize the RandomForestClassifierCustom instance.

        :param n_estimators: The number of trees in the forest (int, default=10).
        :param max_depth: The maximum depth of the tree (int, default=None).
        :param max_features: The number of features to consider when looking for the best split (int, default=None).
        :param random_state: Controls the randomness of the estimator (int, default=None).
        """
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def _fit_tree(self, i: int, random_state: int, max_features: int, max_depth: int, X: np.ndarray, y: np.ndarray) -> tuple:
        """Fit a single decision tree."""
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

    def fit(self, X: np.ndarray, y: np.ndarray, n_jobs: int = 1) -> "RandomForestClassifierCustom":
        """
        Fit the RandomForestClassifierCustom model.

        :param X: The input data array.
        :param y: The target labels.
        :param n_jobs: Number of jobs to run in parallel (int, default=1).

        :return: The fitted model instance.
        """
        self.classes_ = sorted(np.unique(y))

        args = [
            (i, self.random_state, self.max_features, self.max_depth, X, y)
            for i in range(self.n_estimators)
        ]

        with ThreadPoolExecutor(n_jobs) as pool:
            results = pool.map(self._fit_tree, *zip(*args))

        self.trees, self.feat_ids_by_tree = zip(*results)

        return self

    def _predict_proba(self, args: tuple) -> np.ndarray:
        """Predict probabilities for a single tree."""
        tree, feat_idx, X = args
        return tree.predict_proba(X[:, feat_idx])

    def predict_proba(self, X: np.ndarray, n_jobs: int = 1) -> np.ndarray:
        """
        Predict class probabilities for input samples.

        :param X: The input data array.
        :param n_jobs: Number of jobs to run in parallel (int, default=1).

        :return: Class probabilities of the input samples.
        """
        args = []
        for tree, feat_idx in zip(self.trees, self.feat_ids_by_tree):
            args.append((tree, feat_idx, X))

        with ThreadPoolExecutor(n_jobs) as pool:
            probas = list(pool.map(self._predict_proba, args))

        return np.mean(probas, axis=0)

    def predict(self, X: np.ndarray, n_jobs: int = 1) -> np.ndarray:
        """
        Predict the class labels for input samples.

        :param X: The input data array.
        :param n_jobs: Number of jobs to run in parallel (int, default=1).

        :return: Predicted class labels for the input samples.
        """
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)

        return predictions
