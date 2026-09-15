"""Simulated active learning benchmark on pre-docked scores."""

import random
import numpy as np

import src.featurization as featurization
import src.model as model

SCORES_CSV = "data/results/3ert/scores.csv"
LIBRARY_CSV = "data/library/fda.csv"

SEED_SIZE = 300
BATCH_SIZE = 100
ROUNDS = 10
HIT_THRESHOLD = -9.0
N_ENSEMBLE = 5


def select_random(mean_pred, std_pred, n, rng):
    return rng.sample(range(len(mean_pred)), n)


def select_exploit(mean_pred, std_pred, n, rng):
    ranked = np.argsort(mean_pred)
    return ranked[:n]


def select_uncertainty(mean_pred, std_pred, n, rng, k=1.0):
    score = mean_pred - k * std_pred
    uncertain = np.argsort(score)
    return uncertain[:n]


def train_ensemble(X, y, rows, n_models=N_ENSEMBLE):
    """Train n_models XGBoost regressors on the labeled rows."""
    return [
        model.train_model(X[rows], y[rows], random_state=s)
        for s in range(n_models)
    ]


def predict_ensemble(models, X, rows):
    """Return (mean, std) of ensemble predictions for the given rows."""
    preds = np.array([m.predict(X[rows]) for m in models])
    return preds.mean(axis=0), preds.std(axis=0)


def count_hits(y, rows, threshold=HIT_THRESHOLD):
    return int((y[rows] < threshold).sum())


def run(strategy_fn, seed=42):
    X, y, ids = featurization.build_dataset(SCORES_CSV, LIBRARY_CSV)

    rng = random.Random(seed)
    all_rows = list(range(len(ids)))
    rng.shuffle(all_rows)
    labeled = all_rows[:SEED_SIZE]
    unlabeled = all_rows[SEED_SIZE:]

    history = [(len(labeled), count_hits(y, labeled))]

    for _ in range(ROUNDS):
        if not unlabeled:
            break

        models = train_ensemble(X, y, labeled)
        mean_pred, std_pred = predict_ensemble(models, X, unlabeled)

        n = min(BATCH_SIZE, len(unlabeled))
        picks = strategy_fn(mean_pred, std_pred, n, rng)

        chosen = [unlabeled[i] for i in picks]
        labeled.extend(chosen)
        unlabeled = [r for i, r in enumerate(unlabeled) if i not in set(picks)]

        history.append((len(labeled), count_hits(y, labeled)))

    return history


if __name__ == "__main__":
    for name, fn in [("random", select_random),
                     ("exploit", select_exploit),
                     ("uncertainty", select_uncertainty)]:
        runs = [run(fn, seed=s) for s in [1, 2, 3, 4, 5]]
        # runs[i][r] is (n_docked, hits) for seed i, round r
        for r in range(len(runs[0])):
            n = runs[0][r][0]
            hits = [runs[i][r][1] for i in range(len(runs))]
            print(name, n, round(np.mean(hits), 1), round(np.std(hits), 1))