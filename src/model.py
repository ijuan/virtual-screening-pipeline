"""Train and evaluate the docking-score prediction model."""

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error
from xgboost import XGBRegressor


def train_model(X_train, y_train, n_estimators=300, max_depth=6):
    """Train an XGBoost regressor on docking scores."""
    model = XGBRegressor(
        n_estimators=n_estimators,
        max_depth=max_depth,
        random_state=42,
    )
    model.fit(X_train, y_train)
    return model


def evaluate(model, X_test, y_test):
    """Returns (mae, spearman) for the model's predictions on the test set."""
    predictions = model.predict(X_test)
    mae = mean_absolute_error(y_test, predictions)
    spearman = spearmanr(y_test, predictions).correlation
    return (mae, spearman)


def baseline_mean(y_train, y_test):
    """Predict-the-mean baseline. Returns (mae, spearman)."""
    predictions = np.full(len(y_test), y_train.mean())
    mae = mean_absolute_error(y_test, predictions)
    return (mae, float("nan"))