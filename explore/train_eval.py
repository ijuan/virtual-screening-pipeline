import src.split as split
import src.model as model
import src.featurization as featurization
import pandas as pd

scores_csv = "data/results/scores.csv"
library_csv = "data/library/fda.csv"

X, y, ids = featurization.build_dataset(scores_csv, library_csv)
row_of = dict(zip(ids, range(len(ids))))
lib = pd.read_csv(library_csv)
smiles_by_id = dict(zip(lib["chembl_id"], lib["smiles"]))

train_ids, test_ids = split.scaffold_split(ids, smiles_by_id)
train_rows = [row_of[i] for i in train_ids]
test_rows = [row_of[i] for i in test_ids]
X_train, y_train = X[train_rows], y[train_rows]
X_test, y_test = X[test_rows], y[test_rows]
clf = model.train_model(X_train, y_train)

mae, spearman = model.evaluate(clf, X_test, y_test)
base_mae, base_spearman = model.baseline_mean(y_train, y_test)

print(f"model    MAE {mae:.3f}  Spearman {spearman:.3f}")
print(f"baseline MAE {base_mae:.3f}  Spearman {base_spearman}")


