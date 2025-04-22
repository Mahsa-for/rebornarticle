import pandas as pd
import numpy as np
from sklearn import datasets
from scipy.stats import ttest_ind

iris = datasets.load_iris()
iris = pd.DataFrame(data = np.c_[iris["data"], iris["target"]], columns = iris["feature_names"] + ["target"])

target_mapping = {0: "setosa", 1: "versicolor", 2: "virginica"}
iris["species"] = iris["target"].map(target_mapping)

iris = iris[(iris.species == "setosa") | (iris.species == "virginica")].loc[:, ["petal length (cm)", "species"]]

x = iris[iris.species == "setosa"].loc[:, "petal length (cm)"]
y = iris[iris.species == "virginica"].loc[:, "petal length (cm)"]

tt = ttest_ind(x, y, equal_var = False)

df = pd.DataFrame({ "t.statistic": [tt.statistic] , "p.value": [tt.pvalue] })
print(df)
