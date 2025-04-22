from dtreg.load_datatype import load_datatype
from dtreg.to_jsonld import to_jsonld

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

dt1 = load_datatype("https://doi.org/21.T11969/feeb33ad3e4440682a4d") # Data analysis
dt2 = load_datatype("https://doi.org/21.T11969/b9335ce2c99ed87735a6") # Group comparison

instance = dt1.data_analysis(
  is_implemented_by="data-analysis-1.py",
  has_part=dt2.group_comparison(
    label="Group comparison with petal length dependent variable for setosa and virginica irises",
    executes=dt2.software_method(
      label="ttest_ind",
      has_support_url="https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html",
      is_implemented_by="ttest_ind(x, y, equal_var = False)",
      part_of=dt2.software_library(
        label="scipy",
        has_support_url="https://docs.scipy.org/doc/scipy/index.html",
        version_info="1.15.1",
        part_of=dt2.software(
          label="Python",
          version_info="3.12.5",
          has_support_url="https://www.python.org"
        )
      )
    ),
    targets=dt2.component(label="petal length (cm)"),
    has_input=dt2.data_item(
      label="Iris data set with species and petal length",
      source_table=iris,
      has_part=[
        dt2.component(label="petal length (cm)"),
        dt2.component(label="species")
      ],
      has_characteristic=dt2.matrix_size(
        number_of_rows=iris.shape[0],
        number_of_columns=iris.shape[1]
      )
    ),
    has_output=dt2.data_item(
      source_table=df,
      has_part=[
        dt2.component(label="t.statistic"),
        dt2.component(label="df"),
        dt2.component(label="p.value")
      ],
      has_characteristic=dt2.matrix_size(
        number_of_rows=df.shape[0],
        number_of_columns=df.shape[1]
      ),
      has_expression=dt2.figure(source_url="data-analysis-1.png")
    )
  )
)

with open("data-analysis-1.json", "w") as f:
  f.write(to_jsonld(instance))
