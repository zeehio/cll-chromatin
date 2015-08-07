# clone and install pipelines
git clone git@github.com:afrendeiro/pipelines.git
cd pipelines
python setup.py install --user

# run pipelines
pipelines cll-patients ../metadata/sample_annotation.csv
