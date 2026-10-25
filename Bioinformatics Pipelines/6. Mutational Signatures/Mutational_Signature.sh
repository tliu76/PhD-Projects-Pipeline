conda create -n mutational-signatures python=3.11 pip -y
conda activate mutational-signatures
python -m pip install \
  "salamander-learn==0.3.2" \
  "numpy==1.24.4" \
  "scikit-learn==1.3.2" \
  jupyterlab ipykernel \
  matplotlib


python -m ipykernel install --user \
  --name mutational-signatures \
  --display-name "Python (mutational-signatures)"


conda env export --no-builds > environment.yml


prefix: /opt/homebrew/Caskroom/miniforge/base/envs/mutational-signatures
