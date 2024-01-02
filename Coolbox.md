# Coolbox

To install coolbox clone the source code. Enter the Coolbox directory and follow instruction for installation.
```
git clone https://github.com/GangCaoLab/CoolBox.git
cd CoolBox
conda env create --file environment.yml
conda activate coolbox
python setup.py install
```
After installation, you should enable ipywidgets to use the browser in Jupyter notebook:

$ jupyter nbextension enable --py widgetsnbextension

cd /mnt/Jupyterlab
jupyter lab --nobrowser --port 8585
