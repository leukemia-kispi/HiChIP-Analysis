# Coolbox

## Setup of Coolbox and running it with JupyterNotebook

To install coolbox and use it with the Jupyter Notebook web application, clone the source code https://github.com/GangCaoLab/CoolBox.git. Enter the Coolbox directory and follow instruction for installation in conda environment.

```
git clone https://github.com/GangCaoLab/CoolBox.git
cd CoolBox
conda env create --file environment.yml
conda activate coolbox
python setup.py install
```

After installation, you should enable ipywidgets to use the browser in Jupyter notebook (depending on version this may be already enabled):

```
jupyter nbextension enable --py widgetsnbextension
```

Make sure you are in the coolbox conda environment and enter the directory of where you put the working directory. Called JypyterLab in example. (here you put all bigwig files or other data needed for generating visuals) 

```
cd directoryPATH/JupyterLab
jupyter lab --nobrowser --port 8585
```

If working on local machine when the server is initiated, a link with the token for access is made available. 
Copy the link and input it in a internet browser to begin using coolbox in JupyterNotebook.

If you run Jypyter from a virtual machine instance on a cloud you need to SSH Tunnel from Your Local Machine first
to later open it in the internet browser.

On your local machine, open a PowerShell or Command terminal and create an SSH tunnel:

```
ssh -L 8585:localhost:8585 username@remote_server_ip
```
Replace username with your remote server username.

Replace remote_server_ip with the actual IP address of your remote server.

Input http://localhost:8585 in web browser and respond to the password request. The password will corespond to instance password.

Commands and features to visualize genomic data with coolbox API can be found at https://gangcaolab.github.io/CoolBox/quick_start_API.html