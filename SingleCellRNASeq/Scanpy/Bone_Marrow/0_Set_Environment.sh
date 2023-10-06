ls /Users/your_name/.pyenv/versions/
# 3.11.6
pyenv install 3.11.6 
pyenv global 3.11.6 
# Notice that Scanpy cannot be installed on Python version 3.12.0; only versions >=3.8,<3.12 are supported;
# So we switch to a Pyenv-installed Python 3.11.6 globally
pip install scanpy
