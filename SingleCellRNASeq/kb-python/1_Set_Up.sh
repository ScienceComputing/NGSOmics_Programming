work_dir = "/Users/your_name/SingleCell"
cd $work_dir/data/
path_list = ["/Users/your_name/SRATool/sratoolkit.3.0.5-mac64/bin"]
os.environ["PATH"] += os.pathsep + os.pathsep.join(path_list)
print(os.environ.get("PATH"))

pip install --quiet kb-python

# raise UnsupportedOSError
# Reference: https://github.com/ScienceComputing/Python_Programming/blob/main/Utilities/*Error_Type.md#runtime-errors
# 1. kb might have specific version requirements for Python
pyenv which python
pyenv install 3.8
pyenv local 3.8.18
pip install --quiet kb-python
# 2. kb tool is unable to find a compatible kallisto binary on your system
brew install cmake
git clone https://github.com/pachterlab/kallisto.git
cd kallisto
mkdir build
cd build
cmake ..
make
./src/kallisto version
cd ../..
kb ref -d human -i human_index.idx -g human_t2g.txt --verbose --kallisto /your_path/to/kallisto/build/src/kallisto
