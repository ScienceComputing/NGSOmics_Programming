# Activate and deactivate a conda environment
conda activate atac
conda deactivate

# List all conda environments
conda info --envs
# conda environments:
#
#                         /Users/your_name/Library/r-miniconda-arm64
#                         /Users/your_name/Library/r-miniconda-arm64/envs/atac
# base                  *  /Users/your_name/miniconda3
# atac                     /Users/your_name/miniconda3/envs/atac
#                          /opt/homebrew/Caskroom/mambaforge/base/envs/pertpy-env

# Show all channels
conda config --show channels

# Add a channel
conda config --add channels 'channel name or website'

# Remove a channel
conda config --remove channels 'channel name or website'

# Install packages from pip to a conda environment
conda create -n atac
conda activate atac
conda install pip
/Users/your_name/miniconda3/envs/atac/bin/pip install macs2
