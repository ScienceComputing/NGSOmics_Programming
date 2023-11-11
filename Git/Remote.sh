# A remote repo is a repo stored in the cloud through an online repo hosting service such as GitHub, Bitbucket, Gitlab

# Why use remote repos?
# 1. back up everything 2. facilitate the collaboration

# Copy the existing remote repo on GitHub to the local computer
git clone https://github.com/ScienceComputing/NGSOmics_Programming.git new_repo_name

# List the name of a remote repo if we are in that cloned directory
git remote
# Get the URLs for a remote repo
git remote -v

# Name the remote; by default, git name the remote with 'origin'
# Why define the remote name? 
# Use it as a shortcut when merging, instead of listing the URL or path
git remote add NGS https://github.com/ScienceComputing/NGSOmics_Programming.git new_repo_name
