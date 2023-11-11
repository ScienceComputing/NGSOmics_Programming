# To compare the files in a remote repo against the contents of a local repo
git fetch origin main
# origin: name of the remote to fetch from
# main: branch of the local repo to fetch into

git fetch origin specific_local_branch_name

# After fetching, we have the contents of the remote in our local repo
# We need to synchronize contents between the 2 repos
git merge origin main

# The remote repository often contains more up-to-date information than local repositories, making it a key source of truth. 
# Consequently, the process of fetching content from the remote and synchronizing it with the local repo is a common workflow.
# ! Fetch and merge together
git pull origin main
