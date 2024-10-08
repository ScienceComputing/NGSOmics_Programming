# To compare the files in a remote repo against the contents of a local repo
# Gather contents from the remote origin repo into the local main branch
git fetch origin main
# origin: name of the remote to fetch from
# main: branch of the local repo to fetch into

git fetch origin specific_local_branch_name

# [add] Compare the remote branch with its local branch
git diff <remote>/<remote branch> <local branch> 

# After fetching, we have the contents of the remote in our local repo
# We need to synchronize contents between the 2 repos
git merge origin main

# The remote repository often contains more up-to-date information than local repositories, making it a key source of truth. 
# Consequently, the process of fetching content from the remote and synchronizing it with the local repo is a common workflow.
# ! Fetch and merge together
git pull origin main

# When we're working on local changes that haven't been committed, Git prevents us from pulling from a remote.
# Git advises us to commit our changes, informing us that the pull command has been aborted. 
# This highlights the importance of committing our local work before pulling from a remote repo.
