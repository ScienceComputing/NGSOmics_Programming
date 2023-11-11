# Bring our local changes into a remote repo
# Save changes locally first!

# Push into remote from local branch
git push remote_name local_branch_name
git push origin main

# BUT: if there are different commits on the local and remote repos
git pull origin main
git pull --no-edit origin main # Avoid the message
git push origin main
