# Undo changes to an unstaged (no git add) file; this will wipe out all changes made to the unstaged file forever
git checkout -- file_name

# Undo changes to all unstaged files in the current directory and subdirectories
git checkout .

# We move files from the staging area back into the repo and restore their state to versions in the last commit
git reset HEAD
git checkout . 
