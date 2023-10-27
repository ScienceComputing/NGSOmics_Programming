# Add a single file to the staging area
git add file.py

# Add all files in the current directory to the staging area
git add .

# Unstage a single file
git reset HEAD file.py
# Modify this file
vim file.py
# Again add this file to the staging area
git add file.py
# Make the commit
git commit -m "XXX"

# Unstage all files in the current directory to the staging area
git reset HEAD
