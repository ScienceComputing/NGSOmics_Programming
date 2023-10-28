# Restore a repo to the most recent commit
git read HEAD
git checkout .
git add .
git commit -m 'Restore the repo to previous commit'

# Restore a repo to a specific commit
git checkout first_6_8_hash

# Restore a repo to the second-to-most recent commit
git checkout HEAD~1
