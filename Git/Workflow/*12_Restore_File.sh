# Revert to a version from the specific commit
git log --since='Mar 3 2023' # Find the hash
git checkout first_6_8_hash file_name

# Later ...
git add file_name
git commit -m "Restoring version from commit first_6_8_hash"

# Revert to a version from the second-to-most recent commit
git checkout HEAD~1 file_name


